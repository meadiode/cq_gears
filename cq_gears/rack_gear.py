#! /usr/bin/python3

'''
CQ_Gears - CadQuery based involute profile gear generator

Copyright 2021 meadiode@github

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

import numpy as np
import cadquery as cq

from .utils import circle3d_by3points, rotation_matrix, make_shell
from .spur_gear import GearBase


class RackGear(GearBase):

    def __init__(self, module, length, width, height,
                 pressure_angle=20.0, helix_angle=0.0, clearance=0.0,
                 backlash=0.0, **build_params):
        self.m = m = module
        self.a0 = a0 = np.radians(pressure_angle)
        self.clearance = clearance
        self.backlash = backlash
        self.helix_angle = np.radians(helix_angle)
        self.width = width
        self.length = length
        self.height = height
        self.build_params = build_params

        adn = self.ka * m # addendum
        ddn = self.kd * m # dedendum
        
        self.la = la = adn # addendum line
        self.ld = ld = -(ddn + clearance) # dedendum line
        
        # Tooth thickness on the pitch line
        s0 = m * (np.pi / 2.0 - backlash * np.tan(a0)) / 2.0 

        # Calculate rack involute profile points
        p1x = np.tan(a0) * abs(ld)
        p1p2 = (abs(la) + abs(ld)) / np.cos(a0)

        p1 = (-s0 - p1x, ld, 0.0)
        p2 = (np.sin(a0) * p1p2 + p1[0], np.cos(a0) * p1p2 + p1[1], 0.0)
        p3 = (-p2[0], p2[1], 0.0)
        p4 = (-p1[0], p1[1], 0.0)
        p5 = (p4[0] + (np.pi * m - p4[0] * 2.0), p4[1], 0.0)

        self.t_lflank_pts = np.array((p1, p2))
        self.t_tip_pts = np.array((p2, p3))
        self.t_rflank_pts = np.array((p3, p4))
        self.t_root_pts = np.array((p4, p5))

        self.tooth_height = abs(la) + abs(ld)

        # Number of teeth
        self.z = int(np.ceil(self.length / (np.pi * self.m)))



    def tooth_points(self):
        pts = np.concatenate((self.t_lflank_pts, self.t_tip_pts,
                              self.t_rflank_pts, self.t_root_pts))
        return pts


    def gear_points(self):
        tpts = np.concatenate((self.t_lflank_pts, self.t_tip_pts,
                               self.t_rflank_pts, self.t_root_pts))
        pts = tpts.copy()

        for i in range(10):
            ttpts = tpts.copy()
            ttpts[:, 0] += np.pi * self.m * (i + 1)
            pts = np.concatenate((pts, ttpts))

        return pts
    
    
    def _build_tooth_faces(self, helix_angle, x_pos, z_pos, width):
        
        tx = np.tan(helix_angle) * width
        t_faces = []

        for spline in (self.t_lflank_pts, self.t_tip_pts,
                       self.t_rflank_pts, self.t_root_pts):
            face_pts = []

            pts1, pts2 = spline.copy(), spline.copy()

            pts1[:, 0] += x_pos
            pts1[:, 2] += z_pos
            pts2[:, 0] += x_pos + tx
            pts2[:, 2] += z_pos + width

            face_pts.append([cq.Vector(*pt) for pt in pts1])
            face_pts.append([cq.Vector(*pt) for pt in pts2])

            face = cq.Face.makeSplineApprox(face_pts,
                                            tol=self.spline_approx_tol,
                                            minDeg=self.spline_approx_min_deg,
                                            maxDeg=self.spline_approx_max_deg)
            t_faces.append(face)

        return t_faces


    def _build_gear_faces(self):
        t_faces = self._build_tooth_faces(self.helix_angle, 0.0, 0.0,
                                          self.width)

        extra = int(abs(np.ceil(np.tan(self.helix_angle) * \
                                self.width / (np.pi * self.m))))
        
        cp_ext = 10.0
        lt_cut_plane = cq.Face.makePlane(length=self.tooth_height + cp_ext,
                                         width=self.width + cp_ext,
                                         basePnt=(0.0, 0.0, self.width / 2.0),
                                         dir=(-1.0, 0.0, 0.0))
        
        rt_cut_plane = cq.Face.makePlane(length=self.tooth_height + cp_ext,
                                         width=self.width + cp_ext,
                                         basePnt=(self.length, 0.0,
                                                  self.width / 2.0),
                                         dir=(1.0, 0.0, 0.0))
        
        if self.helix_angle != 0.0:  
            if self.helix_angle > 0.0:
                tidx = range(-extra, self.z + 1)
            else:
                tidx = range(self.z + extra + 1)
        else:
            tidx = range(self.z + 1)
        
        def get_xmin(face):
            bb = face.BoundingBox()
            return bb.xmin
        
        def get_xmax(face):
            bb = face.BoundingBox()
            return bb.xmax
        
        faces = []
        for i in tidx:
            for tf in t_faces:
                face = tf.translate((np.pi * self.m * i, 0.0, 0.0))
                
                if i <= extra + 1:
                    cpd = face.split(lt_cut_plane)
                    
                    if isinstance(cpd, cq.Compound):
                        face = max(list(cpd), key=get_xmax)
                    else:
                        face = cpd
                        if face.BoundingBox().xmax < 0.0:
                            continue
                
                if i >= self.z - extra - 1:
                    cpd = face.split(rt_cut_plane)
                    
                    if isinstance(cpd, cq.Compound):
                        face = min(list(cpd), key=get_xmin)
                    else:
                        face = cpd
                        if face.BoundingBox().xmin > self.length:
                            continue
                
                faces.append(face)

        # Build left side face
        wp = cq.Workplane('XY').add(faces)
        pt1 = wp.edges('<X').vertices('<Z').val()
        pt2 = wp.edges('<X').vertices('>Z').val()

        ls_wires = (cq.Workplane('YZ').add(wp)
                    .edges('<X')
                    .toPending()
                    .moveTo(pt1.Y, pt1.Z)
                    .lineTo(self.ld - self.height, pt1.Z)
                    .lineTo(self.ld - self.height, pt2.Z)
                    .lineTo(pt2.Y, pt2.Z)
                    .consolidateWires()).vals()

        ls_face = cq.Face.makeFromWires(ls_wires[0])
        faces.append(ls_face)

        # Build right side face
        pt1 = wp.edges('>X').vertices('<Z').val()
        pt2 = wp.edges('>X').vertices('>Z').val()

        rs_wires = (cq.Workplane('YZ', origin=(self.length, 0.0, 0.0)).add(wp)
                    .edges('>X')
                    .toPending()
                    .moveTo(pt1.Y, pt1.Z)
                    .lineTo(self.ld - self.height, pt1.Z)
                    .lineTo(self.ld - self.height, pt2.Z)
                    .lineTo(pt2.Y, pt2.Z)
                    .consolidateWires()).vals()

        rs_face = cq.Face.makeFromWires(rs_wires[0])
        faces.append(rs_face)

        # Build back face
        bk_wires = (cq.Workplane('XZ',
                                 origin=(0.0, self.ld - self.height, 0.0))
                    .rect(self.length, self.width, centered=False)).vals()

        bk_face = cq.Face.makeFromWires(bk_wires[0])
        faces.append(bk_face)

        # Build top face
        wp = wp.add(ls_face).add(rs_face).add(bk_face)
        tp_wires = (cq.Workplane('XY').add(wp)
                     .edges('>Z')
                     .toPending()
                     .consolidateWires()).vals()

        tp_face = cq.Face.makeFromWires(tp_wires[0])
        faces.append(tp_face)

        # Build bottom face
        bt_wires = (cq.Workplane('XY').add(wp)
                    .edges('<Z')
                    .toPending()
                    .consolidateWires()).vals()

        bt_face = cq.Face.makeFromWires(bt_wires[0])
        faces.append(bt_face)

        return faces


    def _build(self):
        faces = self._build_gear_faces()

        shell = make_shell(faces)
        body = cq.Solid.makeSolid(shell)

        return body


class HerringboneRackGear(RackGear):
    
    def _build_tooth_faces(self, helix_angle, x_pos, z_pos, width):
        tx = np.tan(helix_angle) * (width / 2.0)

        t_faces1 = (super(HerringboneRackGear, self)
                    ._build_tooth_faces(helix_angle, 0.0, 0.0, width / 2.0))

        t_faces2 = (super(HerringboneRackGear, self)
                    ._build_tooth_faces(-helix_angle, tx,
                                        width / 2.0, width / 2.0))

        return t_faces1 + t_faces2
