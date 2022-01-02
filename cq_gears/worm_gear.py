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

from .utils import rotation_matrix, make_shell, make_cross_section_face
from .spur_gear import GearBase


class Worm(GearBase):

    surface_splines = 8 # Number of curve splines to approximate a surface
    wire_comb_tol = 0.1 # Wire combining tolerance    
    t_face_parts = 2

    def __init__(self, module, lead_angle, n_threads, length,
                 pressure_angle=20.0, clearance=0.0,
                 backlash=0.0, **build_params):
        self.m = m = module
        self.a0 = a0 = np.radians(pressure_angle)
        self.clearance = clearance
        self.backlash = backlash
        self.lead_angle = np.radians(lead_angle)
        self.length = length
        self.n_threads = n_threads
        self.build_params = build_params

        # Pitch diameter
        d0 = self.n_threads * m / np.abs(np.tan(self.lead_angle))
        self.r0 = d0 / 2.0 # Pitch radius
        
        adn = self.ka * m # Addendum
        ddn = self.kd * m # Dedendum
        
        self.la = la = adn # Addendum line
        self.ld = ld = -(ddn + clearance) # Dedendum line
        
        self.ra = self.r0 + adn # Addendum radius
        self.rd = self.r0 - ddn # Dedendum radius

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
    
    
    def _build_tooth_faces(self):
        t_faces = []

        ttx = self.m * np.pi * self.n_threads
        start_x = -ttx / 2.0
        step_x = ttx / self.t_face_parts
        part_turn = np.pi * 2.0 / self.t_face_parts * np.sign(self.lead_angle)

        spline_tf = np.linspace((start_x, 0.0),
                                (start_x + step_x, part_turn),
                                self.surface_splines)

        for spline in (self.t_lflank_pts, self.t_tip_pts,
                       self.t_rflank_pts, self.t_root_pts):
            face_pts = []
            
            for tx, alpha in spline_tf:
                r_mat = rotation_matrix((1.0, 0.0, 0.0), alpha)
                pts = (spline + (tx, self.r0, 0.0)) @ r_mat

                face_pts.append([cq.Vector(*pt) for pt in pts])

            face = cq.Face.makeSplineApprox(face_pts,
                                            tol=self.spline_approx_tol,
                                            minDeg=self.spline_approx_min_deg,
                                            maxDeg=self.spline_approx_max_deg)
            t_faces.append(face)

        faces = []
        for n in range(self.t_face_parts):
            for tf in t_faces:
                faces.append(tf
                             .rotate((0.0, 0.0, 0.0),
                                     (1.0, 0.0, 0.0),
                                     np.degrees(-n * part_turn))
                             .translate((step_x * n, 0.0, 0.0)))

        return faces


    def _build_gear_faces(self):
        step = np.pi * self.m * self.n_threads
        turns = int(np.ceil(self.length / step)) + 2
        x_start = -turns * step / 2.0
        tau = np.pi * 2.0 / self.n_threads
        
        t_faces = self._build_tooth_faces()
        faces = []
        
        for th in range(self.n_threads):
            for tf in t_faces:
                faces.append(tf.rotate((0.0, 0.0, 0.0),
                                       (1.0, 0.0, 0.0),
                                       np.degrees(tau * th)))

        nfaces = []
        
        for i in range(turns):
            for tf in faces:
                nfaces.append(tf.translate((step / 2.0 + x_start + i * step,
                                           0.0, 0.0)))

        cp_size = self.ra * 2.0 + 2.0
        cp_x = self.length / 2.0
     
        left_cut_plane = cq.Face.makePlane(length=cp_size, width=cp_size,
                                           basePnt=(-cp_x, 0.0, 0.0),
                                           dir=(-1.0, 0.0, 0.0))
        right_cut_plane = cq.Face.makePlane(length=cp_size, width=cp_size,
                                            basePnt=(cp_x, 0.0, 0.0),
                                            dir=(1.0, 0.0, 0.0))

        lface = make_cross_section_face(nfaces, left_cut_plane,
                                        self.isection_tol, self.wire_comb_tol)
        rface = make_cross_section_face(nfaces, right_cut_plane,
                                        self.isection_tol, self.wire_comb_tol)

        def get_xmin(face):
            bb = face.BoundingBox()
            return bb.xmin

        def get_xmax(face):
            bb = face.BoundingBox()
            return bb.xmax    

        g_faces = []    

        for face in nfaces:
            bb = face.BoundingBox()
            
            if -(self.length / 2.0) < bb.xmin and bb.xmax < (self.length / 2.0):
                g_faces.append(face)
            else:
                cpd = face.split(left_cut_plane)
                if isinstance(cpd, cq.Compound):
                    g_faces.append(max(list(cpd), key=get_xmax))
                    
                cpd = face.split(right_cut_plane)
                if isinstance(cpd, cq.Compound):
                    g_faces.append(min(list(cpd), key=get_xmin))
                            
        g_faces.append(lface)
        g_faces.append(rface)

        return g_faces


    def _make_bore(self, body, bore_d):
        if bore_d is None:
            return body

        body = (cq.Workplane('YZ')
                .add(body)
                .faces('<X')
                .workplane()
                .circle(bore_d / 2.0)
                .cutThruAll())
        
        return body.val()


    def _build(self, bore_d=None):
        faces = self._build_gear_faces()

        shell = make_shell(faces, tol=self.shell_sewing_tol)
        body = cq.Solid.makeSolid(shell)

        body = self._make_bore(body, bore_d)

        return body

