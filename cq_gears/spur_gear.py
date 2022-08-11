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

from __future__ import annotations
import numpy as np
import cadquery as cq
from inspect import signature

from .utils import circle3d_by3points, rotation_matrix, make_shell
from .common import GearBase, PostProcMixin


class SpurGear(GearBase, PostProcMixin):

    def __init__(self, module, teeth_number, width,
                 pressure_angle=20.0, helix_angle=0.0, clearance=0.0,
                 backlash=0.0, **build_params):
        self.m = m = module
        self.z = z = teeth_number
        self.a0 = a0 = np.radians(pressure_angle)
        self.clearance = clearance
        self.backlash = backlash
        self.helix_angle = np.radians(helix_angle)
        self.width = width

        self.d0 = d0 = m * z     # pitch diameter
        self.adn = adn = self.ka / (z / d0) # addendum
        self.ddn = ddn = self.kd / (z / d0) # dedendum
        da = d0 + 2.0 * adn # addendum circle diameter
        dd = d0 - 2.0 * ddn - 2.0 * clearance # dedendum circle diameter
        self.s0 = m * (np.pi / 2.0 - backlash * np.tan(a0)) # tooth thickness on
                                                            # the pitch circle
        self.r0 = r0 = d0 / 2.0 # pitch radius
        self.ra = da / 2.0 # addendum radius
        self.rd = rd = dd / 2.0 # dedendum radius
        self.rb = rb = np.cos(a0) * d0 / 2.0 # base circle radius
        self.rr = max(rb, rd) # tooth root radius
        self.tau = tau = np.pi * 2.0 / z # pitch angle

        if helix_angle != 0.0:
            self.twist_angle = width / \
                                (r0 * np.tan(np.pi / 2.0 - self.helix_angle))
        else:
            self.surface_splines = 2
            self.twist_angle = 0.0

        self.build_params = build_params


    def tooth_profile(self):
        
        # Calculate involute curve points for the left side of the tooth
        r = np.linspace(self.rr, self.ra, self.curve_points)
        cos_a = self.r0 / r * np.cos(self.a0)
        
        a = np.arccos(np.clip(cos_a, -1.0, 1.0))
        inv_a0 = np.tan(self.a0) - self.a0
        inv_a = np.tan(a) - a
        s = r * (self.s0 / self.d0 + inv_a0 - inv_a)
        phi = s / r
        zeros = np.zeros(self.curve_points)

        lflank_pts = np.dstack((np.cos(phi) * r,
                                np.sin(phi) * r,
                                zeros)).squeeze()

        # Calculate tooth tip points - an arc lying on the addendum circle
        b = np.linspace(phi[-1], -phi[-1], self.curve_points)
        tip_pts = np.dstack((np.cos(b) * self.ra,
                             np.sin(b) * self.ra,
                             zeros)).squeeze()

        # Get right side involute curve points by mirroring the left side
        rflank_pts = np.dstack(((np.cos(-phi) * r)[::-1],
                                (np.sin(-phi) * r)[::-1],
                                zeros)).squeeze()

        # Calculate tooth root points - an arc that starts at the right side of
        # the tooth and goes to the left side of the next tooth. The mid-point
        # of that arc lies on the dedendum circle.
        rho = self.tau - phi[0] * 2.0
        # Get the three points defining the arc
        p1 = np.array((rflank_pts[-1][0],
                       rflank_pts[-1][1],
                       0.0))
        p2 = np.array((np.cos(-phi[0] - rho / 2.0) * self.rd,
                       np.sin(-phi[0] - rho / 2.0) * self.rd,
                       0.0))
        p3 = np.array((np.cos(-phi[0] - rho) * self.rr,
                       np.sin(-phi[0] - rho) * self.rr,
                       0.0))

        # Calculate arc's center and radius
        bcr, bcxy = circle3d_by3points(p1, p2, p3)

        # Calculate start and end angles
        t1 = np.arctan2(p1[1] - bcxy[1], p1[0] - bcxy[0])
        t2 = np.arctan2(p3[1] - bcxy[1], p3[0] - bcxy[0])

        if t1 < 0.0:
            t1 += np.pi * 2.0
        if t2 < 0.0:
            t2 += np.pi * 2.0

        t1, t2 = min(t1, t2), max(t1, t2)
        t = np.linspace(t1 + np.pi * 2.0, t2 + np.pi * 2.0, self.curve_points)

        root_pts = np.dstack((bcxy[0] + bcr * np.cos(t),
                              bcxy[1] + bcr * np.sin(t),
                              zeros)).squeeze()

        return (lflank_pts, tip_pts, rflank_pts, root_pts)

    
    def _build_tooth_faces(self, twist_angle_a, twist_angle_b, z_pos, width):
        surf_splines = int(np.ceil(abs(self.twist_angle) / np.pi))
        surf_splines = max(1, surf_splines) * self.surface_splines

        # Spline transformation parameters: (angle around z-axis, z-pos)
        spline_tf = np.linspace((twist_angle_a, z_pos),
                                (twist_angle_b, z_pos + width),
                                surf_splines)
        t_faces = []
        splines = self.tooth_profile()
        
        for spline in splines:
            face_pts = []

            for a, z in spline_tf:
                r_mat = rotation_matrix((0.0, 0.0, 1.0), a)

                pts = spline.copy()
                pts[:, 2] = z
                pts = pts @ r_mat

                face_pts.append([cq.Vector(*pt) for pt in pts])

            face = cq.Face.makeSplineApprox(face_pts,
                                            tol=self.spline_approx_tol,
                                            minDeg=self.spline_approx_min_deg,
                                            maxDeg=self.spline_approx_max_deg)
            t_faces.append(face)

        return t_faces


    def _build_gear_faces(self):
        t_faces = self._build_tooth_faces(0.0, self.twist_angle,
                                          0.0, self.width)
        faces = []
        for i in range(self.z):
            for tf in t_faces:
                faces.append(tf.rotate((0.0, 0.0, 0.0),
                                       (0.0, 0.0, 1.0),
                                       np.degrees(self.tau * i)))
        wp = cq.Workplane('XY').add(faces)
        
        topface_wires = cq.Wire.combine(wp.edges('<Z').vals(),
                                        tol=self.wire_comb_tol)
        topface = cq.Face.makeFromWires(topface_wires[0])
        botface_wires = cq.Wire.combine(wp.edges('>Z').vals(),
                                        tol=self.wire_comb_tol)
        botface = cq.Face.makeFromWires(botface_wires[0])
        wp = wp.add(topface).add(botface)

        return wp.vals()


    def _make_teeth_cutout_wire(self, plane, t1, t2, twist_angle):
        at1 = t1 * self.tau + self.tau / 2.0 + twist_angle
        at2 = t2 * self.tau + self.tau / 2.0 + twist_angle
        p1x = np.cos(at1)
        p1y = np.sin(at1)
        p2x = np.cos((at1 + at2) / 2.0)
        p2y = np.sin((at1 + at2) / 2.0)
        p3x = np.cos(at2)
        p3y = np.sin(at2)
        rc = self.ra + 1.0
        rd = self.rd - 0.01

        res = (plane
               .moveTo(p1x * rd, p1y * rd)
               .lineTo(p1x * rc, p1y * rc)
               .threePointArc((p2x * rc, p2y * rc), (p3x * rc, p3y * rc))
               .lineTo(p3x * rd, p3y * rd)
               .threePointArc((p2x * rd, p2y * rd),
                              (p1x * rd, p1y * rd))
               .close())

        return res


    def _remove_teeth(self, body, t1, t2):
        plane = (cq.Workplane('XY').workplane(offset=-0.1)
                 .add(body))

        if self.twist_angle == 0.0:
            cutout = (self._make_teeth_cutout_wire(plane, t1, t2, 0.0)
                      .extrude(self.width + 0.2, combine=False))
        else:
            cutout = (self._make_teeth_cutout_wire(plane, t1, t2, 0.0)
                      .twistExtrude(self.width + 0.2,
                                    np.degrees(-self.twist_angle),
                                    combine=False))

        body = (cq.Workplane('XY')
                .add(body)
                .cut(cutout)).val()

        return body


    def _make_missing_teeth(self, body, missing_teeth):
        if missing_teeth is None:
            return body

        if isinstance(missing_teeth[0], (list, tuple)):
            for t1, t2 in missing_teeth:
                body = self._remove_teeth(body, t1, t2)
        else:
            t1, t2 = missing_teeth
            body = self._remove_teeth(body, t1, t2)

        return body


    def _build(self, **kv_args):
        faces = self._build_gear_faces()
        shell = make_shell(faces, tol=self.shell_sewing_tol)
        body = cq.Solid.makeSolid(shell)

        body = self._post_process(body, **kv_args)

        return body


class HerringboneGear(SpurGear):


    def _build_tooth_faces(self, twist_angle_a, twist_angle_b, z_pos, width):
        t_faces1 = (super(HerringboneGear, self)
                    ._build_tooth_faces(0.0, self.twist_angle,
                                        0.0, self.width / 2.0))
        t_faces2 = (super(HerringboneGear, self)
                    ._build_tooth_faces(self.twist_angle, 0.0,
                                        self.width / 2.0, self.width / 2.0))
        return t_faces1 + t_faces2


    def _remove_teeth(self, body, t1, t2):
        plane = (cq.Workplane('XY').workplane(offset=-0.1))

        cutout = (self._make_teeth_cutout_wire(plane, t1, t2, 0.0)
                  .twistExtrude(self.width / 2.0 + 0.05,
                                np.degrees(-self.twist_angle))
                  .faces('>Z').workplane())

        cutout = (self._make_teeth_cutout_wire(cutout, t1, t2,
                                               -self.twist_angle)
                  .twistExtrude(self.width / 2.0 + 0.05,
                                np.degrees(self.twist_angle)))
        
        body = (cq.Workplane('XY')
                .add(body)
                .cut(cutout))

        return body
