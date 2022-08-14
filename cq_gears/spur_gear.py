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
from .common import GearBase
from .postprocess import PostProcMixin


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


    def tooth_profile(self) -> tuple(np.ndarray):
        
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


    def tooth_trace_curve(self, t_start: float, t_end: float) -> np.ndarray:
        surf_splines = int(np.ceil(abs(self.twist_angle) / np.pi))
        surf_splines = max(1, surf_splines) * self.surface_splines

        return np.linspace((t_start * self.twist_angle, t_start * self.width),
                           (t_end * self.twist_angle, t_end * self.width),
                           surf_splines)


    def _extrude_faces(self, splines: tuple[np.ndarray]) -> list[cq.Face]:
        t_faces = []
        seg_start = 0.0
        
        for seg_end in self.tooth_trace_curve_segments:
            spline_tf = self.tooth_trace_curve(seg_start, seg_end)

            for spline in splines:
                face_pts = []

                for a, z in spline_tf:
                    r_mat = rotation_matrix((0.0, 0.0, 1.0), a)

                    pts = spline.copy()
                    pts[:, 2] = z
                    pts = pts @ r_mat

                    face_pts.append([cq.Vector(*pt) for pt in pts])

                face = cq.Face.makeSplineApprox(
                                            face_pts,
                                            tol=self.spline_approx_tol,
                                            minDeg=self.spline_approx_min_deg,
                                            maxDeg=self.spline_approx_max_deg)
                t_faces.append(face)
            
            seg_start = seg_end

        return t_faces        


    def _build_tooth_faces(self):
        return self._extrude_faces(self.tooth_profile()) 


    def _build_gear_faces(self):
        t_faces = self._build_tooth_faces()

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


    def _build(self, **kv_args):
        faces = self._build_gear_faces()
        shell = make_shell(faces, tol=self.shell_sewing_tol)
        body = cq.Solid.makeSolid(shell)

        body = self._post_process(body, **kv_args)

        return body


class HerringboneGear(SpurGear):
    
    tooth_trace_curve_segments = (0.5, 1.0)
    
    def tooth_trace_curve(self, t_start: float, t_end: float) -> np.ndarray:
        surf_splines = int(np.ceil(abs(self.twist_angle) / np.pi))
        surf_splines = max(1, surf_splines) * self.surface_splines
        
        if t_start == 0.0:
            return np.linspace((0.0, t_start * self.width),
                               (self.twist_angle, t_end * self.width),
                               surf_splines)
        else:
            return np.linspace((self.twist_angle, t_start * self.width),
                               (0.0, t_end * self.width),
                               surf_splines)

