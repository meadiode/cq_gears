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

from .utils import (circle3d_by3points, rotation_matrix, s_inv, s_arc,
                    angle_between, sphere_to_cartesian, make_shell)

from .spur_gear import GearBase


class BevelGear(GearBase):

    def __init__(self, module, teeth_number, cone_angle, face_width,
                 pressure_angle=20.0, helix_angle=0.0, clearance=0.0,
                 backlash=0.0, **build_params):
        
        self.m = m = module
        self.z = z = teeth_number
        self.a0 = a0 = np.radians(pressure_angle)
        self.clearance = clearance
        self.backlash = backlash
        self.helix_angle = np.radians(helix_angle)
        self.face_width = face_width
        
        # Pitch cone angle
        self.gamma_p = gamma_p = np.radians(cone_angle)

        # base/pitch circle radius
        rp = m * z / 2.0
        
        # great sphere radius, also corresponds to pitch cone flank length
        self.gs_r = gs_r = rp / np.sin(gamma_p)

        assert gs_r > face_width, ('face_width value is too big, '
                                   'it should be < {:0.3f}'.format(gs_r))

        # Polar angles of base, face and root cones 
        self.gamma_b = gamma_b = np.arcsin(np.cos(a0) * np.sin(gamma_p))
        self.gamma_f = gamma_f = gamma_p + np.arctan(self.ka * m / gs_r)
        self.gamma_r = gamma_r = gamma_p - np.arctan(self.kd * m / gs_r)
                
        # Pitch angle
        self.tau = tau = np.pi * 2.0 / z

        if helix_angle != 0.0:
            # Calculate twist(torsion) angle
            beta = np.arctan(face_width * np.tan(self.helix_angle) / \
                             (2.0 * gs_r - face_width))
            self.twist_angle = np.arcsin(gs_r / rp * np.sin(beta)) * 2.0
        else:
            self.surface_splines = 2
            self.twist_angle = 0.0

        self.build_params = build_params

        # The distance between the cone apex and the bottom of the gear
        self.cone_h = np.cos(gamma_r) * gs_r

        # Tooth mirror point azimuth angle
        phi_r = s_inv(gamma_b, gamma_p);
        self.mp_theta = mp_theta = np.pi / z + 2.0 * phi_r

        # Tooth left flank curve points
        gamma_tr = max(gamma_b, gamma_r)
        gamma = np.linspace(gamma_tr, gamma_f, self.curve_points)
        theta = s_inv(gamma_b, gamma)
        self.t_lflank_pts = np.dstack(sphere_to_cartesian(1.0,
                                                          gamma,
                                                          theta)).squeeze()
        # Tooth tip curve points
        theta_tip = np.linspace(theta[-1],
                                mp_theta - theta[-1],
                                self.curve_points)
        self.t_tip_pts = np.dstack(
                            sphere_to_cartesian(1.0,
                                                np.full(self.curve_points,
                                                        gamma_f),
                                                theta_tip)).squeeze()

        # Get the right flank curve points
        self.t_rflank_pts = np.dstack(
                        sphere_to_cartesian(1.0,
                                            gamma[::-1],
                                            mp_theta - theta[::-1])).squeeze()
 
        # Tooth root curve points
        if gamma_r < gamma_b:
            p1 = self.t_rflank_pts[-1]
            p2 = np.array(sphere_to_cartesian(1.0, gamma_b, theta[0] + tau))
            p3 = np.array(sphere_to_cartesian(1.0, gamma_r,
                                              (tau + mp_theta) / 2.0))
            rr, rcc = circle3d_by3points(p1, p2, p3)
            rcc_gamma = np.arccos(np.dot(p3, rcc) / \
                                  (np.linalg.norm(p3) * np.linalg.norm(rcc)))
            p1p3 = angle_between(rcc, p1, p3)
            a_start = (np.pi - p1p3 * 2.0) / 2.0
            a_end = -a_start + np.pi
            self.t_root_pts = np.dstack(s_arc(1.0, gamma_r + rcc_gamma,
                                              (tau + mp_theta) / 2.0,
                                              rcc_gamma,
                                              np.pi / 2.0 + a_start,
                                              np.pi / 2.0 + a_end,
                                              self.curve_points)).squeeze()
        else:
            r_theta = np.linspace(mp_theta - theta[0],
                                  theta[0] + tau, self.curve_points)
            self.t_root_pts = np.dstack(
                                sphere_to_cartesian(1.0,
                                                    np.full(self.curve_points,
                                                            gamma_tr),
                                                    r_theta)).squeeze()


    def tooth_points(self):
        pts = np.concatenate((self.t_lflank_pts, self.t_tip_pts,
                              self.t_rflank_pts, self.t_root_pts))
        return pts

    
    def gear_points(self):
        tpts = np.concatenate((self.t_lflank_pts, self.t_tip_pts,
                               self.t_rflank_pts, self.t_root_pts))
        pts = tpts.copy()
        angle = self.tau
        for i in range(self.z - 1):
            pts = np.concatenate((pts,
                                  tpts @ rotation_matrix((0.0, 0.0, 1.0),
                                                         angle)))
            angle += self.tau
        
        return pts


    def _build_tooth_faces(self):
        pc_h = np.cos(self.gamma_r) * self.gs_r # pitch cone height
        pc_f = pc_h / np.cos(self.gamma_f) # extended pitch cone flank length
        pc_rb = pc_f * np.sin(self.gamma_f) # pitch cone base radius

        # top cone height
        tc_h = np.cos(self.gamma_f) * (self.gs_r - self.face_width)
        tc_f = tc_h / np.cos(self.gamma_r) # top cone flank length
        tc_rb = tc_f * np.sin(self.gamma_f) # top cone base radius

        # start/stop twist angle
        ta1 = -(pc_f - self.gs_r) / self.face_width * self.twist_angle
        ta2 = (self.gs_r - tc_f) / self.face_width * self.twist_angle

        # Transformation parameters: (radius, twist angle)
        spline_tf = np.linspace((pc_f, ta1), (tc_f - 0.01, ta2))

        tcp_size = tc_rb * 1000.0
        top_cut_plane = cq.Face.makePlane(length=tcp_size, width=tcp_size,
                                          basePnt=(0.0, 0.0, tc_h),
                                          dir=(0.0, 0.0, 1.0))
        bcp_size = pc_rb * 1000.0
        bott_cut_plane = cq.Face.makePlane(length=bcp_size, width=bcp_size,
                                           basePnt=(0.0, 0.0, pc_h),
                                           dir=(0.0, 0.0, 1.0))


        def get_zmax(face):
            bb = face.BoundingBox()
            return bb.zmax

        t_faces = []
        for spline in (self.t_lflank_pts, self.t_tip_pts,
                       self.t_rflank_pts, self.t_root_pts):
            face_pts = []

            for r, a in spline_tf:
                r_mat = rotation_matrix((0.0, 0.0, 1.0), a)
                pts = (spline @ r_mat) * r
                face_pts.append([cq.Vector(*pt) for pt in pts])

            # Make faces out of the points
            face = cq.Face.makeSplineApprox(face_pts,
                                            tol=self.spline_approx_tol,
                                            minDeg=self.spline_approx_min_deg,
                                            maxDeg=self.spline_approx_max_deg)

            # Trim top
            cpd = face.split(top_cut_plane)
            if isinstance(cpd, cq.Compound):
                face = max(list(cpd), key=get_zmax)
            else:
                face = cpd

            # Trim bottom
            cpd = face.split(bott_cut_plane)
            if isinstance(cpd, cq.Compound):
                face = min(list(cpd), key=get_zmax)
            else:
                face = cpd

            t_faces.append(face)

        return t_faces


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


    def _trim_bottom(self, body, do_trim=False):
        if not do_trim:
            return body
        r = self.gs_r
        
        p1 = sphere_to_cartesian(r, self.gamma_r * 0.99, np.pi / 2.0)
        p2 = sphere_to_cartesian(r, self.gamma_p, np.pi / 2.0)
        p3 = sphere_to_cartesian(r, self.gamma_f * 1.01, np.pi / 2.0)

        x1 = np.tan(self.gamma_f) * self.cone_h + 1.0

        trimmer = (cq.Workplane('XZ')
                   .moveTo(p1[0], p1[2])
                   .threePointArc((p2[0], p2[2]), (p3[0], p3[2]))
                   .lineTo(x1, p3[2])
                   .lineTo(x1, p1[2])
                   .close()
                   .revolve(combine=False))

        body = cq.Workplane('XY').add(body).cut(trimmer)

        return body.val()


    def _trim_top(self, body, do_trim=False):
        if not do_trim:
            return body
        r = self.gs_r - self.face_width

        p1 = sphere_to_cartesian(r, self.gamma_r, np.pi / 2.0)
        p2 = sphere_to_cartesian(r, self.gamma_p, np.pi / 2.0)
        p3 = sphere_to_cartesian(r, self.gamma_f * 1.01, np.pi / 2.0)

        trimmer = (cq.Workplane('XZ')
                   .moveTo(p1[0], p1[2])
                   .threePointArc((p2[0], p2[2]), (p3[0], p3[2]))
                   .lineTo(0.0, p3[2])
                   .lineTo(0.0, p1[2])
                   .close()
                   .revolve(combine=False))

        body = cq.Workplane('XY').add(body).cut(trimmer)

        return body.val()


    def _make_bore(self, body, bore_d):
        if bore_d is None:
            return body

        body = (cq.Workplane('XY')
                .add(body)
                .faces('<Z')
                .workplane()
                .circle(bore_d / 2.0)
                .cutThruAll())
        
        return body.val()


    def _build(self, bore_d=None, trim_bottom=True, trim_top=True):
        faces = self._build_gear_faces()

        shell = make_shell(faces)
        body = cq.Solid.makeSolid(shell)

        body = self._trim_bottom(body, trim_bottom)
        body = self._trim_top(body, trim_top)

        t_align_angle = -self.mp_theta / 2.0 - np.pi / 2.0 + np.pi / self.z

        # Put the gear on its bottom and align one of the teeth to x axis
        body = (cq.Workplane('XY')
                .add(body)
                .rotate((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 180.0)
                .translate((0.0, 0.0, self.cone_h))
                .rotate((0.0, 0.0, 0.0), (0.0, 0.0, 1.0),
                        np.degrees(t_align_angle))).solids().val()

        body = self._make_bore(body, bore_d)

        return body


class BevelGearPair(GearBase):

    gear_cls = BevelGear
    
    def __init__(self, module, gear_teeth, pinion_teeth, face_width,
                 axis_angle=90.0, pressure_angle=20.0, helix_angle=0.0,
                 clearance=0.0, **build_params):
        
        self.axis_angle = axis_angle = np.radians(axis_angle)
        
        aa_sin = np.sin(axis_angle)
        aa_cos = np.cos(axis_angle)

        # Cone Angle of the Gear
        delta_gear = np.arctan(aa_sin / (pinion_teeth / gear_teeth + aa_cos))
        
        # Cone Angle of the Pinion
        delta_pinion = np.arctan(aa_sin / (gear_teeth / pinion_teeth + aa_cos))
        
        self.gear = self.gear_cls(module, gear_teeth, np.degrees(delta_gear),
                                  face_width, pressure_angle, helix_angle,
                                  clearance)
        
        self.pinion = self.gear_cls(module, pinion_teeth,
                                    np.degrees(delta_pinion), face_width,
                                    pressure_angle, -helix_angle)
        self.build_params = build_params


    def _build(self, gear=True, pinion=True, transform_pinion=True,
               gear_build_args={}, pinion_build_args={}, **kv_args):

        gearset = cq.Workplane('XY')

        if gear:
            args = {**kv_args, **gear_build_args}
            gear = self.gear.build(**args)

            gearset.add(gear)

        if pinion:
            args = {**kv_args, **pinion_build_args}
            pinion = self.pinion.build(**args)

            if transform_pinion:
                dist = -self.pinion.cone_h + self.gear.cone_h
                angle = np.pi / 2.0 - self.axis_angle

                if self.pinion.z % 2 == 0:
                    pinion = pinion.rotate(
                                    (0.0, 0.0, 0.0),
                                    (0.0, 0.0, 1.0),
                                    np.degrees(np.pi / self.pinion.z))

                pinion = (pinion
                          .rotate((0.0, -1.0, self.gear.cone_h),
                                  (0.0, 1.0, self.gear.cone_h),
                                   np.degrees(self.axis_angle))
                          .translate((dist * np.cos(angle),
                                      0.0,
                                      dist * np.sin(angle))))

            gearset.add(pinion)

        gearset = gearset.vals()
        
        if len(gearset) == 1:
            return gearset[0]

        return cq.Compound.makeCompound(gearset)
