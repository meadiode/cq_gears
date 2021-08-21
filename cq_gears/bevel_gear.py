#! /usr/bin/python3

import numpy as np
import cadquery as cq

from .utils import (circle3d_by3points, rotation_matrix, make_spline_approx,
                    s_inv, s_arc, angle_between, sphere_to_cartesian,
                    make_shell, project_to_xy_from_sphere_center,
                    intersection_ray_xy)

class BevelGear:

    ka = 1.0  # addendum coefficient
    kd = 1.25 # dedendum coefficient
    wire_comb_tol = 0.01 # wire combining tolerance

    def __init__(self, module, teeth_number, cone_angle, face_width,
                 pressure_angle=20.0, helix_angle=0.0, clearance=0.0,
                 backlash=0.0, curve_points=20, surface_splines=5):
        
        self.m = m = module
        self.z = z = teeth_number
        self.a0 = a0 = np.radians(pressure_angle)
        self.clearance = clearance
        self.backlash = backlash
        self.curve_points = curve_points
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
            self.surface_splines = surface_splines

            # Calculate twist(torsion) angle
            beta = np.arctan(face_width * np.tan(self.helix_angle) / \
                             (2.0 * gs_r - face_width))
            self.twist_angle = np.arcsin(gs_r / rp * np.sin(beta)) * 2.0
        else:
            self.surface_splines = 2
            self.twist_angle = 0.0

        # The distance between sphere's center and the bottom of the gear
        self.cone_h = np.cos(gamma_r) * gs_r

        # Tooth mirror point azimuth angle
        phi_r = s_inv(gamma_b, gamma_p);
        mp_theta = np.pi / z + 2.0 * phi_r

        # Tooth left flank curve points
        gamma_tr = max(gamma_b, gamma_r)
        gamma = np.linspace(gamma_tr, gamma_f, curve_points)
        theta = s_inv(gamma_b, gamma)
        self.tsidel = np.dstack(sphere_to_cartesian(1.0,
                                                    gamma,
                                                    theta)).squeeze()
        # Tooth tip curve points
        theta_tip = np.linspace(theta[-1], mp_theta - theta[-1], curve_points)
        self.ttip = np.dstack(sphere_to_cartesian(1.0,
                                                  np.full(curve_points,
                                                          gamma_f),
                                                  theta_tip)).squeeze()

        # Get the right flank curve points
        self.tsider = np.dstack(sphere_to_cartesian(1.0,
                                                    gamma[::-1],
                                                    mp_theta - \
                                                        theta[::-1])).squeeze()
 
        # Tooth root curve points
        if gamma_r < gamma_b:
            p1 = self.tsider[-1]
            p2 = np.array(sphere_to_cartesian(1.0, gamma_b, theta[0] + tau))
            p3 = np.array(sphere_to_cartesian(1.0, gamma_r,
                                              (tau + mp_theta) / 2.0))
            rr, rcc = circle3d_by3points(p1, p2, p3)
            rcc_gamma = np.arccos(np.dot(p3, rcc) / \
                                  (np.linalg.norm(p3) * np.linalg.norm(rcc)))
            p1p3 = angle_between(rcc, p1, p3)
            a_start = (np.pi - p1p3 * 2.0) / 2.0
            a_end = -a_start + np.pi
            self.troot = np.dstack(s_arc(1.0, gamma_r + rcc_gamma,
                                         (tau + mp_theta) / 2.0,
                                         rcc_gamma,
                                         np.pi / 2.0 + a_start,
                                         np.pi / 2.0 + a_end,
                                         curve_points)).squeeze()
        else:
            r_theta = np.linspace(mp_theta - theta[0],
                                  theta[0] + tau, curve_points)
            self.troot = np.dstack(sphere_to_cartesian(1.0,
                                                       np.full(curve_points,
                                                               gamma_tr),
                                                       r_theta)).squeeze()


    def tooth_points_xyz(self):
        pts = np.concatenate((self.tsidel, self.ttip,
                              self.tsider, self.troot))
        return pts

    
    def gear_points_xyz(self):
        tpts = np.concatenate((self.tsidel, self.ttip,
                               self.tsider, self.troot))
        pts = tpts.copy()
        angle = self.tau
        for i in range(self.z - 1):
            pts = np.concatenate((pts,
                                  tpts @ rotation_matrix((0.0, 0.0, 1.0),
                                                         angle)))
            angle += self.tau
        
        return pts


    def _build_tooth_profile(self):
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
        spline_tf = np.linspace((tc_f, ta1), (pc_f, ta2))

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
        for spline in (self.tsidel, self.ttip, self.tsider, self.troot):
            face_pts = []

            for r, a in spline_tf:
                r_mat = rotation_matrix((0.0, 0.0, 1.0), a)
                face_pts.append((spline @ r_mat) * r)

            # Make faces out of the points
            face = make_spline_approx(face_pts)

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


    def _build_faces(self):

        faces_wp = cq.Workplane('XY')

        t_faces = self._build_tooth_profile()

        for i in range(self.z):
            for tf in t_faces:
                faces_wp = faces_wp.add(tf.rotate((0.0, 0.0, 0.0),
                                                  (0.0, 0.0, 1.0),
                                                  np.degrees(self.tau * i)))

        topface_wires = cq.Wire.combine(faces_wp.edges('<Z').vals(),
                                        tol=self.wire_comb_tol)
        topface = cq.Face.makeFromWires(topface_wires[0])
        botface_wires = cq.Wire.combine(faces_wp.edges('>Z').vals(),
                                        tol=self.wire_comb_tol)
        botface = cq.Face.makeFromWires(botface_wires[0])

        faces_wp = faces_wp.add(topface).add(botface)

        return faces_wp.vals()


    def _trim_bottom(self, body):
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

        return body


    def _trim_top(self, body):
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

        return body


    def _make_bore(self, body, bore_d):
        if bore_d is None:
            return body

        body = (cq.Workplane('XY')
                .add(body)
                .faces('<Z')
                .workplane()
                .circle(bore_d / 2.0)
                .cutThruAll()
               ).val()
        
        return body


    def build(self, bore_d=None, trim_bottom=True, trim_top=True):
        faces = self._build_faces()

        shell = make_shell(faces)
        body = cq.Solid.makeSolid(shell)

        if trim_bottom:
            body = self._trim_bottom(body)

        if trim_top:
            body = self._trim_top(body)
        body = self._make_bore(body, bore_d)

        # Put the gear on its bottom
        body = (body
                .rotate((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 180.0)
                .translate((0.0, 0.0, self.cone_h)))

        return body