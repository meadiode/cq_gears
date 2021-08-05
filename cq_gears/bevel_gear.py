#! /usr/bin/python3

import numpy as np
import cadquery as cq

from .utils import (circle3d_by3points, rotation_matrix, make_spline_approx,
                    s_inv, s_arc, angle_between, sphere_to_cartesian,
                    make_shell)

class BevelGear:

    ka = 1.0  # addendum coefficient
    kd = 1.25 # dedendum coefficient

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

        # Polar angles of base, face and root cones 
        self.gamma_b = gamma_b = np.arcsin(np.cos(a0) * np.sin(gamma_p))
        self.gamma_f = gamma_f = gamma_p + np.arctan(ka * m / gs_r)
        self.gamma_r = gamma_r = gamma_p - np.arctan(kd * m / gs_r)
                
        # Pitch angle
        self.tau = tau = np.pi * 2.0 / z

        # Tooth mirror point azimuth angle
        phi_r = s_inv(gamma_b, gamma_p);
        mp_theta = np.pi / z + 2.0 * phi_r

        # Tooth left flank curve points
        gamma_tr = max(gamma_b, gamma_r)
        gamma = np.linspace(gamma_tr, gamma_f, curve_points)
        theta = s_inv(gamma_b, gamma)
        self.tsidel = np.dstack(sphere_to_cartesian(gs_r,
                                                    gamma,
                                                    theta)).squeeze()
        # Tooth tip curve points
        theta_tip = np.linspace(theta[-1], mp_theta - theta[-1], curve_points)
        self.ttip = np.dstack(sphere_to_cartesian(gs_r,
                                                  np.full(curve_points,
                                                          gamma_f),
                                                  theta_tip)).squeeze()

        # Get the right flank curve points
        self.tsider = np.dstack(sphere_to_cartesian(gs_r,
                                                    gamma[::-1],
                                                    mp_theta - \
                                                        theta[::-1])).squeeze()
 
        # Tooth root curve points
        if gamma_r < gamma_b:
            p1 = self.tsider[-1]
            p2 = np.array(sphere_to_cartesian(gs_r, gamma_b, theta[0] + tau))
            p3 = np.array(sphere_to_cartesian(gs_r, gamma_r,
                                              (tau + mp_theta) / 2.0))
            rr, rcc = circle3d_by3points(p1, p2, p3)
            rcc_gamma = np.arccos(np.dot(p3, rcc) / \
                                  (np.linalg.norm(p3) * np.linalg.norm(rcc)))
            p1p3 = angle_between(rcc, p1, p3)
            a_start = (np.pi - p1p3 * 2.0) / 2.0
            a_end = -a_start + np.pi
            self.tooth_root = np.dstack(s_arc(gs_r, gamma_r + rcc_gamma,
                                              (tau + mp_theta) / 2.0,
                                              rcc_gamma,
                                              np.pi / 2.0 + a_start,
                                              np.pi / 2.0 + a_end,
                                              curve_points)).squeeze()
        else:
            r_theta = np.linspace(mp_theta - theta[0],
                                  theta[0] + tau, curve_points)
            self.tooth_root = np.dstack(sphere_to_cartesian(
                                                        gs_r,
                                                        np.full(curve_points,
                                                                gamma_tr),
                                                        r_theta)).squeeze()


    def tooth_points_xyz(self):
