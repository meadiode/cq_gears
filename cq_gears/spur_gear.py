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


class GearBase:
    ka = 1.0  # Addendum coefficient
    kd = 1.25 # Dedendum coefficient

    curve_points = 20 # Number of points to approximate a curve
    surface_splines = 5 # Number of curve splines to approximate a surface
    
    wire_comb_tol = 1e-2 # Wire combining tolerance
    spline_approx_tol = 1e-2 # Surface spline approximation tolerance
    shell_sewing_tol = 1e-2 # Tolerance to assembly a shell out of faces
    isection_tol = 1e-7 # Tolerance to find intersections between two surfaces
    spline_approx_min_deg = 3 # Minimum surface spline degree
    spline_approx_max_deg = 8 # Maximum surface spline degree

    
    def __init__(self, *args, **kv_args):
        raise NotImplementedError('Constructor is not defined')

    
    def build(self, **kv_params):
        params = {**self.build_params, **kv_params}
        
        return self._build(**params)



class SpurGear(GearBase):

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

        d0 = m * z         # pitch diameter
        adn = self.ka / (z / d0) # addendum
        ddn = self.kd / (z / d0) # dedendum
        da = d0 + 2.0 * adn # addendum circle diameter
        dd = d0 - 2.0 * ddn - 2.0 * clearance # dedendum circle diameter
        s0 = m * (np.pi / 2.0 - backlash * np.tan(a0)) # tooth thickness on
                                                       # the pitch circle
        inv_a0 = np.tan(a0) - a0

        self.r0 = r0 = d0 / 2.0 # pitch radius
        self.ra = ra = da / 2.0 # addendum radius
        self.rd = rd = dd / 2.0 # dedendum radius
        self.rb = rb = np.cos(a0) * d0 / 2.0 # base circle radius
        self.rr = rr = max(rb, rd) # tooth root radius
        self.tau = tau = np.pi * 2.0 / z # pitch angle

        if helix_angle != 0.0:
            self.twist_angle = width / \
                                (r0 * np.tan(np.pi / 2.0 - self.helix_angle))
        else:
            self.surface_splines = 2
            self.twist_angle = 0.0

        self.build_params = build_params

        # Calculate involute curve points for the left side of the tooth
        r = np.linspace(rr, ra, self.curve_points)
        cos_a = r0 / r * np.cos(a0)
        a = np.arccos(np.clip(cos_a, -1.0, 1.0))
        inv_a = np.tan(a) - a
        s = r * (s0 / d0 + inv_a0 - inv_a)
        phi = s / r
        self.t_lflank_pts = np.dstack((np.cos(phi) * r,
                                       np.sin(phi) * r,
                                       np.zeros(self.curve_points))).squeeze()

        # Calculate tooth tip points - an arc lying on the addendum circle
        b = np.linspace(phi[-1], -phi[-1], self.curve_points)
        self.t_tip_pts = np.dstack((np.cos(b) * ra,
                                    np.sin(b) * ra,
                                    np.zeros(self.curve_points))).squeeze()

        # Get right side involute curve points by mirroring the left side
        self.t_rflank_pts = np.dstack(((np.cos(-phi) * r)[::-1],
                                       (np.sin(-phi) * r)[::-1],
                                       np.zeros(self.curve_points))).squeeze()

        # Calculate tooth root points - an arc that starts at the right side of
        # the tooth and goes to the left side of the next tooth. The mid-point
        # of that arc lies on the dedendum circle.
        rho = tau - phi[0] * 2.0
        # Get the three points defining the arc
        p1 = np.array((self.t_rflank_pts[-1][0],
                       self.t_rflank_pts[-1][1],
                       0.0))
        p2 = np.array((np.cos(-phi[0] - rho / 2.0) * rd,
              np.sin(-phi[0] - rho / 2.0) * rd, 0.0))
        p3 = np.array((np.cos(-phi[0] - rho) * rr,
              np.sin(-phi[0] - rho) * rr, 0.0))

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

        self.t_root_pts = np.dstack((bcxy[0] + bcr * np.cos(t),
                                     bcxy[1] + bcr * np.sin(t),
                                     np.zeros(self.curve_points))).squeeze()


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

    
    def _build_tooth_faces(self, twist_angle_a, twist_angle_b, z_pos, width):
        surf_splines = int(np.ceil(abs(self.twist_angle) / np.pi))
        surf_splines = max(1, surf_splines) * self.surface_splines

        # Spline transformation parameters: (angle around z-axis, z-pos)
        spline_tf = np.linspace((twist_angle_a, z_pos),
                                (twist_angle_b, z_pos + width),
                                surf_splines)
        t_faces = []

        for spline in (self.t_lflank_pts, self.t_tip_pts,
                       self.t_rflank_pts, self.t_root_pts):
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


    def _make_recess(self, body, hub_d, recess_d, recess):
        if recess is None:
            return body

        assert recess_d is not None, 'Recess diameter is not set'

        body = (cq.Workplane('XY')
                .add(body)
                .faces('>Z')
                .workplane())

        if hub_d is not None:
            body = body.circle(hub_d / 2.0)

        body = body.circle(recess_d / 2.0).cutBlind(-recess).val()

        return body


    def _make_hub(self, body, hub_d, hub_length, bore_d):
        if hub_length is None:
            return body

        assert hub_d is not None, 'Hub diameter is not set'

        body = (cq.Workplane('XY')
                .add(body)
                .faces('>Z')
                .workplane())

        if bore_d is not None:
            body = body.circle(bore_d / 2.0)

        body = body.circle(hub_d / 2.0).extrude(hub_length)

        return body.val()


    def _make_spokes(self, body, spokes_id, spokes_od, n_spokes,
                     spoke_width, spoke_fillet):
        if n_spokes is None:
            return body
        assert n_spokes > 1, 'Number of spokes must be > 1'
        assert spoke_width is not None, 'Spoke width is not set'
        assert spokes_od is not None, 'Outer spokes diameter is not set'

        if spokes_id is None:
            r1 = spoke_width / 2.0
        else:
            r1 = max(spoke_width / 2.0, spokes_id / 2.0)

        r2 = spokes_od / 2.0

        r1 += 0.0001
        r2 -= 0.0001

        tau = np.pi * 2.0 / n_spokes
        a1 = np.arcsin((spoke_width / 2.0) / (spokes_id / 2.0))
        a2 = np.arcsin((spoke_width / 2.0) / (spokes_od / 2.0))
        a3 = tau - a2
        a4 = tau - a1

        cutout = (cq.Workplane('XY').workplane(offset=-0.1)
                  .moveTo(np.cos(a1) * r1, np.sin(a1) * r1)
                  .lineTo(np.cos(a2) * r2, np.sin(a2) * r2)
                  .radiusArc((np.cos(a3) * r2, np.sin(a3) * r2), -r2)
                  .lineTo(np.cos(a4) * r1, np.sin(a4) * r1)
                  .radiusArc((np.cos(a1) * r1, np.sin(a1) * r1), r1)
                  .close()
                  .extrude(self.width + 1.0))

        if spoke_fillet is not None:
            cutout = cutout.edges('|Z').fillet(spoke_fillet)

        body = cq.Workplane('XY').add(body)

        for i in range(n_spokes):
            body = body.cut(cutout.rotate((0.0, 0.0, 0.0),
                                          (0.0, 0.0, 1.0),
                                          np.degrees(tau * i)))

        return body.val()


    def _make_chamfer(self, body, chamfer=None, chamfer_top=None,
                      chamfer_bottom=None):    
        E = 0.01
        
        if chamfer is None and chamfer_top is None and chamfer_bottom is None:
            return body
        
        if chamfer is not None:
            if chamfer_top is None:
                chamfer_top = chamfer
            if chamfer_bottom is None:
                chamfer_bottom = chamfer
                
        if chamfer_top is not None:
            if isinstance(chamfer_top, (list, tuple)):
                wx, wy = chamfer_top
            else:
                wx, wy = chamfer_top, chamfer_top
                
            cutter = (cq.Workplane('XZ')
                      .moveTo(self.ra - wx, self.width + E)
                      .hLine(wx + E)
                      .vLine(-wy - E)
                      .close()
                      .revolve())

            body = (cq.Workplane('XY')
                    .add(body)
                    .cut(cutter))
            
        if chamfer_bottom is not None:
            if isinstance(chamfer_bottom, (list, tuple)):
                wx, wy = chamfer_bottom
            else:
                wx, wy = chamfer_bottom, chamfer_bottom
                
            cutter = (cq.Workplane('XZ')
                      .moveTo(self.ra + E, wy)
                      .vLine(-wy - E)
                      .hLine(-wx - E)
                      .close()
                      .revolve())

            body = (cq.Workplane('XY')
                    .add(body)
                    .cut(cutter))
            
            
        return body.val()


    def _build(self, bore_d=None, missing_teeth=None,
               hub_d=None, hub_length=None, recess_d=None, recess=None,
               n_spokes=None, spoke_width=None, spoke_fillet=None,
               spokes_id=None, spokes_od=None, chamfer=None, chamfer_top=None,
               chamfer_bottom=None, *args, **kv_args):
            faces = self._build_gear_faces()

            shell = make_shell(faces, tol=self.shell_sewing_tol)
            body = cq.Solid.makeSolid(shell)

            body = self._make_chamfer(body, chamfer, chamfer_top, chamfer_bottom)
            body = self._make_bore(body, bore_d)
            body = self._make_missing_teeth(body, missing_teeth)
            body = self._make_recess(body, hub_d, recess_d, recess)
            body = self._make_hub(body, hub_d, hub_length, bore_d)
            
            if spokes_id is None:
                spokes_id = hub_d

            if spokes_od is None:
                spokes_od = recess_d

            body = self._make_spokes(body, spokes_id, spokes_od, n_spokes,
                                     spoke_width, spoke_fillet)


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
