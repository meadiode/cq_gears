#! /usr/bin/python3

import numpy as np
import cadquery as cq

from .utils import (circle3d_by3points, rotation_matrix, make_spline_approx,
                    make_shell)


class SpurGear:

    def __init__(self, module, teeth_number, width,
                 pressure_angle=20.0, helix_angle=0.0, clearance=0.0,
                 curve_points=20, surface_splines=5):
        self.m = m = module
        self.z = z = teeth_number
        self.a0 = a0 = np.radians(pressure_angle)
        self.c = c = clearance
        self.curve_points = curve_points
        self.helix_angle = np.radians(helix_angle)
        self.width = width

        d0 = m * z         # pitch diameter
        da = m * (z + 2.0) # addendum circle diameter
        dd = m * (z - 2.0) - 2.0 * c # dedendum circle diameter
        s0 = m * (np.pi / 2.0 - c * np.tan(a0)) # tooth thickness on
                                                # the pitch circle
        inv_a0 = np.tan(a0) - a0

        self.r0 = r0 = d0 / 2.0 # pitch radius
        self.ra = ra = da / 2.0 # addendum radius
        self.rd = rd = dd / 2.0 # dedendum radius
        self.rb = rb = np.cos(a0) * d0 / 2.0 # base circle radius
        self.rr = rr = max(rb, rd) # tooth root radius
        self.tau = tau = np.pi * 2.0 / z # pitch angle

        if helix_angle != 0.0:
            self.surface_splines = surface_splines
            self.twist_angle = width / \
                                (r0 * np.tan(np.pi / 2.0 - self.helix_angle))
        else:
            self.surface_splines = 2
            self.twist_angle = 0.0

        # Calculate involute curve points for the left side of the tooth
        r = np.linspace(rr, ra, curve_points)
        cos_a = r0 / r * np.cos(a0)
        a = np.arccos(np.clip(cos_a, -1.0, 1.0))
        inv_a = np.tan(a) - a
        s = r * (s0 / d0 + inv_a0 - inv_a)
        phi = s / r
        self.tsidel_x = np.cos(phi) * r
        self.tsidel_y = np.sin(phi) * r


        # Calculate tooth tip points - an arc lying on the addendum circle
        b = np.linspace(phi[-1], -phi[-1], curve_points)
        self.ttip_x = np.cos(b) * ra
        self.ttip_y = np.sin(b) * ra


        # Get right side involute curve points by mirroring the left side
        self.tsider_x = (np.cos(-phi) * r)[::-1]
        self.tsider_y = (np.sin(-phi) * r)[::-1]


        # Calculate tooth root points - an arc starting at the right side of
        # the tooth and goes to the left side of the next tooth. The mid-point
        # of that arc lies on the dedendum circle.
        rho = tau - phi[0] * 2.0
        # Get the three points defining the arc
        p1 = np.array((self.tsider_x[-1], self.tsider_y[-1], 0.0))
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
        t = np.linspace(t1 + np.pi * 2.0, t2 + np.pi * 2.0, curve_points)
        
        self.troot_x = bcxy[0] + bcr * np.cos(t)
        self.troot_y = bcxy[1] + bcr * np.sin(t)
    

    def tooth_points_xy(self):
        pts_x = np.concatenate((self.tsidel_x, self.ttip_x,
                                self.tsider_x, self.troot_x))
        pts_y = np.concatenate((self.tsidel_y, self.ttip_y,
                                self.tsider_y, self.troot_y))
        return pts_x, pts_y


    def gear_points_xy(self):
        tx, ty = self.tooth_points_xy()
        pts_x, pts_y = tx.copy(), ty.copy()

        angle = -self.tau
        for i in range(self.z - 1):
            pts_x = np.concatenate((pts_x,
                                    tx * np.cos(angle) - ty * np.sin(angle)))
            pts_y = np.concatenate((pts_y,
                                    tx * np.sin(angle) + ty * np.cos(angle)))
            # Going clock-wise
            angle -= self.tau

        return pts_x, pts_y


    def _build_horizontal_face(self):
        tooth_splines = []

        for x, y in zip((self.tsidel_x, self.ttip_x,
                         self.tsider_x, self.troot_x),
                        (self.tsidel_y, self.ttip_y,
                         self.tsider_y, self.troot_y)):
            pts = []

            for pt in np.dstack((x, y)).squeeze():
                pts.append(cq.Vector(tuple(pt)))

            tooth_splines.append(cq.Edge.makeSpline(pts))

        tooth_wire = cq.Wire.assembleEdges(tooth_splines)

        wires = []
        for i in range(self.z): 
            wires.append(tooth_wire.rotate((0.0, 0.0, 0.0),
                                           (0.0, 0.0, 1.0),
                                           np.degrees(self.tau) * i))

        wr = cq.Wire.combine(wires)
        face = cq.Face.makeFromWires(wr[0])

        return face


    def _build_profile(self, twist_angle_a, twist_angle_b, z_pos, width):
        splines = []

        for x, y in zip((self.tsidel_x, self.ttip_x,
                         self.tsider_x, self.troot_x),
                        (self.tsidel_y, self.ttip_y,
                         self.tsider_y, self.troot_y)):
            pts = np.dstack((x, y, np.full(x.shape[0], 0.0))).squeeze()
            splines.append(pts)

        # Spline transformation parameters: (angle around z-axis, z-pos)
        spline_tf = np.linspace((twist_angle_a, z_pos),
                                (twist_angle_b, z_pos + width),
                                self.surface_splines)
        t_faces = []

        for spline in splines:
            face_pts = []

            for a, z in spline_tf:
                r_mat = rotation_matrix((0.0, 0.0, 1.0), a)

                pts = spline.copy()
                pts[:, 2] = z
                pts = pts @ r_mat

                face_pts.append(pts)

            face_pts = np.stack(face_pts)
            face = make_spline_approx(face_pts, maxDeg=15)
            t_faces.append(face)


        faces = []
        for i in range(self.z):
            for tf in t_faces:
                faces.append(tf.rotate((0.0, 0.0, 0.0),
                                       (0.0, 0.0, 1.0),
                                       np.degrees(self.tau * i)))

        return faces


    def _build_faces(self):
        faces = self._build_profile(0.0, self.twist_angle, 0.0, self.width)

        bottom = self._build_horizontal_face()
        top = bottom.rotate((0.0, 0.0, 0.0),
                            (0.0, 0.0, 1.0),
                            -np.degrees(self.twist_angle))
        top = top.translate((0.0, 0.0, self.width))
        faces.append(bottom)
        faces.append(top)

        return faces


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
                .cut(cutout))

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

        body = body.circle(recess_d / 2.0).cutBlind(-recess)

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

        return body


    def _make_spokes(self, body, hub_d, recess_d, n_spokes,
                     spoke_width, spoke_fillet):
        if n_spokes is None:
            return body
        assert n_spokes > 0, 'Number of spokes should be > 0'
        assert spoke_width is not None, 'Spoke width is not set'
        assert recess_d is not None, 'Recess diameter is not set'

        if hub_d is None:
            r1 = spoke_width / 2.0
        else:
            r1 = max(spoke_width / 2.0, hub_d / 2.0)

        r2 = recess_d / 2.0

        r1 += 0.0001
        r2 -= 0.0001

        tau = np.pi * 2.0 / n_spokes
        a1 = np.arcsin((spoke_width / 2.0) / (hub_d / 2.0))
        a2 = np.arcsin((spoke_width / 2.0) / (recess_d / 2.0))
        a3 = tau - a2
        a4 = tau - a1

        cutout = (cq.Workplane('XY').workplane(offset=-0.1)
                  .moveTo(np.cos(a1) * r1, np.sin(a1) * r1)
                  .lineTo(np.cos(a2) * r2, np.sin(a2) * r2)
                  .radiusArc((np.cos(a3) * r2, np.sin(a3) * r2), -r2)
                  .lineTo(np.cos(a4) * r1, np.sin(a4) * r1)
                  .radiusArc((np.cos(a1) * r1, np.sin(a1) * r1), r1)
                  .close()
                  .extrude(self.width + 1.0)
                  )

        if spoke_fillet is not None:
            cutout = cutout.edges('|Z').fillet(spoke_fillet)

        for i in range(n_spokes):
            body = body.cut(cutout.rotate((0.0, 0.0, 0.0),
                                          (0.0, 0.0, 1.0),
                                          np.degrees(tau * i)))

        return body


    def build(self, bore_d=None, missing_teeth=None,
              hub_d=None, hub_length=None, recess_d=None, recess=None,
              n_spokes=None, spoke_width=None, spoke_fillet=None,
              spokes_id=None, spokes_od=None):
        faces = self._build_faces()

        shell = make_shell(faces)
        body = cq.Solid.makeSolid(shell)

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

    def _build_faces(self):
        faces1 = self._build_profile(0.0, self.twist_angle,
                                     0.0, self.width / 2.0)
        faces2 = self._build_profile(self.twist_angle, 0.0,
                                     self.width / 2.0, self.width / 2.0)
        faces = faces1 + faces2

        bottom = self._build_horizontal_face()
        top = bottom.translate((0.0, 0.0, self.width))
        faces.append(bottom)
        faces.append(top)

        return faces


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
