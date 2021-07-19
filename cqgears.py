#! /usr/bin/python3

import numpy as np
import cadquery as cq

from OCP.Precision import Precision
from OCP.TColgp import TColgp_HArray2OfPnt, TColgp_HArray1OfPnt
from OCP.GeomAPI import (GeomAPI_PointsToBSplineSurface,
                         GeomAPI_PointsToBSpline)
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeFace,
                                BRepBuilderAPI_Sewing)


#
# Math utility functions
#

def sphere_to_cartesian(r, gamma, theta):
    '''Convert spherical coordinates to cartesian
       r - sphere radius
       gamma - polar angle
       theta - azimuth angle
       return - cartesian coordinates
    '''
    return (r * np.sin(gamma) * np.sin(theta),
            r * np.sin(gamma) * np.cos(theta),
            r * np.cos(gamma))


def s_arc(sr, c_gamma, c_theta, r_delta, start, end, n=32):
    '''Get arc points plotted on a sphere's surface
       sr - sphere radius
       c_gamma - polar angle of the circle's center
       c_theta - azimuth angle of the circle's center
       r_delta - angle between OC and OP, where O is sphere's center,
                 C is the circle's center, P is any point on the circle
       start - arc start angle
       end - arc end angle
       n - number of points
       return - circle points
    '''
    t = np.expand_dims(np.linspace(start, end, n), axis=1)
    a = sphere_to_cartesian(1.0, c_gamma + r_delta, c_theta)
    k = sphere_to_cartesian(1.0, c_gamma, c_theta)
    c = np.cos(t) * a + np.sin(t) * np.cross(k, a) + \
        np.dot(k, a) * (1.0 - np.cos(t)) * k
    c = c * sr
    return [dim.squeeze() for dim in np.hsplit(c, 3)]


def s_inv(gamma0, gamma):
    phi = np.arccos(np.tan(gamma0) / np.tan(gamma))    
    return np.arccos(np.cos(gamma) / np.cos(gamma0)) / np.sin(gamma0) - phi


def circle3d_by3points(a, b, c):
    u = b - a
    w = np.cross(c - a, u)
    u = u / np.linalg.norm(u)
    w = w / np.linalg.norm(w)
    v = np.cross(w, u)
    
    bx = np.dot(b - a, u)
    cx, cy = np.dot(c - a, u), np.dot(c - a, v)
    
    h = ((cx - bx / 2.0) ** 2 + cy ** 2 - (bx / 2.0) ** 2) / (2.0 * cy)
    cc = a + u * (bx / 2.0) + v * h
    r = np.linalg.norm(a - cc)

    return r, cc


def rotation_matrix(axis, alpha):
    ux, uy, uz = axis
    sina, cosa = np.sin(alpha), np.cos(alpha)
    r_mat = np.array((
        (cosa + (1.0 - cosa) * ux ** 2,
         ux * uy * (1.0 - cosa) - uz * sina,
         ux * uz * (1.0 - cosa) + uy * sina),
        (uy * ux * (1.0 - cosa) + uz * sina,
         cosa + (1.0 - cosa) * uy ** 2,
         uy * uz * (1.0 - cosa) - ux * sina),
        (uz * ux * (1.0 - cosa) - uy * sina,
         uz * uy * (1.0 - cosa) + ux * sina,
         cosa + (1.0 - cosa) * uz ** 2),
    ))

    return r_mat


def angle_between(o, a, b):
    p = a - o
    q = b - o
    return np.arccos(np.dot(p, q) / (np.linalg.norm(p) * np.linalg.norm(q)))


def project_to_xy_from_sphere_center(pts, sphere_r):
    gammas = np.arccos(np.dot(pts, np.array((0.0, 0.0, 1.0))) / \
                       np.linalg.norm(pts, axis=1))
    thetas = np.arctan2(pts[:, 0], pts[:, 1])
    radiuses = sphere_r /  np.cos(gammas)
    proj_pts = np.dstack(sphere_to_cartesian(radiuses,
                                             gammas, thetas)).squeeze()
    
    return proj_pts


def intersection_ray_xy(points_a, points_b):
    ''' Get intersection point between a ray defined by points a and b
        and the plane XY
    '''
    ab = points_b - points_a
    t = points_a[:, 2] / -ab[:, 2]
    res = points_a + ab * np.expand_dims(t, 1)
    return res


#
# OpenCascade utility functions
#


def make_spline_approx(points, tol=1e-2, smoothing=None, minDeg=1, maxDeg=3):
    '''
    Approximate a surface through the provided points.
    '''
    points_ = TColgp_HArray2OfPnt(1, len(points), 1, len(points[0]))

    for i, vi in enumerate(points):
        for j, v in enumerate(vi):
            v = cq.Vector(*v)
            points_.SetValue(i + 1, j + 1, v.toPnt())

    if smoothing:
        spline_builder = GeomAPI_PointsToBSplineSurface(
            points_, *smoothing, DegMax=maxDeg, Tol3D=tol
        )
    else:
        spline_builder = GeomAPI_PointsToBSplineSurface(
            points_, DegMin=minDeg, DegMax=maxDeg, Tol3D=tol
        )

    if not spline_builder.IsDone():
        raise ValueError("B-spline approximation failed")

    spline_geom = spline_builder.Surface()

    return cq.Face(BRepBuilderAPI_MakeFace(spline_geom,
                                           Precision.Confusion_s()).Face())


def make_spline_approx_1d(points, tol=1e-3, smoothing=None, minDeg=1, maxDeg=6):
    '''
    Approximate a spline through the provided points.
    '''
    points_ = TColgp_HArray1OfPnt(1, len(points))
    
    for ix, v in enumerate(points):
        v = cq.Vector(*v)
        points_.SetValue(ix + 1, v.toPnt())

    if smoothing:
        spline_builder = GeomAPI_PointsToBSpline(
                            points_, *smoothing, DegMax=maxDeg, Tol3D=tol
                         )
    else:
        spline_builder = GeomAPI_PointsToBSpline(
                            points_, DegMin=minDeg, DegMax=maxDeg, Tol3D=tol
                         )

    if not spline_builder.IsDone():
        raise ValueError("B-spline approximation failed")

    spline_geom = spline_builder.Curve()

    return cq.Edge(BRepBuilderAPI_MakeEdge(spline_geom).Edge())


def make_shell(faces):
    shell_builder = BRepBuilderAPI_Sewing(tolerance=1e-2)

    for face in faces:
        shell_builder.Add(face.wrapped)

    shell_builder.Perform()
    s = shell_builder.SewedShape()

    return cq.Shell(s)



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


    def build(self, bore_d=None, missing_teeth=None):
        faces = self._build_faces()

        shell = make_shell(faces)
        body = cq.Solid.makeSolid(shell)

        body = self._make_bore(body, bore_d)
        body = self._make_missing_teeth(body, missing_teeth)

        return body



class HerringBoneGear(SpurGear):

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
