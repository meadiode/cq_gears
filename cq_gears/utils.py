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
