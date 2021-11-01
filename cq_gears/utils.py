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

from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing
from OCP.GeomAPI import GeomAPI_IntSS
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeWire,
                                BRepBuilderAPI_MakeEdge,
                                BRepBuilderAPI_MakeFace)
from OCP.TopTools import TopTools_HSequenceOfShape, TopTools_ListOfShape
from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds
from OCP.ShapeFix import ShapeFix_Face
from OCP.TopoDS import TopoDS
from OCP.BRepCheck import BRepCheck_Analyzer


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
       c_gamma - polar angle of the arc's center
       c_theta - azimuth angle of the arc's center
       r_delta - angle between OC and OP, where O is sphere's center,
                 C is the arc's center, P is any point on the arc
       start - arc start angle
       end - arc end angle
       n - number of points
       return - arc points
    '''
    t = np.expand_dims(np.linspace(start, end, n), axis=1)
    a = sphere_to_cartesian(1.0, c_gamma + r_delta, c_theta)
    k = sphere_to_cartesian(1.0, c_gamma, c_theta)
    c = np.cos(t) * a + np.sin(t) * np.cross(k, a) + \
        np.dot(k, a) * (1.0 - np.cos(t)) * k
    c = c * sr
    return [dim.squeeze() for dim in np.hsplit(c, 3)]


def s_inv(gamma0, gamma):
    '''Spherical involute curve function
       gamma0 - base polar angle that corresponds to the pitch cone 
       gamma - polar angle(s) of point(s) lying on the curve, for which
               corresponding azimuth angle(s) to be calculated
       return - azimuth angle(s) of the point(s) lying on the curve

    ''' 
    phi = np.arccos(np.tan(gamma0) / np.tan(gamma))    
    return np.arccos(np.cos(gamma) / np.cos(gamma0)) / np.sin(gamma0) - phi


def circle3d_by3points(a, b, c):
    '''Find a circle in 3d space defined by 3 points
       a, b, c - 3d points a,b and c which are lying on the circle to find
       return - circle's radius, circle's center 
    '''
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
    '''Construct a 3d rotation transform matrix
       axis - a 3d-axis to rotate about
       alpha - an angle to rotate to
       return - an array of the shape (3, 3) representing resulting matrix
    '''
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
    '''Find an angle between two vectors - OA and OB
       o, a, b - 3d-points defining the vectors OA and OB
       return - an angle between those two vectors 
    '''
    p = a - o
    q = b - o
    return np.arccos(np.dot(p, q) / (np.linalg.norm(p) * np.linalg.norm(q)))


#
# OpenCascade utility functions
#

def make_shell(faces, tol=1e-2):
    '''
    This is similar to cq.Shell.makeShell, but allows to pass the tolerance
    parameter.
    '''
    shell_builder = BRepBuilderAPI_Sewing(tol)

    for face in faces:
        shell_builder.Add(face.wrapped)

    shell_builder.Perform()
    s = shell_builder.SewedShape()

    return cq.Shell(s)


def make_cross_section_face(faces, cut_plane, int_tol=1e-7, wire_con_tol=1e-3):
    ss = GeomAPI_IntSS()
    cps = BRepAdaptor_Surface(cut_plane.wrapped).Surface().Surface()

    curves = []

    for f in faces:
        gfs = BRepAdaptor_Surface(f.wrapped).Surface().Surface()
        ss.Perform(cps, gfs, int_tol)

        if ss.NbLines():
            for i in range(ss.NbLines()):
                curves.append(ss.Line(i + 1))

    edges = []

    for curve in curves:
        eb = BRepBuilderAPI_MakeEdge(curve)
        edges.append(eb.Edge())


    wb = BRepBuilderAPI_MakeWire()
    elist = TopTools_ListOfShape()

    for edge in edges:
        elist.Append(edge)
    
    wb.Add(elist)
    
    if not wb.IsDone():
        edges_in = TopTools_HSequenceOfShape()
        wires_out = TopTools_HSequenceOfShape()

        for edge in edges:
            edges_in.Append(edge)

        ShapeAnalysis_FreeBounds.ConnectEdgesToWires_s(edges_in, wire_con_tol,
                                                       False, wires_out)

        wire = TopoDS.Wire_s(wires_out.First())
    else:
        wire = wb.Wire()

    fb = BRepBuilderAPI_MakeFace(wire, True)
    face = fb.Face()

    if not cq.Face(face).isValid():
        # Erroneous face construction could be traced/analyzed here
        # an = BRepCheck_Analyzer(face)
        # anres = an.Result(face).Status().First()        
        
        fix = ShapeFix_Face(face)
        fix.FixOrientation()
        fix.Perform()

        face = fix.Face()

    return cq.Face(face)
