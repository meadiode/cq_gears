#! /usr/bin/python3

import numpy as np
import cadquery as cq
import warnings

from .utils import (circle3d_by3points, rotation_matrix, make_spline_approx,
                    make_shell)

from .spur_gear import SpurGear, HerringboneGear



class RingGear(SpurGear):

    def __init__(self, module, teeth_number, width, rim_width,
                 pressure_angle=20.0, helix_angle=0.0, clearance=0.0,
                 backlash=0.0, curve_points=20, surface_splines=5):
        self.m = m = module
        self.z = z = teeth_number
        self.a0 = a0 = np.radians(pressure_angle)
        self.clearance = clearance
        self.backlash = backlash
        self.curve_points = curve_points
        self.helix_angle = np.radians(helix_angle)
        self.width = width
        self.rim_width = rim_width

        d0 = m * z         # pitch diameter
        adn = self.ka / (z / d0) # addendum
        ddn = self.kd / (z / d0) # dedendum
        da = d0 - 2.0 * adn # addendum circle diameter
        dd = d0 + 2.0 * ddn + 2.0 * clearance # dedendum circle diameter
        s0 = m * (np.pi / 2.0 + backlash * np.tan(a0)) # tooth thickness on
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

        self.rim_r = rd + rim_width

        # Calculate involute curve points for the left side of the tooth
        r = np.linspace(ra, rr, curve_points)
        cos_a = r0 / r * np.cos(a0)
        a = np.arccos(np.clip(cos_a, -1.0, 1.0))
        inv_a = np.tan(a) - a
        s = r * (s0 / d0 + inv_a0 - inv_a)
        phi = s / r
        self.tsidel_x = np.cos(phi) * r
        self.tsidel_y = np.sin(phi) * r


        # Calculate tooth tip points - an arc lying on the addendum circle
        b = np.linspace(phi[-1], -phi[-1], curve_points)
        self.ttip_x = np.cos(b) * rd
        self.ttip_y = np.sin(b) * rd


        # Get right side involute curve points by mirroring the left side
        self.tsider_x = (np.cos(-phi) * r)[::-1]
        self.tsider_y = (np.sin(-phi) * r)[::-1]


        # Calculate tooth root points - an arc that starts at the right side of
        # the tooth and goes to the left side of the next tooth. The mid-point
        # of that arc lies on the dedendum circle.
        rho = tau - phi[0] * 2.0
        # Get the three points defining the arc
        p1 = np.array((self.tsider_x[-1], self.tsider_y[-1], 0.0))
        p2 = np.array((np.cos(-phi[0] - rho / 2.0) * ra,
              np.sin(-phi[0] - rho / 2.0) * ra, 0.0))
        p3 = np.array((np.cos(-phi[0] - rho) * ra,
              np.sin(-phi[0] - rho) * ra, 0.0))

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
        t = np.linspace(t2 + np.pi * 2.0, t1 + np.pi * 2.0, curve_points)
        
        self.troot_x = bcxy[0] + bcr * np.cos(t)
        self.troot_y = bcxy[1] + bcr * np.sin(t)


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

        rim = cq.Wire.makeCircle(self.rim_r,
                                 cq.Vector(0.0, 0.0, 0.0),
                                 cq.Vector(0.0, 0.0, 1.0))
        
        wr = cq.Wire.combine(wires)
        face = cq.Face.makeFromWires(rim, wr)

        return face


    def _build_rim_face(self):
        w1 = cq.Wire.makeCircle(self.rim_r,
                                cq.Vector(0.0, 0.0, 0.0),
                                cq.Vector(0.0, 0.0, 1.0))
        w2 = cq.Wire.makeCircle(self.rim_r,
                                cq.Vector(0.0, 0.0, self.width),
                                cq.Vector(0.0, 0.0, 1.0))

        face = cq.Face.makeRuledSurface(w1, w2)

        return face


    def _build_faces(self):
        faces = super(RingGear, self)._build_faces()
        faces.append(self._build_rim_face())

        return faces


    def build(self):
        faces = self._build_faces()

        shell = make_shell(faces)
        body = cq.Solid.makeSolid(shell)

        return body


class HerringboneRingGear(RingGear):


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

        faces.append(self._build_rim_face())

        return faces


class PlanetaryGearset:

    gear_cls = SpurGear
    ring_gear_cls = RingGear

    def __init__(self, module, sun_teeth_number, planet_teeth_number, width,
                 rim_width, n_planets, pressure_angle=20.0,
                 helix_angle=0.0, clearance=0.0, backlash=0.0, curve_points=20,
                 surface_splines=5):

        ring_z = sun_teeth_number + planet_teeth_number * 2

        self.sun = self.gear_cls(module, sun_teeth_number, width,
                                 pressure_angle=pressure_angle,
                                 helix_angle=helix_angle,
                                 clearance=clearance,
                                 backlash=backlash,
                                 curve_points=curve_points,
                                 surface_splines=surface_splines)

        self.planet = self.gear_cls(module, planet_teeth_number, width,
                                    pressure_angle=pressure_angle,
                                    helix_angle=-helix_angle,
                                    clearance=clearance,
                                    backlash=backlash,
                                    curve_points=curve_points,
                                    surface_splines=surface_splines)

        self.ring = self.ring_gear_cls(module, ring_z, width, rim_width,
                                       pressure_angle=pressure_angle,
                                       helix_angle=-helix_angle,
                                       clearance=clearance,
                                       backlash=backlash,
                                       curve_points=curve_points,
                                       surface_splines=surface_splines)

        self.orbit_r = self.sun.r0 + self.planet.r0
        self.n_planets = n_planets

        if ((sun_teeth_number + planet_teeth_number) % n_planets):
           if (sun_teeth_number % n_planets) or (ring_z % n_planets):
            warnings.warn('Planet gears being spaced evenly probably won\'t '
                          'mesh properly (if at all) with the given number of '
                          'teeth for the sun/planet gears and the number '
                          'of planets')



    def build(self, sun=True, planets=True, ring=True,
              sun_build_args={}, planet_build_args={}, **kv_args):

        gearset = cq.Workplane('XY')

        if sun:
            args = {**kv_args, **sun_build_args}
            sun = self.sun.build(**args)

            if (self.planet.z % 2) != 0:
                sun = sun.rotate((0.0, 0.0, 0.0), (0.0, 0.0, 1.0),
                                   np.degrees(self.sun.tau / 2.0))
            gearset.add(sun)


        if planets and self.n_planets > 0:
            planet_a = np.pi * 2.0 / self.n_planets

            if isinstance(planets, (list, tuple)):
                tobuild = planets
            else:
                tobuild = [True, ] * self.n_planets
            args = {**kv_args, **planet_build_args}

            for i in range(self.n_planets):
                if not tobuild[i]:
                    continue

                planet = (self.planet.build(**args)
                          .rotate((0.0, 0.0, 0.0), (0.0, 0.0, 1.0),
                                  np.degrees(self.planet.tau / 2.0))
                          .translate((np.cos(i * planet_a) * self.orbit_r,
                                      np.sin(i * planet_a) * self.orbit_r,
                                      0.0)))
                gearset.add(planet)

        if ring:
            ring = (self.ring.build()
                    .rotate((0.0, 0.0, 0.0), (0.0, 0.0, 1.0),
                            np.degrees(self.ring.tau / 2.0)))
            gearset.add(ring)

        return gearset


class HerringbonePlanetaryGearset(PlanetaryGearset):
    gear_cls = HerringboneGear
    ring_gear_cls = HerringboneRingGear
