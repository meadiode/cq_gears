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
import warnings

from .utils import circle3d_by3points, rotation_matrix, make_shell
from .spur_gear import GearBase, SpurGear, HerringboneGear


class RingGear(SpurGear):

    def __init__(self, module, teeth_number, width, rim_width,
                 pressure_angle=20.0, helix_angle=0.0, clearance=0.0,
                 backlash=0.0, **build_params):
        self.m = m = module
        self.z = z = teeth_number
        self.a0 = a0 = np.radians(pressure_angle)
        self.clearance = clearance
        self.backlash = backlash
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
            self.twist_angle = width / \
                                (r0 * np.tan(np.pi / 2.0 - self.helix_angle))
        else:
            self.surface_splines = 2
            self.twist_angle = 0.0

        self.rim_r = rd + rim_width
        self.build_params = build_params

        # Calculate involute curve points for the left side of the tooth
        r = np.linspace(ra, rr, self.curve_points)
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
        self.t_tip_pts = np.dstack((np.cos(b) * rd,
                                    np.sin(b) * rd,
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
        t = np.linspace(t2 + np.pi * 2.0, t1 + np.pi * 2.0, self.curve_points)
        
        self.t_root_pts = np.dstack((bcxy[0] + bcr * np.cos(t),
                                     bcxy[1] + bcr * np.sin(t),
                                     np.zeros(self.curve_points))).squeeze()


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
        topface_rim_wire = cq.Wire.makeCircle(self.rim_r,
                                              cq.Vector(0.0, 0.0, 0.0),
                                              cq.Vector(0.0, 0.0, 1.0))
        topface = cq.Face.makeFromWires(topface_rim_wire, topface_wires)
        
        botface_wires = cq.Wire.combine(wp.edges('>Z').vals(),
                                        tol=self.wire_comb_tol)
        botface_rim_wire = cq.Wire.makeCircle(self.rim_r,
                                              cq.Vector(0.0, 0.0, self.width),
                                              cq.Vector(0.0, 0.0, 1.0))
        botface = cq.Face.makeFromWires(botface_rim_wire, botface_wires)
        wp = wp.add(topface).add(botface)
        wp = wp.add(self._build_rim_face())

        return wp.vals()


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
                      .moveTo(self.ra - E, self.width - wy)
                      .vLine(wy + E)
                      .hLine(wx + E)
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
                      .moveTo(self.ra + wx, -E)
                      .hLine(-wx - E)
                      .vLine(wy + E)
                      .close()
                      .revolve())

            body = (cq.Workplane('XY')
                    .add(body)
                    .cut(cutter))

        return body.val()


    def _build(self, chamfer=None, chamfer_top=None,
               chamfer_bottom=None, *args, **kv_args):
        faces = self._build_gear_faces()

        shell = make_shell(faces)
        body = cq.Solid.makeSolid(shell)
        
        body = self._make_chamfer(body, chamfer, chamfer_top, chamfer_bottom)
        
        return body


class HerringboneRingGear(RingGear):


    def _build_tooth_faces(self, twist_angle_a, twist_angle_b, z_pos, width):
        t_faces1 = (super(HerringboneRingGear, self)
                    ._build_tooth_faces(0.0, self.twist_angle,
                                        0.0, self.width / 2.0))
        t_faces2 = (super(HerringboneRingGear, self)
                    ._build_tooth_faces(self.twist_angle, 0.0,
                                        self.width / 2.0, self.width / 2.0))
        return t_faces1 + t_faces2


class PlanetaryGearset(GearBase):

    gear_cls = SpurGear
    ring_gear_cls = RingGear

    asm_sun_color = 'gold'
    asm_planet_color = 'lightsteelblue'
    asm_ring_color = 'goldenrod'

    def __init__(self, module, sun_teeth_number, planet_teeth_number, width,
                 rim_width, n_planets, pressure_angle=20.0,
                 helix_angle=0.0, clearance=0.0, backlash=0.0, **build_params):

        ring_z = sun_teeth_number + planet_teeth_number * 2

        self.sun = self.gear_cls(module, sun_teeth_number, width,
                                 pressure_angle=pressure_angle,
                                 helix_angle=helix_angle,
                                 clearance=clearance,
                                 backlash=backlash)

        self.planet = self.gear_cls(module, planet_teeth_number, width,
                                    pressure_angle=pressure_angle,
                                    helix_angle=-helix_angle,
                                    clearance=clearance,
                                    backlash=backlash)

        self.ring = self.ring_gear_cls(module, ring_z, width, rim_width,
                                       pressure_angle=pressure_angle,
                                       helix_angle=-helix_angle,
                                       clearance=clearance,
                                       backlash=backlash)

        self.orbit_r = self.sun.r0 + self.planet.r0
        self.n_planets = n_planets

        if ((sun_teeth_number + planet_teeth_number) % n_planets):
           if (sun_teeth_number % n_planets) or (ring_z % n_planets):
            warnings.warn('Planet gears being spaced evenly probably won\'t '
                          'mesh properly (if at all) with the given number of '
                          'teeth for the sun/planet gears and the number '
                          'of planets')

        self.build_params = build_params


    def _build(self, *args, **kv_args):
        asm = self.assemble(*args, **kv_args)
        return asm.toCompound()


    def assemble(self, build_sun=True, build_planets=True, build_ring=True,
                 sun_build_args={}, planet_build_args={},
                 ring_build_args={}, **kv_args):

        gearset = cq.Assembly(name='planetary')

        if build_sun:
            if 'sun_build_args' in self.build_params:
                in_args = self.build_params['sun_build_args']
            else:
                in_args = {}

            args = {**self.build_params,
                    **in_args,
                    **kv_args,
                    **sun_build_args}

            sun = self.sun.build(**args)
            
            if (self.planet.z % 2) != 0:
                loc = cq.Location(cq.Vector(0.0, 0.0, 0.0),
                                  cq.Vector(0.0, 0.0, 1.0),
                                  np.degrees(self.sun.tau / 2.0))
            else:
                loc = cq.Location()
            
            gearset.add(sun, name='sun', loc=loc,
                        color=cq.Color(self.asm_sun_color))


        if build_planets and self.n_planets > 0:
            planet_a = np.pi * 2.0 / self.n_planets

            if isinstance(build_planets, (list, tuple)):
                tobuild = build_planets
            else:
                tobuild = [True, ] * self.n_planets

            if 'planet_build_args' in self.build_params:
                in_args = self.build_params['planet_build_args']
            else:
                in_args = {}

            args = {**self.build_params,
                    **in_args,
                    **kv_args,
                    **planet_build_args}

            planet_body = self.planet.build(**args)

            planets = cq.Assembly(name='planets')
            
            for i, bld in enumerate(tobuild):
                if not bld:
                    continue

                loc = cq.Location(cq.Vector(np.cos(i * planet_a) * self.orbit_r,
                                            np.sin(i * planet_a) * self.orbit_r,
                                            0.0),
                                  cq.Vector(0.0, 0.0, 1.0),
                                  np.degrees(self.planet.tau / 2.0))
                
                planets.add(planet_body, name=f'planet_{i:02}', loc=loc,
                            color=cq.Color(self.asm_planet_color))
            gearset.add(planets)

        if build_ring:
            if 'ring_build_args' in self.build_params:
                in_args = self.build_params['ring_build_args']
            else:
                in_args = {}

            args = {**self.build_params,
                    **in_args,
                    **kv_args,
                    **ring_build_args}
            
            ring = self.ring.build(**args)
            loc = cq.Location(cq.Vector(0.0, 0.0, 0.0),
                              cq.Vector(0.0, 0.0, 1.0),
                              np.degrees(self.ring.tau / 2.0))
            gearset.add(ring, name='ring', loc=loc,
                        color=cq.Color(self.asm_ring_color))
            
        return gearset


class HerringbonePlanetaryGearset(PlanetaryGearset):
    gear_cls = HerringboneGear
    ring_gear_cls = HerringboneRingGear
