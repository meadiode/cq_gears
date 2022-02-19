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
from .spur_gear import GearBase, SpurGear
from .utils import circle3d_by3points

class CrossedHelicalGear(SpurGear):
    
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


        at0 = np.arctan(a0 / np.cos(self.helix_angle))
        mt = m / np.cos(self.helix_angle)
        
        d0 = mt * z         # pitch diameter
        adn = self.ka / (z / d0) # addendum
        ddn = self.kd / (z / d0) # dedendum
        da = d0 + 2.0 * adn # addendum circle diameter
        dd = d0 - 2.0 * ddn - 2.0 * clearance # dedendum circle diameter
        inv_a0 = np.tan(at0) - at0

        self.r0 = r0 = d0 / 2.0 # pitch radius
        self.ra = ra = da / 2.0 # addendum radius
        self.rd = rd = dd / 2.0 # dedendum radius
        self.rb = rb = np.cos(at0) * r0 # base circle radius
        self.rr = rr = max(rb, rd) # tooth root radius
        self.tau = tau = np.pi * 2.0 / z # pitch angle

        s0 =  r0 * np.pi / z # tooth thickness on the pitch circle
        
        
        if helix_angle != 0.0:
            self.twist_angle = width / \
                                (r0 * np.tan(np.pi / 2.0 - self.helix_angle))
        else:
            self.surface_splines = 2
            self.twist_angle = 0.0

        self.build_params = build_params

        # Calculate involute curve points for the left side of the tooth
        r = np.linspace(rr, ra, self.curve_points)
        cos_a = r0 / r * np.cos(at0)
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



class CrossedGearPair(GearBase):

    gear1_cls = CrossedHelicalGear
    gear2_cls = CrossedHelicalGear
    asm_gear1_color = 'goldenrod'
    asm_gear2_color = 'lightsteelblue'
    
    def __init__(self, module, gear1_teeth_number, gear2_teeth_number,
                 gear1_width, gear2_width, pressure_angle=20.0,
                 shaft_angle=90.0, gear1_helix_angle=None, clearance=0.0,
                 backlash=0.0, **build_params):
        
        if gear1_helix_angle is None:
            g1_helix = shaft_angle / 2.0
            g2_helix = shaft_angle / 2.0
        else:            
            g1_helix = gear1_helix_angle
            g2_helix = shaft_angle - gear1_helix_angle
            
        
        self.gear1 = self.gear1_cls(module, gear1_teeth_number, gear1_width,
                                    pressure_angle=pressure_angle,
                                    helix_angle=g1_helix,
                                    clearance=clearance,
                                    backlash=backlash)
        
        self.gear2 = self.gear2_cls(module, gear2_teeth_number, gear2_width,
                                    pressure_angle=pressure_angle,
                                    helix_angle=g2_helix,
                                    clearance=clearance,
                                    backlash=backlash)
        self.shaft_angle = np.radians(shaft_angle)
        self.build_params = build_params


    def assemble(self, build_gear1=True, build_gear2=True, transform_gear2=True,
                 gear1_build_args={}, gear2_build_args={}, **kv_args):
        
        gearset = cq.Assembly(name='crossed_pair')

        if build_gear1:
            args = {**self.build_params, **kv_args, **gear1_build_args}
            gear1 = self.gear1.build(**args)

            gearset.add(gear1, name='gear1', loc=cq.Location(),
                        color=cq.Color(self.asm_gear1_color))

        if build_gear2:
            args = {**self.build_params, **kv_args, **gear2_build_args}
            gear2 = self.gear2.build(**args)

            if transform_gear2:
                ratio = self.gear1.z / self.gear2.z
                align_angle = 0.0 if self.gear2.z % 2 else 180.0 / self.gear2.z            
                align_angle += np.degrees(self.gear2.twist_angle + \
                                          self.gear1.twist_angle * ratio) / 2.0
            
                loc = cq.Location(cq.Vector(self.gear1.r0 + self.gear2.r0,
                                            0.0,
                                            self.gear1.width / 2.0))            
                loc *= cq.Location(cq.Vector(0.0, 0.0, 0.0),
                                   cq.Vector(1.0, 0.0, 0.0),
                                   np.degrees(self.shaft_angle))
                loc *= cq.Location(cq.Vector(0.0, 0.0, -self.gear2.width / 2.0))
                loc *= cq.Location(cq.Vector(0.0, 0.0, 0.0),
                                   cq.Vector(0.0, 0.0, 1.0),
                                   align_angle)
                
            else:
                loc = cq.Location()

            gearset.add(gear2, name='gear2', loc=loc,
                        color=cq.Color(self.asm_gear2_color))

        return gearset
        
        
    def _build(self, *args, **kv_args):
        asm = self.assemble(*args, **kv_args)
        return asm.toCompound()


class HyperbolicGear(SpurGear):

    surface_splines = 2
    
    def __init__(self, module, teeth_number, width, twist_angle,
                 pressure_angle=20.0, clearance=0.0,
                 backlash=0.0, **build_params):
        
        super(HyperbolicGear, self).__init__(module=module,
                                             teeth_number=teeth_number,
                                             width=width,
                                             pressure_angle=pressure_angle,
                                             helix_angle=0.0,
                                             clearance=clearance,
                                             backlash=backlash,
                                             **build_params)
        self.twist_angle = np.radians(twist_angle)
        
        ln = np.cos(self.twist_angle) * self.r0
        ht = np.sin(self.twist_angle) * self.r0
        rpx = (self.r0 + ln) / 2.0
        rpy = ht / 2.0
        self.throat_r = np.sqrt(rpx ** 2 + rpy ** 2)


class HyperbolicGearPair(GearBase):

    gear_cls = HyperbolicGear
    asm_gear1_color = 'goldenrod'
    asm_gear2_color = 'lightsteelblue'
    
    def __init__(self, module, gear1_teeth_number, width, shaft_angle,
                 gear2_teeth_number=None, pressure_angle=20.0,
                 clearance=0.0, backlash=0.0, **build_params):
        
        if gear2_teeth_number is None:
            gear2_teeth_number = gear1_teeth_number
        
        g1_r0 = module * gear1_teeth_number / 2.0
        g2_r0 = module * gear2_teeth_number / 2.0
        
        alpha = np.radians(shaft_angle / 2.0)
        hh = (width / 2.0) * np.tan(alpha)
        
        gear1_twist_angle = np.arcsin(hh / g1_r0) * 2.0
        gear2_twist_angle = np.arcsin(hh / g2_r0) * 2.0
        
        if np.isnan(gear1_twist_angle) or np.isnan(gear2_twist_angle):
            raise ValueError('Impossible to calculate the twist angle for the '
                             'given shaft angle / teeth number / gear width')
        
        self.shaft_angle = np.radians(shaft_angle)
                    
        self.gear1 = self.gear_cls(module, gear1_teeth_number, width,
                                   twist_angle=np.degrees(gear1_twist_angle),
                                   pressure_angle=pressure_angle,
                                   clearance=clearance,
                                   backlash=backlash)
        
        self.gear2 = self.gear_cls(module, gear2_teeth_number, width,
                                   twist_angle=np.degrees(gear2_twist_angle),
                                   pressure_angle=pressure_angle,
                                   clearance=clearance,
                                   backlash=backlash)
        self.build_params = build_params
        

    def assemble(self, build_gear1=True, build_gear2=True, transform_gear2=True,
                 gear1_build_args={}, gear2_build_args={}, **kv_args):

        gearset = cq.Assembly(name='hyperbolic_pair')

        if build_gear1:
            args = {**self.build_params, **kv_args, **gear1_build_args}
            gear1 = self.gear1.build(**args)
            gearset.add(gear1, name='gear1', loc=cq.Location(),
                        color=cq.Color(self.asm_gear1_color))

        if build_gear2:
            args = {**self.build_params, **kv_args, **gear2_build_args}
            gear2 = self.gear2.build(**args)

            if transform_gear2:

                ratio = self.gear1.z / self.gear2.z
                align_angle = 0.0 if self.gear2.z % 2 else 180.0 / self.gear2.z
                align_angle += np.degrees(self.gear2.twist_angle + \
                                          self.gear1.twist_angle * ratio) / 2.0
                
                loc = cq.Location(cq.Vector(self.gear1.throat_r + \
                                            self.gear2.throat_r,
                                            0.0,
                                            self.gear1.width / 2.0))
                loc *= cq.Location(cq.Vector(0.0, 0.0, 0.0),
                                   cq.Vector(1.0, 0.0, 0.0),
                                   np.degrees(self.shaft_angle))
                loc *= cq.Location(cq.Vector(0.0, 0.0, -self.gear2.width / 2.0))
                loc *= cq.Location(cq.Vector(0.0, 0.0, 0.0),
                                   cq.Vector(0.0, 0.0, 1.0),
                                   align_angle)
            else:
                loc = cq.Location()

            gearset.add(gear2, name='gear2', loc=loc,
                        color=cq.Color(self.asm_gear2_color))

        return gearset
        

    def _build(self, *args, **kv_args):
        asm = self.assemble(*args, **kv_args)
        return asm.toCompound()
