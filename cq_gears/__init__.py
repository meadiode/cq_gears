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


__version__='0.62'

import cadquery as cq
from .spur_gear import SpurGear, HerringboneGear, CustomTwistGear
from .ring_gear import (RingGear, HerringboneRingGear, PlanetaryGearset,
                        HerringbonePlanetaryGearset)
from .bevel_gear import BevelGear, BevelGearPair
from .rack_gear import RackGear, HerringboneRackGear
from .worm_gear import Worm, GloboidWorm
from .crossed_helical_gear import (CrossedHelicalGear, CrossedGearPair,
                                   HyperbolicGear, HyperbolicGearPair)

__all__ = [
    'SpurGear',
    'HerringboneGear',
    'CustomTwistGear',
    'RingGear',
    'HerringboneRingGear',
    'PlanetaryGearset',
    'HerringbonePlanetaryGearset',
    'BevelGear',
    'BevelGearPair',
    'RackGear',
    'HerringboneRackGear',
    'Worm',
    'GloboidWorm',
    'CrossedHelicalGear',
    'CrossedGearPair',
    'HyperbolicGear',
    'HyperbolicGearPair',
]


def gear(self, gear_, *build_args, **build_kv_args):
    '''Build a gear body from a gear object. This function to be monkey-patched
    into cadquery.Workplane, so cq_gears could be used as a plugin.
    self - cadquery.Workplane instance
    gear_ - an instance of one of the gear class - SpurGear, RingGear... etc
    build_args - positional arguments to pass to the gear's build function
    build_kv_args - key-value arguments to pass to the gear's build function
    return - modified cadquer.Workplane instance(self)
    '''    
    gear_body = gear_.build(*build_args, **build_kv_args)
    gears = self.eachpoint(lambda loc: gear_body.located(loc), True) 
    
    return gears


def addGear(self, gear_, *build_args, **build_kv_args):
    '''Same as Workplane.gear(...), but combines the resulting gear(s) with
       the other existing solids on the stack (if any). 
    '''    
    return self.union(gear(self, gear_, *build_args, **build_kv_args))


# Patch the functions into the Workplane
cq.Workplane.gear = gear
cq.Workplane.addGear = addGear
