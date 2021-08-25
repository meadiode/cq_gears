
__version__='0.45'

import cadquery as cq
from .spur_gear import SpurGear, HerringboneGear
from .ring_gear import (RingGear, HerringboneRingGear, PlanetaryGearset,
                        HerringbonePlanetaryGearset)
from .bevel_gear import BevelGear, BevelGearPair

__all__ = [
    'SpurGear',
    'HerringboneGear',
    'RingGear',
    'HerringboneRingGear',
    'PlanetaryGearset',
    'HerringbonePlanetaryGearset',
    'BevelGear',
    'BevelGearPair',
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


def add_gear(self, gear_, *build_args, **build_kv_args):
    '''Same as Workplane.gear(...), but combines the resulting gear(s) with
       the other existing solids on the stack (if any). 
    '''    
    return self.union(gear(self, gear_, *build_args, **build_kv_args))


# Patch the functions into the Workplane
cq.Workplane.gear = gear
cq.Workplane.add_gear = add_gear
