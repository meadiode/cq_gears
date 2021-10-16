
import numpy as np
import cadquery as cq
from cq_gears import SpurGear, HerringboneGear

from utils import (BuildFailure, SolidBuildFailure, VolumeCheckFailure,
                   BBoxZCheckFailure, BBoxXYCheckFailure)
from utils import _TestGear, test_timeout


MODULE_MIN = 0.05
MODULE_MAX = 10.0

TEETH_NUMBER_MIN = 3
TEETH_NUMBER_MAX = 200

FACE_WIDTH_MIN = 0.1
FACE_WIDTH_MAX = 1000.0

PRESSURE_ANGLE_MIN = 0.5
PRESSURE_ANGLE_MAX = 30.0

HELIX_ANGLE_MIN = -70.0
HELIX_ANGLE_MAX = 70.0

BBOX_CHECK_TOL = 0.5


class TestSpurGear(_TestGear):

    gear_cls = SpurGear
    argnames = ['module', 'teeth_number', 'width', 'pressure_angle',
                'helix_angle']

    @staticmethod
    def gen_params(seed, n):
        rng = np.random.default_rng(seed)

        module = MODULE_MIN + (MODULE_MAX - MODULE_MIN) * rng.random(n)
        teeth_number = rng.integers(TEETH_NUMBER_MIN, TEETH_NUMBER_MAX + 1, n)
        width = FACE_WIDTH_MIN + \
                (FACE_WIDTH_MAX - FACE_WIDTH_MIN) * rng.random(n)
        pressure_angle = PRESSURE_ANGLE_MIN + \
                    (PRESSURE_ANGLE_MAX - PRESSURE_ANGLE_MIN) * rng.random(n)
        helix_angle = HELIX_ANGLE_MIN + \
                      (HELIX_ANGLE_MAX - HELIX_ANGLE_MIN) * rng.random(n)
        
        params = []

        for vals in zip(module, teeth_number, width,
                        pressure_angle, helix_angle):
            params.append({k : v.item() for k, v in zip(TestSpurGear.argnames, vals)})

        return params


    @staticmethod
    def gear_build(out_q, gear):
        try:
            body = gear.build()
        except:
            out_q.put(BuildFailure())
            return

        # The result should be a solid body
        if not isinstance(body, cq.Solid):
            out_q.put(SolidBuildFailure())
            return

        # Check roughly the resulting volume of the gear
        vmax = gear.width * np.pi * gear.ra ** 2
        vmin = gear.width * np.pi * gear.rd ** 2

        if not (vmin < body.Volume() < vmax):
            out_q.put(VolumeCheckFailure())
            return

        # Bounding box check
        bb = body.BoundingBox()

        if abs(gear.width - (bb.zmax - bb.zmin)) > BBOX_CHECK_TOL:
            out_q.put(BBoxZCheckFailure())
            return

        maxd = max((bb.xmax - bb.xmin), (bb.ymax - bb.ymin))
        
        if abs(gear.ra * 2.0 - maxd) > BBOX_CHECK_TOL:
            out_q.put(BBoxXYCheckFailure())
            return

        out_q.put(None)


class TestHerringboneGear(TestSpurGear):
    gear_cls = HerringboneGear
