
import numpy as np
import cadquery as cq
from cq_gears import RackGear, HerringboneRackGear

from utils import (BuildFailure, SolidBuildFailure, VolumeCheckFailure,
                   BBoxZCheckFailure, BBoxXYCheckFailure)
from utils import _TestGear, test_timeout


MODULE_MIN = 0.05
MODULE_MAX = 10.0

LENGTH_MIN = 1.0
LENGTH_MAX = 2000.0

FACE_WIDTH_MIN = 0.1
FACE_WIDTH_MAX = 1000.0

HEIGHT_MIN = 0.1
HEIGHT_MAX = 100.0

PRESSURE_ANGLE_MIN = 0.5
PRESSURE_ANGLE_MAX = 30.0

HELIX_ANGLE_MIN = -70.0
HELIX_ANGLE_MAX = 70.0

BBOX_CHECK_TOL = 0.5


class TestRackGear(_TestGear):

    gear_cls = RackGear
    argnames = ['module', 'length', 'width', 'height', 'pressure_angle',
                'helix_angle']

    @staticmethod
    def gen_params(seed, n):
        rng = np.random.default_rng(seed)

        module = MODULE_MIN + (MODULE_MAX - MODULE_MIN) * rng.random(n)
        length = LENGTH_MIN + (LENGTH_MAX - LENGTH_MIN) * rng.random(n)
        width = FACE_WIDTH_MIN + \
                (FACE_WIDTH_MAX - FACE_WIDTH_MIN) * rng.random(n)
        height = HEIGHT_MIN + (HEIGHT_MAX - HEIGHT_MIN) * rng.random(n)
        pressure_angle = PRESSURE_ANGLE_MIN + \
                    (PRESSURE_ANGLE_MAX - PRESSURE_ANGLE_MIN) * rng.random(n)
        helix_angle = HELIX_ANGLE_MIN + \
                      (HELIX_ANGLE_MAX - HELIX_ANGLE_MIN) * rng.random(n)
        
        params = []

        for vals in zip(module, length, width, height,
                        pressure_angle, helix_angle):
            params.append({k : v.item() for k, v in \
                                    zip(TestRackGear.argnames, vals)})

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
        vmax = gear.length * gear.width * (gear.height + gear.tooth_height)
        vmin = gear.length * gear.width * gear.height

        if not (vmin < body.Volume() < vmax):
            out_q.put(VolumeCheckFailure())
            return

        # Bounding box check
        bb = body.BoundingBox()

        if abs(gear.width - (bb.zmax - bb.zmin)) > BBOX_CHECK_TOL:
            out_q.put(BBoxZCheckFailure())
            return

        
        ydim = gear.height + gear.tooth_height

        if abs(ydim - (bb.ymax - bb.ymin)) > BBOX_CHECK_TOL or \
           abs(gear.length - (bb.xmax - bb.xmin)) > BBOX_CHECK_TOL:        
            out_q.put(BBoxXYCheckFailure())
            return

        out_q.put(None)


class TestHerringboneRackGear(TestRackGear):
    gear_cls = HerringboneRackGear
