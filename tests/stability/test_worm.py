
import numpy as np
import cadquery as cq
from cq_gears import Worm

from utils import (BuildFailure, SolidBuildFailure, VolumeCheckFailure,
                   BBoxZCheckFailure, BBoxXYCheckFailure)
from utils import _TestGear, test_timeout


MODULE_MIN = 0.05
MODULE_MAX = 10.0

THREADS_MIN = 1
THREADS_MAX = 16

LENGTH_MIN = 10.0
LENGTH_MAX = 2000.0

LEAD_ANGLE_MIN = 0.5
ROOT_RADIUS_MIN = 1.0 

PRESSURE_ANGLE_MIN = 0.5
PRESSURE_ANGLE_MAX = 30.0

BBOX_CHECK_TOL = 0.5


class TestWorm(_TestGear):

    gear_cls = Worm
    argnames = ['module', 'lead_angle', 'n_threads', 'length', 'pressure_angle']

    @staticmethod
    def gen_params(seed, n):
        rng = np.random.default_rng(seed)

        module = MODULE_MIN + (MODULE_MAX - MODULE_MIN) * rng.random(n)
        n_threads = rng.integers(THREADS_MIN, THREADS_MAX + 1, n)
        lead_angle_max = (n_threads * module) / \
                         (2.0 * Worm.kd * module + 2.0 * ROOT_RADIUS_MIN)
        lead_angle = LEAD_ANGLE_MIN + \
                        (lead_angle_max - LEAD_ANGLE_MIN) * rng.random(n)
        length = LENGTH_MIN + (LENGTH_MAX - LENGTH_MIN) * rng.random(n)
        pressure_angle = PRESSURE_ANGLE_MIN + \
                    (PRESSURE_ANGLE_MAX - PRESSURE_ANGLE_MIN) * rng.random(n)
        
        params = []

        for vals in zip(module, lead_angle, n_threads, length, pressure_angle):
            params.append({k : v.item() for k, v in \
                            zip(TestWorm.argnames, vals)})

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
        vmax = gear.length * np.pi * gear.ra ** 2
        vmin = gear.length * np.pi * gear.rd ** 2

        if not (vmin < body.Volume() < vmax):
            out_q.put(VolumeCheckFailure())
            return

        # Bounding box check
        bb = body.BoundingBox()

        if abs(gear.length - (bb.xmax - bb.xmin)) > BBOX_CHECK_TOL:
            out_q.put(BBoxZCheckFailure())
            return

        maxd = max((bb.zmax - bb.zmin), (bb.ymax - bb.ymin))
        
        if abs(gear.ra * 2.0 - maxd) > BBOX_CHECK_TOL:
            out_q.put(BBoxXYCheckFailure())
            return

        out_q.put(None)
