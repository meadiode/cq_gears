
import pytest
import numpy as np
import cadquery as cq
from cq_gears import SpurGear


MODULE_MIN = 0.05
MODULE_MAX = 10.0

TEETH_NUMBER_MIN = 3
TEETH_NUMBER_MAX = 80

FACE_WIDTH_MIN = 0.1
FACE_WIDTH_MAX = 1000.0

PRESSURE_ANGLE_MIN = 0.5
PRESSURE_ANGLE_MAX = 30.0

HELIX_ANGLE_MIN = -70.0
HELIX_ANGLE_MAX = 70.0

BBOX_CHECK_TOL = 0.5


class TestSpurGear:

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
            params.append({k : v for k, v in zip(TestSpurGear.argnames, vals)})

        return params


    def test_spur_gear(self, module, teeth_number, width,
                       pressure_angle, helix_angle):
        try:
            gear = SpurGear(module, teeth_number, width,
                            pressure_angle, helix_angle)
        except:
            raise Exception('PRECALCULATION')
        
        try:
            body = gear.build()
        except:
            raise Exception('BUILDING')

        # The result should be a solid body
        assert isinstance(body, cq.Solid), 'NOT_SOLID'

        # Check roughly the resulting volume of the gear
        vmax = gear.width * np.pi * gear.ra ** 2
        vmin = gear.width * np.pi * gear.rd ** 2
        assert vmin < body.Volume() < vmax, 'VOLUME'

        # Bounding box check
        bb = body.BoundingBox()
        assert abs(gear.width - (bb.zmax - bb.zmin)) <= BBOX_CHECK_TOL, 'BBOX_Z'

        maxd = max((bb.xmax - bb.xmin), (bb.ymax - bb.ymin))
        assert abs(gear.ra * 2.0 - maxd) <= BBOX_CHECK_TOL, 'BBOX_XY'
