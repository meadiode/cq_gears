
import pytest
import numpy as np
import cadquery as cq
from cq_gears import SpurGear
import multiprocessing as mp
import queue


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


class CQGearsFailure(Exception):
    tag = 'UNDEFINED'    


class PrecalcFailure(CQGearsFailure):
    tag = 'PRECALC'


class BuildFailure(CQGearsFailure):
    tag = 'BUILDING'


class SolidBuildFailure(CQGearsFailure):
    tag = 'NOT_SOLID'


class VolumeCheckFailure(CQGearsFailure):
    tag = 'VOLUME'


class BBoxZCheckFailure(CQGearsFailure):
    tag = 'BBOX_Z'


class BBoxXYCheckFailure(CQGearsFailure): 
    tag = 'BBOX_XY'


class TimeoutFailure(CQGearsFailure): 
    tag = 'TIMEOUT'


@pytest.fixture
def test_timeout(request):
    return request.config.getoption("--test_timeout")


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
            params.append({k : v.item() for k, v in zip(TestSpurGear.argnames, vals)})

        return params


    def test_spur_gear(self, gear_params, json_metadata, test_timeout):
        json_metadata['gear_params'] = gear_params

        try:
            gear = SpurGear(**gear_params)
        except:
            json_metadata['failure_tag'] = PrecalcFailure.tag
            raise PrecalcFailure()
        
        
        def gear_build(out_q):
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


        fmp = mp.get_context(method='fork')
        
        proc_res = fmp.Queue()
        proc = fmp.Process(target=gear_build, args=(proc_res,))

        proc.start()

        try:
            res = proc_res.get(timeout=test_timeout)
        except (fmp.TimeoutError, queue.Empty):
            proc.terminate()
            json_metadata['failure_tag'] = TimeoutFailure.tag
            raise TimeoutFailure()

        proc.join()

        if isinstance(res, Exception):
            json_metadata['failure_tag'] = res.tag
            raise res
