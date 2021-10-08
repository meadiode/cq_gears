import pytest
import numpy as np
import cadquery as cq
from cq_gears import BevelGear
import multiprocessing as mp
import queue


from utils import (BuildFailure, SolidBuildFailure, VolumeCheckFailure,
                    BBoxZCheckFailure, BBoxXYCheckFailure)
from utils import _TestGear, test_timeout


MODULE_MIN = 0.05
MODULE_MAX = 10.0

TEETH_NUMBER_MIN = 3
TEETH_NUMBER_MAX = 200

CONE_ANGLE_MIN = 1.0
CONE_ANGLE_MAX = 89.0

PRESSURE_ANGLE_MIN = 0.5
PRESSURE_ANGLE_MAX = 30.0

HELIX_ANGLE_MIN = -70.0
HELIX_ANGLE_MAX = 70.0

BBOX_CHECK_TOL = 0.5


class TestBevelGear(_TestGear):

    gear_cls = BevelGear
    
    argnames = ['module', 'teeth_number', 'cone_angle', 'face_width',
                 'pressure_angle', 'helix_angle']

    @staticmethod
    def gen_params(seed, n):
        rng = np.random.default_rng(seed)

        module = MODULE_MIN + (MODULE_MAX - MODULE_MIN) * rng.random(n)
        teeth_number = rng.integers(TEETH_NUMBER_MIN, TEETH_NUMBER_MAX + 1, n)

        cone_angle =  CONE_ANGLE_MIN + \
                      (CONE_ANGLE_MAX - CONE_ANGLE_MIN) * rng.random(n)


        # Generating random face_width'
        gamma_p = np.radians(cone_angle)
        rp = module * teeth_number / 2.0        
        gs_r = rp / np.sin(gamma_p)

        face_width =  gs_r * rng.random(n)

        pressure_angle = PRESSURE_ANGLE_MIN + \
                    (PRESSURE_ANGLE_MAX - PRESSURE_ANGLE_MIN) * rng.random(n)
        helix_angle = HELIX_ANGLE_MIN + \
                      (HELIX_ANGLE_MAX - HELIX_ANGLE_MIN) * rng.random(n)
        
        params = []

        for vals in zip(module, teeth_number, cone_angle, face_width,
                        pressure_angle, helix_angle):
            params.append({k : v.item() for k, v in zip(TestBevelGear.argnames, vals)})

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
        # vmax = gear.width * np.pi * gear.ra ** 2
        # vmin = gear.width * np.pi * gear.rd ** 2

        # if not (vmin < body.Volume() < vmax):
        #     out_q.put(VolumeCheckFailure())
        #     return

        # # Bounding box check
        # bb = body.BoundingBox()

        # if abs(gear.width - (bb.zmax - bb.zmin)) > BBOX_CHECK_TOL:
        #     out_q.put(BBoxZCheckFailure())
        #     return

        # maxd = max((bb.xmax - bb.xmin), (bb.ymax - bb.ymin))
        
        # if abs(gear.ra * 2.0 - maxd) > BBOX_CHECK_TOL:
        #     out_q.put(BBoxXYCheckFailure())
        #     return

        out_q.put(None)




    # def test_bevel_gear(self, gear_params, json_metadata, test_timeout):
    #     json_metadata['gear_params'] = gear_params

    #     try:
    #         gear = BevelGear(**gear_params)
    #     except:
    #         json_metadata['failure_tag'] = PrecalcFailure.tag
    #         raise PrecalcFailure()
        
        
    #     def gear_build(out_q):
    #         try:
    #             body = gear.build()
    #         except:
    #             out_q.put(BuildFailure())
    #             return

    #         # The result should be a solid body
    #         if not isinstance(body, cq.Solid):
    #             out_q.put(SolidBuildFailure())
    #             return

    #         # Check roughly the resulting volume of the gear
    #         # vmax = gear.width * np.pi * gear.ra ** 2
    #         # vmin = gear.width * np.pi * gear.rd ** 2

    #         # if not (vmin < body.Volume() < vmax):
    #         #     out_q.put(VolumeCheckFailure())
    #         #     return

    #         # # Bounding box check
    #         # bb = body.BoundingBox()

    #         # if abs(gear.width - (bb.zmax - bb.zmin)) > BBOX_CHECK_TOL:
    #         #     out_q.put(BBoxZCheckFailure())
    #         #     return

    #         # maxd = max((bb.xmax - bb.xmin), (bb.ymax - bb.ymin))
            
    #         # if abs(gear.ra * 2.0 - maxd) > BBOX_CHECK_TOL:
    #         #     out_q.put(BBoxXYCheckFailure())
    #         #     return

    #         out_q.put(None)


    #     fmp = mp.get_context(method='fork')
        
    #     proc_res = fmp.Queue()
    #     proc = fmp.Process(target=gear_build, args=(proc_res,))

    #     proc.start()

    #     try:
    #         res = proc_res.get(timeout=test_timeout)
    #     except (fmp.TimeoutError, queue.Empty):
    #         proc.terminate()
    #         json_metadata['failure_tag'] = TimeoutFailure.tag
    #         raise TimeoutFailure()

    #     proc.join()

    #     if isinstance(res, Exception):
    #         json_metadata['failure_tag'] = res.tag
    #         raise res


if __name__ == '__main__':
    print(TestBevelGear.gen_params(42, 3))

