
import pytest
import multiprocessing as mp
import queue


@pytest.fixture
def test_timeout(request):
    return request.config.getoption("--test_timeout")


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



class _TestGear:

    argnames = []
    gear_cls = None


    @staticmethod
    def gen_params(seed, n):
        raise NotImplementedError


    @staticmethod
    def gear_build(out_q, gear):
        raise NotImplementedError


    def test(self, gear_params, json_metadata, test_timeout):
        json_metadata['gear_params'] = gear_params

        try:
            gear = self.gear_cls(**gear_params)
        except:
            json_metadata['failure_tag'] = PrecalcFailure.tag
            raise PrecalcFailure()


        fmp = mp.get_context(method='fork')

        proc_res = fmp.Queue()
        proc = fmp.Process(target=self.gear_build, args=(proc_res, gear))

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
