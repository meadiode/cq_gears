
import cq_gears

TOL = 0.001


class TestRegression:


    def test(self, cls_name, expected, args):
        cls = getattr(cq_gears, cls_name)
        
        gear = cls(**args)
        body = gear.build()
    
        vol = body.Volume()
        bb = body.BoundingBox()

        assert vol - TOL <= expected['volume'] <= vol + TOL,\
              'Volume check failed'
        assert bb.xlen - TOL <= expected['bbox'][0] <= bb.xlen + TOL,\
              'Bounding box x-dim check failed'
        assert bb.ylen - TOL <= expected['bbox'][1] <= bb.ylen + TOL,\
              'Bounding box y-dim check failed'
        assert bb.zlen - TOL <= expected['bbox'][2] <= bb.zlen + TOL,\
              'Bounding box z-dim check failed'
