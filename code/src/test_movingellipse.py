from movingellipse import *

p = Parameters('test.ini')
s = Simulator(p)
res = s.run()

def test_projections_have_proper_dimensions():
    assert res.projections.shape[0] == p.Nt
    assert res.projections.shape[1] == p.imageSize