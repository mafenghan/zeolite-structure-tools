import numpy as np
from zeolite_structure_tools.geometry import xyz_to_rsf
def test_z_axis():
    r,s,f=xyz_to_rsf([0,0,1]); assert np.isclose(r,1); assert np.isclose(s,0); assert np.isclose(f,0)
def test_negative_x():
    _,s,_=xyz_to_rsf([-1,0,0]); assert np.isclose(abs(s),np.pi)
