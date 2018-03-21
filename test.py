import numpy as np
import math
from pandas import DataFrame

rad = np.deg2rad
deg = np.rad2deg

pi = np.pi
R = lambda x: np.array([[np.cos(rad(x)) , -np.sin(rad(x))], [np.sin(rad(x)), np.cos(rad(x))]])
s0 = np.array([0, 1])

arg = -np.array([45, 135, 225, 315, -180])
s1 = R(arg[0]).dot(s0)
s2 = R(arg[1]).dot(s0)
s3 = R(arg[2]).dot(s0)
s4 = R(arg[3]).dot(s0)
s5 = R(arg[4]).dot(s0)
s = np.array([s1,s2,s3,s4,s5])

s00 = np.array([s0,s0,s0,s0,s0])
cos_psy = s.dot(s0)
sin_phi = np.cross(s, s00)

for i, s_i in enumerate(s):
    sc = np.sign(cos_psy[i])
    c = sc*cos_psy[i]
    s = sc*sin_phi[i]
    Ry = np.array([[c, -s], [s, c]])
    #print((s_i, Ry.dot(Ry.dot(s_i))))
    print(Ry.dot(s_i))

arg_cos = deg(np.arccos(cos_psy))
arg_sin = deg(np.arcsin(sin_phi))
cw = np.sign(cos_psy*sin_phi)
print(DataFrame({'Atheta':arg, #'Ctheta':arg_cos, 'Stheta':arg_sin,
                 'cos':cos_psy, 'sin':sin_phi, 'Ttheta':cw}) )
