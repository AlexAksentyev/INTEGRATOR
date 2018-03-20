import numpy as np
import math

pi = np.pi
R = lambda x: np.array([[math.cos(x) , -math.sin(x)], [math.sin(x), math.cos(x)]])
s0 = np.array([0, 1])
s1 = R(.3).dot(s0)
s1[0] += .02
s2 = R(pi/2-.23).dot(s0)
s3 = R(pi-.05).dot(s0)
s4 = R(-pi + .4).dot(s0)
s5 = R(-pi/2 + .24).dot(s0)
s = np.array([s1,s2,s3,s4,s5])

s00 = np.array([s0,s0,s0,s0,s0])
cos_psy = s.dot(s0)
sin_phi = np.cross(s, s00)
# sin_phi *= np.sign(cos_psy)
for i, s_i in enumerate(s):
    Ry = np.array([[cos_psy[i], -sin_phi[i]], [sin_phi[i], cos_psy[i]]])
    print(Ry.dot(s_i))
