import pyximport; pyximport.install(pyimport=True)

import CParticle as PCL
import CElement as ENT
import utilFunc as U

from time import clock

#%%

# hardware parameters
Lq = 5
Ls = .15

GSFP = 0 
GSDP = 0

Lw = 361.55403e-2
B = .082439761
E = -120e5
#%%
# lattice elements

DL_25  = ENT.Drift(.25)
DL_15 = ENT.Drift(.15)
DL2_2 = ENT.Drift(2.2)
BPM = ENT.Drift(15)

QDS = ENT.MQuad(Lq, -.86)
QFS = ENT.MQuad(Lq, .831)

QDA = ENT.MQuad(Lq, -1.023)
QFA = ENT.MQuad(Lq, 1.364)

Sf = ENT.MSext(Ls, GSFP)
Sd = ENT.MSext(Ls, GSDP)

R3 = ENT.Wien(Lw,5e-2,PCL.Particle([0]),E,B)

#%%
# lattice definition

SS1H2 = [QDS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS , QFS , DL_25 , DL_15 , ENT.Element(0,0) , #  RF ,
                                     DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QDS , QDS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS]


StateList = U.form_state_list((0e-3,0e-3),(-0e-3,0e-3),3,1)
E = PCL.Ensemble.from_state(StateList)
E.setReference(0)
n = E.count()-1
ddk = 2e-4/n
for i in range(1,E.count()):
    E[i].set(dK=3e-4-(i-1)*ddk)

start = clock()
E.track(SS1H2,1000,inner=True)
print('Time passed: {} sec'.format(clock()-start))