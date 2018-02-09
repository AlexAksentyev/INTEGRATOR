#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 09:57:37 2018

@author: alexa
"""

import element as ent
import lattice as ltc
from tracker import Tracker
from particle import Particle

trkr = Tracker()
deu = Particle()

#%%
## lattice elements

OD1 = ent.Drift(25e-2, "OD1")
OD2 = ent.Drift(25e-2, "OD2")
ORE = ent.Drift(2.17, "ORE")
ORB = ent.Drift(2.2, "ORB")

QDA1 = ent.MQuad(5e-2, -10.23, 'QDA1')
QFA1 = ent.MQuad(5e-2,  13.64, 'QFA1')
QDA2 = ent.MQuad(5e-2, -8.60,  'QDA2')
QFA2 = ent.MQuad(5e-2,  8.31, 'QFA2')


OSF = ent.MSext(15e-2, 0, "OSF")
OSD = ent.MSext(15e-2, 0, "OSD")

# here set grads to 0 for now
SDP = ent.MSext(15e-2, -3.39597809*0,"SDP")
SFP = ent.MSext(15e-2,  2.7695769*0, "SFP")
SDN = ent.MSext(15e-2,  3.79310524*0,"SDN")
SFN = ent.MSext(15e-2, -2.09836542*0,"SFN")

BPM = ent.Drift(15e-2, "BPM")

RBE = ent.Wien(180.77969e-2, 5e-2, deu, -120e5, 0.46002779, name="RBE")

#%%
## definition of lattice subsections

SS1H2 = [QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2,
         QFA2, OD1, OSF, OD2, ORB, OD2, BPM, OD1, QDA2,
         QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2]

ARC1 = [QFA1, OD1, SFP, OD2, RBE, OD2, BPM, OD1, QDA1,
        QDA1, OD1, SDP, OD2, RBE, OD2, BPM, OD1, QFA1]
ARC1 = ARC1*8

SS2H1 = [QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDP, OD2, ORB, OD2, BPM, OD1, QFA2,
         QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2]

SS2H2 = [QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2,
         QFA2, OD1, OSF, OD2, ORB, OD2, BPM, OD1, QDA2,
         QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2]

ARC2 = [QFA1, OD1, SFP, OD2, RBE, OD2, BPM, OD1, QDA1,
        QDA1, OD1, SDP, OD2, RBE, OD2, BPM, OD1, QFA1]
ARC2 = ARC2*8

SS1H1 = [QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDP, OD2, ORB, OD2, BPM, OD1, QFA2,
         QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2]

BNL_segments = dict(SS1H2=SS1H2, ARC1=ARC1, SS2H1=SS2H2,
                    SS2H2=SS2H2, ARC2=ARC2, SS1H1=SS1H1)

#%%
if __name__ == '__main__':
     ## prepping lattice segments
    segments = list()
    for name, segment in  BNL_segments.items():
        segments.append(ltc.Lattice(segment, name))

    ## creating the E+B lattice
    lattice = ltc.Lattice(BNL_segments['SS1H2'],'SS1H2')
    for segment in segments[1:]:
        lattice = lattice + segment

    lattice.name = 'BNL'
    lattice.insert_RF(0, 0, deu)
    #%%
    from particle_log import StateList

    n_turns = int(5)

    bunch = StateList(Sz=1, x=(-1e-3,1e-3,3), dK=(0, 1e-4, 2))
    log = trkr.track(deu, bunch, lattice, n_turns)