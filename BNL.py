#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 09:57:37 2018

@author: alexa
"""

def make_lattice(particle):
    import element as ent
    import lattice as ltc

    #%%
    ## lattice elements

    OD1 = ent.Drift(particle, 25e-2, "OD1")
    OD2 = ent.Drift(particle, 25e-2, "OD2")
    ORE = ent.Drift(particle, 2.17, "ORE")
    ORB = ent.Drift(particle, 2.2, "ORB")
    BPM = ent.Drift(particle, 15e-2, "BPM")

    QDA1 = ent.MQuad(particle, 5e-2, -10.23, 'QDA1')
    QFA1 = ent.MQuad(particle, 5e-2,  13.64, 'QFA1')
    QDA2 = ent.MQuad(particle, 5e-2, -8.60,  'QDA2')
    QFA2 = ent.MQuad(particle, 5e-2,  8.31, 'QFA2')


    OSF = ent.MSext(particle, 15e-2, 0, "OSF")
    OSD = ent.MSext(particle, 15e-2, 0, "OSD")

    # here set grads to 0 for now
    SDP = ent.MSext(particle, 15e-2, -3.39597809*0,"SDP")
    SFP = ent.MSext(particle, 15e-2,  2.7695769*0, "SFP")
    SDN = ent.MSext(particle, 15e-2,  3.79310524*0,"SDN")
    SFN = ent.MSext(particle, 15e-2, -2.09836542*0,"SFN")

    Obs = ent.Observer(particle, 'OUT')

    RBE = ent.CylWien(particle, 180.77969e-2, 120e5, 0.46002779, name="RBE")

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
             QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2,
             Obs] #TESTING OBSERVER

    SS2H2 = [QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2,
             QFA2, OD1, OSF, OD2, ORB, OD2, BPM, OD1, QDA2,
             QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2]

    ARC2 = [QFA1, OD1, SFP, OD2, RBE, OD2, BPM, OD1, QDA1,
            QDA1, OD1, SDP, OD2, RBE, OD2, BPM, OD1, QFA1]
    ARC2 = ARC2*8

    SS1H1 = [QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2,
             QDA2, OD1, SDP, OD2, ORB, OD2, BPM, OD1, QFA2,
             QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2]

    BNL_segments = dict(SS1H2=SS1H2, ARC1=ARC1, SS2H1=SS2H1,
                        SS2H2=SS2H2, ARC2=ARC2, SS1H1=SS1H1)

    #%%
    ## prepping lattice segments
    segments = list()
    for name, segment in  BNL_segments.items():
        segments.append(ltc.Lattice(segment, name))

    ## creating the E+B lattice
    lattice = ltc.Lattice(BNL_segments['SS1H2'],'SS1H2')
    for segment in segments[1:]:
        lattice = lattice + segment


    lattice = ltc.Lattice([ent.RF(particle, lattice.length, 75e3)], 'RF') + lattice
    lattice.name = 'BNL'

    return lattice

if __name__ == '__main__':
    from particle import Particle
    from state_list import StateList
    from lattice import track
    import matplotlib.pyplot as plt

    deu = Particle()
    lattice = make_lattice(deu)
    state = StateList(x = [-1e-3, 1e-3], d = [-.5e-4, 1e-4]).array

    log = track(state, lattice.TM('SS1H2'), 10)
