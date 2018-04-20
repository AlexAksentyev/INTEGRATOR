#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 13:41:04 2018

@author: alexa
"""

import element as ent
import lattice as ltc
from tracker import Tracker
from particle import Particle
from time import clock
from particle_log import StateList
from lattice import Lattice


def make_lattice(particle):
    QDA2 = ent.MQuad(5e-2, -8.2)
    QFA2 = ent.MQuad(5e-2, 7.36)
    OD1 = ent.Drift(25e-2)
    ORB = ent.Drift(220e-2)
    OD2 = ent.Drift(25e-2)
    OSF = ent.Drift(15e-2)
    OSD = ent.Drift(15e-2)
    BPM = ent.Drift(15e-2)

    lattice = Lattice([QDA2, OD1, ORB, OD2, BPM, OD1, QFA2,
                       QFA2, OD1, OSF, OD2, ORB, OD2, BPM, OD1, QDA2,
                       QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2], name='FODO_lattice')

    # lattice.insert_RF(0, 0, particle, E_field=15e7)
    return lattice

if __name__ == '__main__':
    trkr = Tracker()
    deu = Particle()

    #%%

    ensemble = StateList(x=(-1e-3, 1e-3, 5), dK=(0, 1e-4, 5), Theta=(0, .52, 3), Sz=1)

    start = clock()
    log = trkr.track(deu, ensemble, lattice, 10)
    print("tracking took {} sec".format(clock()-start))


