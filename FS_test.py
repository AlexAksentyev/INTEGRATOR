#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 11:58:17 2018

@author: alexa
"""

from particle import Particle, CLIGHT
import element as ent
from tracker import Tracker
from lattice import Lattice
from particle_log import StateList



deu = Particle()
bunch = StateList(Sz=1, x=(-1e-3, 1e-3, 3), dK=(0, 1e-4, 3))

gamma, beta = deu.GammaBeta()
B0 = 1
E0 = CLIGHT*beta*gamma**2*deu.G*B0/(beta**2*gamma**2*deu.G - 1)
element = ent.Wien(25e-2, 5e-2, deu, E0, B0)
element = ent.MDipole(25e-2, deu, B_field=B0)

trkr = Tracker()
lattice = Lattice([element], 'test_FS')
lattice.insert_RF(0, 5e-4, deu, E_field=15e7)

log = trkr.track(deu, bunch, lattice, int(1e3))
log.plot()
