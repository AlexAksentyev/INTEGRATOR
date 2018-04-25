#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 16:16:39 2018

@author: alexa
"""
## importing modules
## defining constant objects

from particle import Particle
from tracker import Tracker
from particle_log import StateList
from lattice import Lattice
import element as ent # import elements
import matplotlib.pyplot as plt

trkr = Tracker()
particle = Particle()   # by default sets up a deuteron


#%%
## definition of the initial state ensemble

# in order to define an ensemble of particles do
ensemble = StateList(x=(-1e-3, 1e-3, 5), dK=(0, 1e-4, 5), Theta=(0, .52, 3), Sz=1)
    # the argument format: variable (from, to, number of intermediate values)
    # to set variable value shared accross all states, specify a single value:
    # variable = shared value


#%%
## creation of lattice elements

element_0 = ent.MDipole(25e-2, particle, B_field=1)
element_1 = ent.MQuad(25e-2, 10.11)


#%%
## setting up a lattice object

lattice = Lattice([element_0, element_1], name='Example_lattice')
                                        # the name argument is used in Tracker.track
                                        # to set up a file to output tracking data into

#%%
## tracking

log = trkr.track(particle, ensemble, lattice, 100)

log.plot('Sx', 's', pids=[1,11,65])
plt.title('absolute values plot')
log.plot('-D Sx', 's', pids=[1,11,65])
plt.title('relative values plot')
#%%
## setting up an RF element inside the lattice

lattice.insert_RF(0, 0, particle)

log.plot(pids=[3,17,29,67]) # by default, log.plot() plots the fish plot

#%%
## tilting lattice elements

lattice.tilt('sxy', (.1,-.01, 0), (.3, .05, .02))

# prints the order and angles of tilt for the 0-th lattice element to the console
print(lattice[0].tilt_)

log_tilted = trkr.track(particle, ensemble, lattice, 100)

log_tilted.plot('x', 's', pids=[5, 11, 17])
plt.title('absolute values plot')
log_tilted.plot('-D x', 's', pids=[5, 11, 17])
plt.title('relative values plot')
