# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 20:08:45 2017

@author: Аксентьев
"""
from copy import deepcopy
import re
import collections as cln
import numpy as np
from element import RF
from utilities import MutableNamedTuple

class Segment(MutableNamedTuple):
    __slots__ = ['_EIDs', '_TM']

    def __init__(self, eid_list, transfer_matrix=np.eye(6)):
        self._EIDs = eid_list[:]
        self._TM = transfer_matrix

    @property
    def TM(self):
        return self._TM
    @property
    def EIDs(self):
        return self._EIDs

    def shift(self, places):
        self._EIDs = [places + e for e in self._EIDs]

class Lattice:
    def __init__(self, element_sequence, name, segment_map=None):
        ## creating _sequence of elements from pointers
        self._sequence = deepcopy(element_sequence)
        self._transfer_matrix = np.eye(6) ## ALSO, construct the segment transfer matrix
        self._length = 0
        self._count = 0
        ids = set() # element id's
        for index, element in enumerate(self._sequence):
            self._transfer_matrix = element.M*self._transfer_matrix
            self._length += element.length
            self._count += 1
            eid = id(element)
            if eid in ids: # if this id has been encountered before
                self._sequence[index] = deepcopy(element) # replace this position with a new object
            ids.add(eid)
        ##
        self.name = name
        if segment_map is None:
            self.segment_map = {self.name: Segment(list(range(self._count)), self._transfer_matrix)}
        else:
            self.segment_map = segment_map

        self._state = 0

    @property
    def state(self):
        return self._state

    @property
    def length(self):
        return self._length

    @property
    def count(self):
        return self._count

    def __add__(self, other):
        if not isinstance(other, Lattice):
            print('Cannot add a non-lattice object')
            return # but think about adding an element sequence

        el_sequence = self._sequence + other._sequence
        name = self.name + "+" + other.name
        for segment in other.segments():
            segment.shift(self.count)
        segment_map = dict(**self.segment_map, **{other.name: other.segment_map})

        return Lattice(el_sequence, name, segment_map)

    def __getitem__(self, idx):
        return self._sequence[idx]

    def __repr__(self):
        return self._sequence.__repr__()

    def elements(self, from_=0, to_=None):
        """Generator for iterating through lattice elements.
        """
        if to_ is None or to_ > self.count:
            to_ = self.count
        for eid in range(from_, to_):
            yield self[eid]

    def segments(self):
        for segment in self.segment_map.values():
            yield segment

    def list_names(self, full=False):
        names = [e.name for e in self]
        if full:
            return names
        return np.unique([re.sub('_.*', '', e) for e in names])

    def s_tilt(self, mean_angle=0, sigma=0):
        """Tilts lattice elements about the s-axis. The tilt angles are taken from 
        a normal distribution with mean mean_angle and standard deviation sigma.

        Arguments
        ________________
        mean_angle : float
            sets up a systematic angle deviation;
            each tuple element corresponds to an *order* letter;
            if the number of elements in the tuple > than that in the order
            string, the remaining angles are ignored

        sigma : float
            sets up the rms of the angle distributions
        """

        # pick angles from a the distribution
        angle = np.random.normal(mean_angle, sigma, size=self.count)
        
        # tilting proper
        ## updates the transfer matrix
        self._transfer_matrix = np.eye(6)
        for i, element in enumerate(self.elements()):
            element.s_tilt(angle[i])
            self._transfer_matrix = element.M*self._transfer_matrix

        self._state += 1  # the lattice is in a new state => write into a new group

    def clear(self):
        """Resets the lattice to the original state."""
        self._transfer_matrix = np.eye(6)
        for element in self.elements():
            element.clear_tilt()
            self._transfer_matrix = element.M*self._transfer_matrix
            
        self._state = 0


if __name__ == '__main__':
    from particle import Particle, EZERO, CLIGHT
    from plog import PLog
    from element import Drift, MQuad, RF, MDipoleSect, CylWien
    from lattice import Lattice
    from state_list import StateList
    from time import clock
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    p = Particle()
    O = Drift(p, 25e-2)
    F = MQuad(p, 25e-2, 8.6)
    D = MQuad(p, 25e-2, -8.11)
    RF = RF(p, 25e-2*3, 75000)
    
    FODO = Lattice([O,F,O,D,O,RF], 'FODO')
    
    DIP = MDipoleSect(p, 25e-2, 1)
    
    CWF = CylWien(p, 361e-2, 120e-5, .46)

    ARC = Lattice([DIP, CWF], 'ARC')

    lattice = FODO + ARC
