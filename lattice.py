# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 20:08:45 2017

@author: Аксентьев

TODO: 
* rethink the Segment-Lattice paradigm; make Segment a full-fledged class, and Lattice 
as a collection of segments
"""
from copy import deepcopy
import re
import collections as cln
import numpy as np
from element import RF
from utilities import MutableNamedTuple
from plog import PLog

def compose2(f, g):
    """Composition f o g."""
    return lambda *a, **kw: f(g(*a, **kw))

###*** obsolete 
def track(state, transfer_map, acc_len, n_trn, n_rec = None):
    n_trn = int(n_trn)
    n_rec = int(n_rec) if n_rec is not None else n_trn
    n_skip = int(n_trn/n_rec)
    print("N turn {}, N rec {}, N skip {}".format(n_trn, n_rec, n_skip))
    log = PLog(state, n_rec+1) # +1 for initial state

    TM_n = transfer_map**n_skip

    for i in range(1, n_rec+1):
        trn = i*n_skip
        state = TM_n*state # state is matrix
        log[i] = ((trn, -1, trn*acc_len), state.A) # retrieve array

    return log

def track_each(state, map_sequence, n_trn):
    if n_trn > 100:
        print("You don't want that many turns")
        return 
    n_el = len(map_sequence)
    print(n_el)
    log = PLog(state, n_trn*n_el+1) # +1 for initial state

    i = 1; s=0
    for trn in range(n_trn):
        # print("turn {}".format(trn))
        for el in range(n_el):
            element = map_sequence[el]
            # print("element: {} ({})".format(element.name, el))
            state = element.TM*state
            s += element.length
            log[i] = ((trn, el, s), state.A)
            i += 1

    return log
###***

class Segment(MutableNamedTuple):
    __slots__ = ['_EIDs', '_TM', '_count']

    def __init__(self, eid_list, transfer_map=np.eye(6)):
        self._EIDs = eid_list[:]
        self._TM = transfer_map
        self._count = len(self._EIDs)

    @property
    def TM(self):
        return self._TM

    @property
    def count(self):
        return self._count

    @property
    def EIDs(self):
        return self._EIDs

    def shift(self, places):
        self._EIDs = [places + e for e in self._EIDs]

    def merge(self, other):
        if self._TM.__class__ != other._TM.__class__:
            print("Incompatible transfer maps!")
            return
        if isinstance(self._TM, np.ndarray):
            return Segment(self._EIDs+other._EIDs, other._TM*self._TM)
        else:
            return Segment(self._EIDs+other._EIDs, compose2(other._TM, self._TM))

class Lattice:
    def __init__(self, element_sequence, name, segment_map=None):
        self._sequence = deepcopy(element_sequence)
        self._count = len(self._sequence)
        self._length = sum([element.length for element in self._sequence])
        self.name = name
        self._state = 0 # this keeps track of the lattice tilt state

        print("Dealing with segment name: {}".format(self.name))

        ## elements in _sequence are pointers, and hence multiple pointers will point to the same
        ## element. Therefore when a _sequence element it tilted, the underlying element will be tilted
        ## several times, instead of several places in the lattice being misaligned.
        ## to avoid this, deep copy repeating element pointers
        ids = set() # element id's
        rep_cnt = 0 # repeat (element) counter for check
        for index, element in enumerate(self._sequence):
            eid = id(element)
            if eid in ids: # if this id has been encountered before
                self._sequence[index] = deepcopy(element) # replace this position with a new object
                rep_cnt += 1
            ids.add(eid)
        print("Repeat elements: {}".format(rep_cnt))

        ## now deal with the segment map structure
        if segment_map is None:
            ## this is to split the element sequence into segments:
            flag = np.array([element.call for element in self.elements()])
            split_i = np.array(list(range(self.count)))[flag]
            split_i = np.concatenate((split_i, split_i+1))
            split_i.sort()
            segment_split = [seg for seg in np.split(self._sequence, split_i) if len(seg) > 0] # remove emtpy segments
            ## now construct the segment transfer maps
            self._compute_segment_maps(segment_split)
        else:
            self.segment_map = segment_map


    def _compute_segment_maps(self, segment_split):
        self.segment_map = {}
        eid0 = 0
        for cnt, seg in enumerate(segment_split):
            if seg[0].call:
                tm = lambda x: x
                for eid, el in enumerate(seg):
                    tm = compose2(el.TM, tm)
            else:
                tm = np.eye(6)
                for eid, el in enumerate(seg):
                    tm = el.TM*tm
            eid1 = eid0+eid+1
            seg_name = self.name+'_'+str(cnt)
            self.segment_map.update({seg_name: Segment(list(range(eid0, eid1)), tm)})
            eid0 = eid1

    def _segment_split(self):
        """Produces a list of [element: element in segment] for segment in Lattice::Segments."""
        segment_split = []
        for _, seg in self.segments():
            segment_split.append(np.array(self._sequence)[seg.EIDs])

        return segment_split

    @property
    def state(self):
        return self._state

    @property
    def length(self):
        return self._length

    @property
    def count(self):
        return self._count

    def TM(self, segment_name=None):
        if segment_name is not None:
            tm = self.segment_map[segment_name].TM
        else:
            tm = self._transfer_map
            
        return tm

    def merge_segments(self, *names):
        merged_name = '+'.join(names)
        merged = self.segment_map.pop(names[0])
        for name in names[1:]:
            new = self.segment_map.pop(name)
            merged = merged.merge(new)
        self.segment_map.update({merged_name: merged})

    def __call__(self, state, n_trn, n_rec=None):
        print("Lattice.__call__ method:")
        n_trn = int(n_trn)
        n_rec = int(n_rec) if n_rec is not None else n_trn
        n_skip = int(n_trn/n_rec)
        log = PLog(state, n_rec+1) # +1 for initial state

        print("Turn: {}, Rec: {}, Skip: {}".format(n_trn, n_rec, n_skip))

        i = 1; s = 0
        for trn in range(1, n_trn+1):
            for key, segment in self.segments():
                TM = segment.TM
                if isinstance(TM, np.ndarray):
                    state = TM*state
                else:
                    state = TM(state)
            if trn%n_skip == 0:
                log[i] = ((trn, -1, trn*self._length), state.A)
                i += 1

        return log

    def __add__(self, other):
        if not isinstance(other, Lattice):
            print('Cannot add a non-lattice object')
            return # but think about adding an element sequence

        el_sequence = self._sequence + other._sequence
        name = self.name + "+" + other.name
        for _, segment in other.segments():
            segment.shift(self.count)
        segment_map = dict(**self.segment_map, **other.segment_map)

        return Lattice(el_sequence, name, segment_map) # False for now for compatibility with track/track_each

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
        for name, segment in self.segment_map.items():
            yield name, segment

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
        mean_angle : radians
            sets up a systematic angle deviation;
            each tuple element corresponds to an *order* letter;
            if the number of elements in the tuple > than that in the order
            string, the remaining angles are ignored

        sigma : radians
            sets up the rms of the angle distributions
        """

        # pick angles from a the distribution
        angle = np.random.normal(mean_angle, sigma, size=self.count)
        # tilt elements
        for i, element in enumerate(self.elements()):
            element.s_tilt(angle[i])

        segment_split = self._segment_split()
        self._compute_segment_maps(segment_split)

        self._state += 1  # the lattice is in a new state => write into a new group

    def clear(self):
        """Resets the lattice to the original state."""
        self._transfer_matrix = np.eye(6)
        for element in self.elements():
            element.clear_tilt()

        segment_split = self._segment_split()
        self._compute_segment_maps(segment_split)
            
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
    F = MQuad(p, 25e-2, 8.6, 'F')
    D = MQuad(p, 25e-2, -8.11, 'D')
    RF = RF(p, 25e-2*3, 75000)
    
    FODO_0 = Lattice([O,F,O,D], 'FODO')
    FO = Lattice([O, F], 'FO')
    DO = Lattice([O, D], 'DO')
    FODO_1 = FO + DO; FODO_1.name='FODO'
    
    DIP = MDipoleSect(p, 25e-2, 1)
    CWF = CylWien(p, 361e-2, 120e-5, .46)

    ARC = Lattice([DIP, CWF, RF], 'ARC')

    lattice = FODO_0 + ARC
