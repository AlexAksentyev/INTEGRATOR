# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 20:08:45 2017

@author: Аксентьев
"""
import copy
import re
import collections as cln
import numpy as np
from element import RF as ERF#, Tilt, Shift
from utilities import MutableNamedTuple

class RF(MutableNamedTuple):
    """Container class to keep information about the lattice's RF element.
    """
    __slots__ = ['index', 'count']

class Lattice:
    def __init__(self, element_sequence, name, segment_map=None):
        sequence = copy.deepcopy(element_sequence)
        ## ensuring the lattice doesn't have more than one RF element
        self.RF = RF(None, 0)
        self.length = 0
        remove_indices = []
        for ind, element in enumerate(sequence):
            if isinstance(element, ERF):
                print('Found RF at {}'.format(ind))
                if self.RF.count < 1:
                    self.RF.count += 1
                    self.RF.index = ind
                else:
                    print('Second RF element encountered; not adding to lattice')
                    remove_indices.append(ind)
                    continue
            self.length += element.length
        remove_indices.reverse() # so that when i start popping, RF elements don't shift places
        for ind in remove_indices:
            sequence.pop(ind)

        self.count = len(sequence)

        ## creating _sequence of elements from pointers
        self._sequence = sequence
        self._transfer_matrix = np.eye(6) ## ALSO, construct the segment transfer matrix
        ids = set()
        for index, element in enumerate(self._sequence):
            self._transfer_matrix = element.M*self._transfer_matrix
            eid = id(element)
            if eid in ids:
                self._sequence[index] = copy.deepcopy(element)
            ids.add(eid)
            
        ##

        self.name = name
        if segment_map is None:
            self.segment_map = {name: list(range(self.count))}
        else:
            self.segment_map = segment_map

        self._state = 0

    def unique_ids(self):
        ids = set()
        for element in self.elements():
            eid = id(element)
            if eid in ids:
                continue
            ids.add(eid)
        return ids

    @property
    def state(self):
        return self._state

    @property
    def M(self):
        return self._transfer_matrix

    def __add__(self, other):
        if not isinstance(other, Lattice):
            print('Cannot add a non-lattice object')
            return # but think about adding an element sequence

        # remove RF's of both
        rf0 = self.pop(self.RF.index)
        if rf0 is not None:
            self.RF.count -= 1
            self.RF.index = None
            print('Achtung! segment {} was updated.'.format(self.name))
        rf1 = other.pop(other.RF.index)
        if rf1 is not None:
            other.RF.count -= 1
            other.RF.index = None
            print('Achtung! segment {} was updated.'.format(other.name))
        print('New lattice w/o RF elements.')


        # make a sections map:
        # (segment_name, {indices in lattice})
        seg0_num = self.count
        seg1_num = other.count
        segment_map = copy.deepcopy(self.segment_map)
        segmap1 = other.segment_map
        segmap1.update({other.name: list(range(seg0_num, seg0_num + seg1_num))})
        segment_map.update(segmap1)

        # make the common sequence, compound lattice name
        el_sequence = self._sequence + other._sequence
        name = self.name + '+' + other.name

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

    def insert_RF(self, position, length, reference_particle, **ERF_pars):
        """
        Arguments:
            position : integer
                the index inside the lattice element sequence (starts from 0)

            length : float
                at length 0 the RF is simulated by an energy kick;
                internally, for the computation of the potential energy jump,
                assumes a default length of 5e-4, i.e. V = E_field * 5e-4.

        Keyword arguments (ERF_pars):
            E_field : float
                the amplitude of the electric field wave (default 15e5)

            phase : float
                the time-domain phase of the reference particle (default 1.5*pi)

            H_number : integer
                the harmonic huimber (default 50)

            name : string
                defaults to "RF"

        """
        if self.RF.index is not None and self.RF.index != position:
            print("""Trying to add a second RF element;
                  current RF position is {}""".format(self.RF.index))
            return
        if position == self.RF.index:
            print('Replacing RF {}'.format(self.pop(position)))
            self.RF.count -= 1

        # creating an RF element
        full_acc_len = self.length + length
        rf = ERF(reference_particle, length, full_acc_len, **ERF_pars)

        self.insert(position, rf)
        self.RF.index = position
        self.RF.count += 1

    def get_RF(self):
        try:
            return self[self.RF.index]
        except TypeError:
            return None

    def _update_segment_map(self, index, element, action='add'):
        action = action.lower()
        if action[:3] == 'add':
            action = 1
            self.segment_map.update({element.name: [index]})
        elif action[:3] == 'sub':
            action = -1
            try:
                self.segment_map.pop(element.name)
            except KeyError:
                pass
        else:
            raise Exception('Unknown action')

        for name, map_ in self.segment_map.items():
            if name == element.name:
                continue
            self.segment_map[name] = [j for j in map_ if j < index] \
                + [j + action for j in map_ if j >= index]


    ## these two will require recomputation fo the transfer map i think
    ## better yet, make the segment transfer maps items of Lattice.segment_map
    def insert(self, index, element):
        # inserting the element
        self._sequence.insert(index, element)
        self.count += 1
        self.length += element.length
        # updating segment map
        self._update_segment_map(index, element, 'add')

    def pop(self, ind):
        if ind is None:
            print('Passed index is None; no element popped.')
            return None

        try:
            # popping the element
            element = self._sequence.pop(ind)
            self.count -= 1
            self.length -= element.length
            if isinstance(element, ERF):
                self.RF.count -= 1
                self.RF.index = None
            # updating segment map
            self._update_segment_map(ind, element, 'subtract')
        except IndexError:
            print('Wrong index; no element popped.')
            return None

        return element

    def make_segment(self, name):
        """Groups all elements with *name* into a separate segment.
        Doesn't remove the elements from other segments.

        Arguments:
            name: string
                for the present, the exact name of the element; i.e, names
                in the form name_index won't be understood.
        """
        # find all elements with name
        element_names = self.list_names(True)
        element_indices = [i for i, x in enumerate(element_names) if x == name]
        self.segment_map.update({name: element_indices})


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
