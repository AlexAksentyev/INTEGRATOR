# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 20:08:45 2017

@author: Аксентьев
"""
import copy
import re
import collections as cln
import numpy as np
from element import ERF
from utilities import MutableNamedTuple

class RF(MutableNamedTuple):
    """Container class to keep information about the lattice's RF element
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

        self._sequence = sequence
        self.count = len(sequence)
        self.name = name
        if segment_map is None:
            self.segment_map = {name: list(range(self.count))}
        else:
            self.segment_map = segment_map

    def __add__(self, other):
        if not isinstance(other, Lattice):
            print('Cannot add a non-lattice object')
            return # but think about adding an element sequence

        # remove RF's of both
        rf0 = self.pop(self.RF.index)
        if rf0 is not None:
            self.RF.count -= 1
            self.RF.index = None
        rf1 = other.pop(other.RF.index)
        if rf1 is not None:
            other.RF.count -= 1
            other.RF.index = None
        print('New lattice w/o RF elements.')
        print('Achtung! the segments {}, {} were updated.'.format(self.name, other.name))

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
        name = self.name + '_' + other.name

        return Lattice(el_sequence, name, segment_map)

    def __getitem__(self, idx):
        return self._sequence[idx]

    def __repr__(self):
        return self._sequence.__repr__()

    def __iter__(self):
        self.__current_id = 0
        return self

    def __next__(self):
        last_id = self.count - 1
        if self.__current_id <= last_id:
            result = self[self.__current_id]
            self.__current_id += 1
            return result
        else:
            raise StopIteration

    def insertRF(self, position, length, ensemble, **ERF_pars):
        if self.RF.index is not None and self.RF.index != position:
            print("""Trying to add a second RF element;
                  current RF position is {}""".format(self.RF.index))
            return
        if position == self.RF.index:
            print('Replacing RF {}'.format(self.pop(position)))
            self.RF.count -= 1

        # creating an RF element
        full_acc_len = self.length + length
        rf = ERF(length, ensemble, full_acc_len, **ERF_pars)

        self.insert(position, rf)
        self.RF.index = position
        self.RF.count += 1

    def getRF(self):
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

    def insert(self, index, element):
        # inserting the element
        self._sequence.insert(index, element)
        self.count += 1
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
            if isinstance(element, ERF):
                self.RF.count -= 1
                self.RF.index = None
            # updating segment map
            self._update_segment_map(ind, element, 'subtract')
        except IndexError:
            print('Wrong index; no element popped.')
            return None

        return element

    def list_names(self, full=False):
        names = [e.name for e in self]
        if full:
            return names
        return np.unique([re.sub('_.*', '', e) for e in names])

    def tilt(self, order='S', mean_angle=(0,), sigma=(0,), append=False):
        n = len(order)
        if not isinstance(mean_angle, cln.Sequence): mean_angle = (mean_angle,)
        if not isinstance(sigma, cln.Sequence): sigma = (sigma,)
        nmean = len(mean_angle)
        nsig = len(sigma)

        try:
            angle = np.random.normal(mean_angle, sigma, size=(self.count, n))
        except ValueError:
            print('Dimension mismatch: order {}, mean {}, sigma {}'.format(n, nmean, nsig))
            return

        nuids = len(np.unique([id(e.tilt_) for e in self]))
        if nuids != self.count:
            print('Non-unique elements ({}/{}) in lattice. Smart tilting.'.format(self.count-nuids, self.count))
            print('\t Not tilting:')

        i = 0
        ids = set()
        cnt = 1
        for element in self:
            eid = id(element.tilt_)
            if eid in ids:
                print('\t\t {} element {} at lattice index {}'.format(cnt, element.name, i))
                cnt += 1
                i += 1
                continue
            element.tilt(order, *angle[i], append=append)
            ids.add(eid)
            i += 1
            
    def plot_segment(self, segment_name, log, Ylab='-D dK', Xlab='-D Theta', **kwargs):
        try:
            eids = self.segment_map[segment_name]
        except KeyError:
            print('Wrong segment name.')
            return
        
        ii = [eid in eids for eid in log[:,0]['EID']]
        segment_log = log[ii]
        
        segment_log.plot(Ylab=Ylab, Xlab=Xlab, **kwargs)
        
#%%
if __name__ == '__main__':
    from ensemble import Ensemble
    from particle import Particle
    from BNL import SSb1H2, SSb1H1

    bunch = Ensemble.populate(Particle(), Sz=1, x=(-1e-3, 1e-3, 3))

    RF0 = ERF(5, bunch, 5)
    RF1 = ERF(5, bunch, 5)

    SSb1H2.insert(5, RF1)
    SSb1H2.insert(0, RF0)

    segment_0 = Lattice(SSb1H1, 'SSb1H1')
    segment_1 = Lattice(SSb1H2, 'SSb1H2')
    section = segment_0 + segment_1
    section.insertRF(section.segment_map['SSb1H2'][0], 0, bunch, E_field=15e7)

    from tracker import Tracker
    trkr = Tracker()
    trkr.set_controls(inner=False,breaks=3)
    bunch.log = trkr.track(bunch, section, 10)

    bunch.plot('Sx', 's')
