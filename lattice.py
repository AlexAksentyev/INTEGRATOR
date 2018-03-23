# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 20:08:45 2017

@author: Аксентьев
"""
import copy
import re
import collections as cln
import numpy as np
from element import ERF, Tilt, Shift
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
        ids = set()
        for index, element in enumerate(self._sequence):
            eid = id(element)
            if eid in ids:
                #print('Repeated element; copying...')
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

    def tilt(self, order='S', mean_angle=(0,), sigma=(0,), append=False):
        """
        Arguments
        ________________

        order : string
            order of axises (x, y, s; case-insensitive) about which
            rotations are performed, as in: 'xsxy', meaning
            "rotate about the x axis, then s, then x again, then y"

        mean_angle : tuple
            sets up a systematic angle deviation;
            each tuple element corresponds to an *order* letter;
            if the number of elements in the tuple > than that in the order
            string, the remaining angles are ignored

        sigma : tuple
            sets up the rms of the angle distributions

        append : boolean
            if the lattice has already been rotated, append=True adds new
            rotations on top of the old ones

        """
        # prepping tilt sigma, mean
        tilt_num = len(order)
        if not isinstance(mean_angle, cln.Sequence): mean_angle = (mean_angle,)
        if not isinstance(sigma, cln.Sequence): sigma = (sigma,)
        nmean = len(mean_angle)
        nsig = len(sigma)

        # make tilt angles
        try:
            angle = np.random.normal(mean_angle, sigma, size=(self.count, tilt_num))
        except ValueError:
            print('Dimension mismatch: order {}, mean {}, sigma {}'.format(tilt_num, nmean, nsig))
            return

        # tilting proper
        # i = 0
        for i, element in enumerate(self.elements()):
            element.tilt(order, *angle[i], append=append)
            # i += 1

        self._state += 1  # the lattice is in a new state => write into a new group

    def shift(self, x=(0,0), y=(0,0), append=False):
        """
        Arguments
        ________________
        
        x, y : pair (mean shift, sigma shift)

        append : boolean
            if the lattice has already been rotated, append=True adds new
            shifts on top of the old ones

        """
        mu, sigma = x
        x_shift = np.random.normal(mu, sigma, self.count)
        mu, sigma = y
        y_shift = np.random.normal(mu, sigma, self.count)

        for i, element in enumerate(self.elements()):
            element.shift(x_shift[i], y_shift[i], append=append)

        self._state += 1

    def clear(self):
        """Resets the lattice to the original state."""
        for element in self.elements():
            element.tilt_ = Tilt()
            element.shift_ = Shift()
            
        self._state = 0

        

    def plot_segment(self, segment_name, log,
                     Ylab='-D dK', Xlab='-D Theta', **kwargs):
        try:
            eids = self.segment_map[segment_name]
        except KeyError:
            print('Wrong segment name.')
            return

        ii = [eid in eids for eid in log[:, 0]['EID']]
        segment_log = log[ii]

        segment_log.plot(Ylab=Ylab, Xlab=Xlab, **kwargs)

    def segment_edges(self, log):
        s_list = list()
        for ind, key_val in enumerate(self.segment_map.items()):
            name, seg_ids = key_val
            i0 = [eid in seg_ids for eid in log[:, 0]['EID']].index(True)
            ii = np.arange(i0, len(log), self.count)
            s_list += log[ii, 0]['s'].tolist()

        return s_list

#%%
if __name__ == '__main__':
    from particle_log import StateList
    from particle import Particle
    import BNL
    from matplotlib import pyplot as plt
    from tracker import Tracker

    import element as ent

    trkr = Tracker()
    trkr.set_controls(inner=False, breaks=3)

    #%%
    ini_states = StateList(Sz=1, x=1e-3)
    deuteron = Particle()

    segment_0 = Lattice(BNL.SSb1H2, 'SS')
    segment_1 = Lattice(BNL.ARCb2H2, 'ARC')
    segment_2 = Lattice(BNL.SSe1H1, 'SS')

    #%%

    section = segment_0 + segment_1 + segment_2
    section.insert_RF(0, 0, deuteron, E_field=15e7)

#%%
    R3 = ent.Wien(361.55403e-2, 5e-2, deuteron, -120e5, .082439761*0, name="R3")

    section = Lattice([R3], 'Wien')

    log = trkr.track(deuteron, ini_states, section, 1000)

    #%%
#    section.plot_segment('RF', log, 'Sx', 's')
    log.plot('x', 's', pids=[0])
#    log.plot('Sx', 's', new_plot=False)
#    sec_edges = section.segment_edges(log)
#    for edge in sec_edges:
#        plt.axvline(x=edge, linewidth=.8, color='b')
#    plt.title(section.name)
