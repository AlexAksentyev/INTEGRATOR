#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:21:59 2017

@author: alexa

"""
from collections import namedtuple
from time import clock
from importlib import reload
from tables import open_file, StringCol

import numpy as np
from scipy.integrate import odeint

import rhs
reload(rhs)

from particle_log import PLog


TrackerControls = namedtuple('TrackerControls', ['fwd', 'inner', 'breaks', 'ncut'])

class StopTracking(Exception):
    """ Raise upon ValueError in RHS() to stop Tracker from tracking."""
    pass

class Tracker:
    """Handles tracking of ensembles thru lattices."""
    def __init__(self):
        self.controls = None
        self._last_pnt = -1 # marker of the state vector upon exiting an element

        ## created in Tracker::track
        self.lattice = None
        self.ensemble = None
        self.rhs = None # needs to know lattice and ensemble data to set up
        self.file_handle = None # for data backup

    def set_controls(self, fwd=True, inner=False, breaks=101, ncut=0):
        """Choose whether to track a lattice forward/backward (fwd),
        log state variable values inside elements (inner),
        at how many nodes to compute state vars inside an element (breaks),
        data for how many turns to keep in RAM/how often to do file backup (ncut).
        If ncut == 0, keep all data in RAM; in this case I still back up every
        10 % of the total number of turns (see Tracker::track).
        """
        self.controls = TrackerControls(fwd, inner, breaks, ncut)

    def _create_logs(self, nturns):
        """Returns the first index from which to fill the log.
        A filled log will have a number of nan values
        if the lattice contains an RF element with bool_skip = True
        and the control inner = True. That's normal; it's because when counting 
        the number nrow at case inner, we assume ALL elements will return 
        a brks number of states, whereas in _run_turn we explicitly say that if 
        bool_skip = True, only the last state is written to log.
        Hence we preallocate a brks*nturns more rows than necessary.
        """

        ## the number of records in a p-log
        brks = self.controls.breaks
        el_num = self.lattice.count
        inner = self.controls.inner
        ncut = self.controls.ncut
        if inner:
            nrow = el_num*brks
            if ncut == 0:
                nrow *= nturns
            else:
                nrow *= ncut
            ind = 0 # odeint will return injection values
        else:
            nrow = el_num
            if ncut == 0:
                nrow *= nturns
            else:
                nrow *= ncut
            nrow += 1 # +1 for injection values
            ind = 1 # odeint won't return injection values; set them manually

        self.log = PLog(self.ensemble, nrow)

        return ind

    def _setup_file(self, file_handle):
        ## setting up file
        self.file_handle = file_handle
        # write particle parameters (Mass0, KinEn0, G)
        ppar = self.ensemble.particle.get_params()
        self.file_handle.create_table('/', 'Particle', ppar)
        # separate group to write p-logs in
        self.file_handle.create_group('/', 'Logs', 'Particle logs')

        ## creating data tables to fill
        for p in self.ensemble:
            self.file_handle.create_table(self.file_handle.root.Logs,
                                          'P'+str(p.pid), self.log.dtype)

    def _run_turn(self, current_turn, log_index, state):
        brks = self.controls.breaks

        el_num = self.lattice.count
        el_seq = self.lattice._sequence

        n_ics = self.rhs.n_ics
        n_var = rhs.VAR_NUM

        for i in range(el_num): # element loop
            # pick element
            if self.controls.fwd: element = el_seq[i]
            else: element = el_seq[len(el_num)-1-i]

            # choose if integrate or otherwise advance the state vector
            rf_el = element.bool_skip

            #integrate at these points
            at = np.linspace(0, element.length, brks)

            try:
                element.front_kick(state)
                state = state.reshape(n_ics*n_var) # flat [x0,y0,...,x1,y1,...,x2,y2]
                if not rf_el:
                    vals = odeint(self.rhs, state, at, args=(element, brks))
                    state = vals[brks-1] # [x0,y0,...,x1,y1,...]
                else:
                    element.advance(state)
                state = state.reshape(n_ics, n_var) # [[x0,y0,...],
                                                    # [x1,y1,...],
                                                    # [x2,y2,...]]
                element.rear_kick(state)
                if not rf_el and self.controls.inner:
                    for k in range(brks-1):
                        self.log[log_index] = ((current_turn, element.name, i, k), vals[k])
                        log_index += 1
                self.log[log_index] = ((current_turn, element.name, i, self._last_pnt), state.flatten())
                log_index += 1
            except ValueError:
                print('NAN error: Element {}, turn {}, log index {}'.format(element.name, current_turn, log_index))
                raise StopTracking
        # end element loop

        return state, log_index

    def _write_log(self, from_, to_):
        old_ind, ind = from_, to_
        # write data to file
        for p in self.ensemble:
            tbl = getattr(self.file_handle.root.Logs, 'P'+str(p.pid))
            tbl.append(self.log[old_ind:ind, p.pid])
            tbl.flush()

    def track(self, ensemble, lattice, n_turns):
        """Track ensemble through lattice for n_turns.
        Returns the particle log of type PLog.
        """

        self.ensemble = ensemble
        self.lattice = lattice
        ## check for controls; if they haven't been set, use defaults
        if getattr(self, 'controls', None) is None:
            self.set_controls()
            print('Setting default controls.')

        ## check if there's more than one RF element
        if lattice.RF.count > 1:
            print('\t\t More than one ({}) RF;\n aborting.'.format(lattice.RF.count))
            return

        log_ind = self._create_logs(n_turns)

        ncut = self.controls.ncut
        cut = True # reset log index after writing to file
        if ncut == 0:
            ncut = np.floor(n_turns/10) # ncut is used below to decide whether to write data
                        # if we keep all data in RAM, still backup every 10% of turns
                        # 10% is arbitrary
            if ncut < 101: # unless it's less than  100,
                ncut = 100 # if so, backup every 100 turns
            cut = False

        latname = lattice.name
        print('Saving data to file {} every {} turns'.format(latname, ncut))

         # initial state vector
        ics = list()
        for ic in ensemble.ics.values(): ics.append(list(ic.values()))
        state = np.array(ics) # [[x0,y0,...], [x1,y1,...], [x2,y2,...]]

        # create the RHS
        self.rhs = rhs.RHS(ensemble, lattice.getRF()) # setting up the RHS

        # opening hdf5 file to output data
        # write the used particle parameters (Mass0, KinEn0, G)
        # write particle logs
        filename = './data/{}.h5'.format(latname)
        with open_file(filename, 'w', latname) as file_handle:
            self._setup_file(file_handle)
            old_ind = 0 # log is written into file from old_ind:log_ind
            old_turn = 0 # every ncut turns
            old_percent = -1 # to print 0 %
            print('\t\t LATTICE: {} '.format(latname))
            for turn in range(1, n_turns+1): # turn loop

                # print computation progress
                percent = int((turn-1)/n_turns*100)
                if percent%10 == 0 and percent != old_percent:
                    print('Complete {} %'.format(percent))
                    old_percent = percent

                try:
                    state, log_ind = self._run_turn(turn, log_ind, state)
                except StopTracking:
                    return

                if (turn-old_turn)%ncut == 0:
                    print('turn {}, writing data ...'.format(turn))
                    start = clock()
                    self._write_log(old_ind, log_ind)
                    if cut: log_ind = 0 # if cut, logs can only keep data for one turn => overwrite
                                    # otherwise, keep writing log
                    old_ind = log_ind
                    old_turn = turn
                    print('current index {}'.format(old_ind))
                    print('writing took: {:04.2f} secs'.format(clock()-start))

            # end turn loop

            ## write any remainig data
            print('writing remaining ({}) data ...'.format(log_ind-old_ind))
            start = clock()
            self._write_log(old_ind, log_ind)
            print('writing took: {:04.2f} secs'.format(clock()-start))

        print('Complete 100 %')
        
        return self.log
        del self.log

#%%

if __name__ == '__main__':
    from ensemble import Ensemble
    from particle import Particle
    from element import MQuad, Drift
    from lattice import Lattice
    import copy

#    from BNL import BDA
    t = Tracker()

    #%%
    E = Ensemble.populate(Particle(), x=(-5e-3, 5e-3, 3), dK=(0, 1e-4, 3), Sz=1)
    E_tilted = copy.deepcopy(E)

    mqf = MQuad(5e-2, 8.6, 'QF')
    mqd = MQuad(5e-2, -8.3, 'QD')
    dft = Drift(25e-2)
    dft2 = Drift(25e-2)
    FODO = [mqf, dft, mqd, dft]
    MAGNET = [dft, dft, dft]

    LAT = Lattice(FODO, 'test')
    LAT_tilted = copy.deepcopy(LAT)
    LAT.insertRF(0, 0, E, E_field=15e7)
    LAT_tilted.insertRF(0, 0, E_tilted, E_field=15e7)
    
    LAT_tilted.tilt('s', 0, .003)

    #%%
    start = clock()
    log = t.track(E_tilted, LAT_tilted, 10)
    print('time passed {:04.2f}'.format(clock()-start))
    #%%
    E_tilted.plot('Sx', 'x')
