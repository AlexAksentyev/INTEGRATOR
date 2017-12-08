#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:21:59 2017

@author: alexa

TODO:
    * Tracker should have Particle Logs, a pointer copy to which is given to Ensemble
    at the end of tracking, and deleted

"""
from collections import namedtuple
from time import clock
from importlib import reload

import numpy as NP
from scipy.integrate import odeint

import RHS
reload(RHS)


TrackerControls = namedtuple('TrackerControls', ['inner', 'fwd', 'breaks', 'ncut'])

class Tracker:
    """ Handles tracking of ensembles thru lattices
    """
    def __init__(self, lattice, ensemble):
        self.lattice = lattice
        self.ensemble = ensemble

        self.controls = None
        self._last_pnt = -1 # marker of the state vector upon exiting an element

    def set_controls(self, inner=True, fwd=True, breaks=101, ncut=0):
        self.controls = TrackerControls(inner, fwd, breaks, ncut)

    def _create_logs(self, nturns):
        """ returns first index from which to fill log
        """
        ## p-log data type
        from tables import StringCol
        names = ['START']+[e.fName for e in self.lattice]
        max_len_name = len(names[NP.argmax(names)])
        el_field_type = StringCol(max_len_name) #otherwise problems writing into hdf5 file
        vartype = [('Turn', int), ('Element', el_field_type), ('Point', int)]
        vartype += list(zip(RHS.varname, NP.repeat(float, RHS.varnum)))

        ## the number of records in a p-log
        brks = self.controls.breaks
        el_num = self.lattice.ElCount
        inner = self.controls.inner
        ncut = self.controls.ncut
        if inner:
            nrow = el_num*brks
            if ncut == 0: nrow *= nturns
            else: nrow *= ncut
            ind = 0 # odeint will return injection values
        else:
            nrow = el_num
            if ncut == 0: nrow *= nturns
            else: nrow *= ncut
            nrow += 1 # +1 for injection values
            ind = 1 # odeint won't return injection values; set them manually

        # creating particle logs
        for pid, ic in self.ensemble.ics.items():
            setattr(self.ensemble.Log, 'P'+str(pid), NP.recarray(nrow, dtype=vartype))
            self.ensemble[pid].Log.fill(NP.nan) # in case we later get a NaN error in tracking,
                                       # pre-fill the log with nans
                                       # Turn,Point fill fill with a big random integer
            ic = list(ic.values())
            self.ensemble[pid].Log[0] = 0, names[0], self._last_pnt, *ic # saving injection values
                                                # (will be overwritten if inner is true)

        return ind

    def _setup_file(self, file_handle):
        ## setting up file
        self.file_handle = file_handle
        # write particle parameters (Mass0, KinEn0, G)
        ppar = self.ensemble.Particle.getParams()
        self.file_handle.create_table('/', 'Particle', ppar)
        # separate group to write p-logs in
        self.file_handle.create_group('/', 'Logs', 'Particle logs')

        ## creating data tables to fill
        for p in self.ensemble:
            self.file_handle.create_table(self.file_handle.root.Logs,
                                          'P'+str(p.PID), p.Log.dtype)

    def _write_log(self, from_, to_):
        old_ind, ind = from_, to_
        # write data to file
        for p in self.ensemble:
            tbl = getattr(self.file_handle.root.Logs, 'P'+str(p.PID))
            tbl.append(p.Log[old_ind:ind])
            tbl.flush()

    def track(self, ntimes):
        ## check for controls; if they haven't been set, use defaults
        if getattr(self, 'Controls', None) is None:
            self.setControls()
            print('Setting default controls.')

        lattice = self.lattice
        ## check if there's more than one RF element
        if lattice.RFCount > 1:
            print('\t\t More than one ({}) RF;\n aborting.'.format(lattice.RFCount))
            return

        ind = self._create_logs(ntimes)

        ncut = self.controls.ncut
        cut = True # reset log index after writing to file
        if ncut == 0:
            ncut = NP.floor(ntimes/10) # ncut is used below to decide whether to write data
                        # if we keep all data in RAM, still backup every 10% of turns
                        # 10% is arbitrary
            cut = False

        latname = lattice.Name
        print('Saving data to file {} every {} turns'.format(latname, ncut))

         # current state vector
        ics = list()
        ensemble = self.ensemble
        for ic in ensemble.ics.values(): ics.append(list(ic.values()))
        state = NP.array(ics) # [[x0,y0,...], [x1,y1,...], [x2,y2,...]]
        n_ics = len(state)
        n_var = len(state[0])

        # create the RHS
        rhs = RHS.RHS(ensemble, lattice.getRF()) # setting up the RHS

        # opening hdf5 file to output data
        # write the used particle parameters (Mass0, KinEn0, G)
        # write particle logs
        filename = './data/{}.h5'.format(latname)
        el_seq = lattice.Sequence
        el_num = lattice.ElCount
        import tables as TBL
        with TBL.open_file(filename, 'w', latname) as f:
            self._setup_file(f)
            old_ind = 0 # log is written into file from old_ind:ind
            old_turn = 0 # every ncut turns
            old_percent = -1 # to print 0 %
            cnt = 0; cnt_tot = el_num*ntimes # for progress bar
            print('\t\t LATTICE: {} '.format(latname))
            fwd = self.controls.fwd
            brks = self.controls.breaks
            inner = self.controls.inner
            for n in range(1, ntimes+1): # turn loop
                for i in range(el_num): # element loop
                    # pick element
                    if fwd: element = el_seq[i]
                    else: element = el_seq[len(el_num)-1-i]

                    # print computation progress
                    percent = int(cnt/cnt_tot*100)
                    if percent%10 == 0 and percent != old_percent:
                        print('Complete {} %'.format(percent))
                        old_percent = percent
                    cnt += 1

                    # choose if integrate or otherwise advance the state vector
                    rf_el = element.bSkip

                    #integrate at these points
                    at = NP.linspace(0, element.fLength, brks)

                    try:
                        element.frontKick(state)
                        state = state.reshape(n_ics*n_var) # flat [x0,y0,...,x1,y1,...,x2,y2]
                        if not rf_el:
                            vals = odeint(rhs, state, at, args=(element, brks))
                            state = vals[brks-1] # [x0,y0,...,x1,y1,...]
                        else:
                            element.advance(state)
                        state = state.reshape(n_ics, n_var) # [[x0,y0,...],
                                                            # [x1,y1,...],
                                                            # [x2,y2,...]]
                        element.rearKick(state)
                        if not rf_el and inner:
                            for k in range(brks-1):
                                valsk = vals[k].reshape(n_ics, n_var)
                                for pid in ensemble.ics.keys():
                                    ensemble[pid].Log[ind] = n, element.fName, k, *valsk[pid]
                                ind += 1
                        for pid in ensemble.ics.keys():
                            ensemble[pid].Log[ind] = n, element.fName, self._last_pnt, *state[pid]
                        ind += 1
                    except ValueError:
                        print('NAN error: Element {}, turn {}, log index {}'.format(element.fName, n, ind))
                        return
                # end element loop

                if (n-old_turn)%ncut == 0:
                    print('turn {}, writing data ...'.format(n))
                    start = clock()
                    self._write_log(old_ind, ind)
                    if cut: ind = 0 # if cut, logs can only keep data for one turn => overwrite
                                    # otherwise, keep writing log
                    old_ind = ind
                    old_turn = n
                    print('current index {}'.format(old_ind))
                    print('writing took: {:04.2f} secs'.format(clock()-start))

            # end turn loop

            ## write any remainig data
            print('writing remaining ({}) data ...'.format(ind-old_ind))
            start = clock()
            self._write_log(old_ind, ind)
            print('writing took: {:04.2f} secs'.format(clock()-start))

        print('Complete 100 %')

#%%

if __name__ == '__main__':
    import Ensemble as ENS
    import Particle as PCL
    import Element as ENT

    from BNL import SSb1H2
    E = ENS.Ensemble.populate(PCL.Particle(), x=(-1e-3, 1e-3, 3), dK=(0, 1e-4, 3), Sz=1)

    lat = ENT.Lattice(SSb1H2, 'SSb1H2')
    t = Tracker(lat, E)
    start = clock()
    t.track(10)
    print('time passed {:04.2f}'.format(clock()-start))
