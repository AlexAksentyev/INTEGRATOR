#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 09:01:42 2017

@author: alexa
"""
import numpy as np
import tables as tbl
import rhs
import copy
import particle as pcl

from utilities import Bundle

#%%
class StateList:
    """Create an ensemble of initial conditions.
    """
    def __init__(self, **kwargs):

        keys = kwargs.keys()

        # create defined variables
        ntot = 1
        arg_dict = dict()
        for key, val in kwargs.items():
            try:
                lb, ub, num = val
                val = np.linspace(lb, ub, num)
            except TypeError: # if key = value, set value for all pcls
                num = 1
            ntot *= num
            arg_dict.update({rhs.IMAP[key]: val})

        # make mesh
        mesh = dict(zip(keys, np.meshgrid(*list(arg_dict.values()))))

        vartype = list(zip(rhs.VAR_NAME, np.repeat(float, rhs.VAR_NUM)))
        self.state_list = np.zeros(ntot+1, dtype=vartype) # +1 for genuine refrence particle

        #write data
        for key, value in mesh.items():
            self.state_list[key] = np.array([0]+value.reshape(ntot).tolist())

        self.state_list[0]['Sz'] = 1

        #here i should remove any duplicate reference particles
#        self.state_list = np.unique(self.state_list) # this works but messes up the order

        # convert to list of dicts for use with ensemble
        self.state_list = [dict(zip(self.state_list.dtype.names, x)) for x in self.state_list]

    @classmethod
    def from_list(cls, state_list):
        if not isinstance(state_list, list):
            print('Need a list!')
            return
        if len(state_list[0]) != rhs.VAR_NUM:
            print('Wrong state vector length')
            return
        if isinstance(state_list[0], list):
            state_list = [dict(zip(rhs.VAR_NAME, x)) for x in state_list]
            print('Converted to a list of dicts.')

        result = cls()
        result.state_list = state_list

        return result

    def __len__(self):
        return len(self.state_list)

    def __getitem__(self, pid):
        return self.state_list[pid]

    def __repr__(self):
        from pandas import DataFrame
        return str(DataFrame(self.state_list))

    def as_list(self):
        states = list()
        for d in self.state_list:
            states.append(list(d.values()))
        return states

#%%


#%%


class PLog(np.recarray):
    """Particle log is a 3D recarray, where the 1st index refers to the particle (ics),
    the 2nd to the row (record), the 3rd to the variable (like x,y,s, &c).
    """

    max_len_name = 10
    el_field_type = tbl.StringCol(max_len_name) # otherwise problems writing into hdf5 file
    metadata_type = [('Turn', int), ('Element', el_field_type),
                     ('EID', int), ('Point', int)]
    variable_type = list(zip(rhs.VAR_NAME, np.repeat(float, rhs.VAR_NUM)))
    record_type = metadata_type + [('PID', int)] + variable_type
                ## PID field is weird like this to enable automatic
                ## pid setting in setitem (if it were metadata, i'd have to supply it)
                ## i could actually do without it, now that i know that subsetting
                ## a recarray makes it 1D, and I have to manually reshape it.

    _state_var_1st_ind = len(metadata_type) + 1 # len(metadata) + len(PID) - 1 + 1
    last_pnt_marker = -1 # marker of the state vector upon exiting an element

    def __new__(cls, ini_states, particle, n_records):

        if not len(ini_states[0]) == rhs.VAR_NUM:
            print('Wrong number of state variables')
            return
        if not isinstance(ini_states[0], dict):
            ini_states = [dict(zip(rhs.VAR_NAME, state)) for state in ini_states]


        n_ics = len(ini_states)

        obj = np.recarray((n_records, n_ics), dtype=cls.record_type, order='C').view(cls)
        super(PLog, obj).fill(np.nan)
        obj.n_ics = n_ics
        obj.ics = ini_states
        obj.particle = particle

        ics = list()
        for ic in ini_states:
            ics.append(list(ic.values()))

        ics = np.array(ics).flatten()

        obj[0] = ((0, 'START', -1, cls.last_pnt_marker), ics)

        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

        self.particle = getattr(obj, 'particle', None)
        self.n_ics = getattr(obj, 'n_ics', None)
        self.ics = getattr(obj, 'ics', None)

    @classmethod
    def from_file(cls, filename, directory='./data/'):
        """Returns the multi-D PLog, and a dictionary of initial conditions
        for initializing a corresponding ensemble
        """


        filename = '{}{}.h5'.format(directory, filename)
        with tbl.open_file(filename) as file_handle:
            try:
                log_tables = file_handle.root.logs
            except tbl.NoSuchNodeError:
                print('File contains no logs')
                return

            npids = len(file_handle.list_nodes('/logs'))

            tbl0 = file_handle.root.logs.P0
            nrow = tbl0.nrows
            rec_type = tbl0.dtype

            Log = np.recarray((nrow, npids), dtype=rec_type)

            ind_x = cls._state_var_1st_ind
            ics = list() # ics and particle data are required to build _host
            for ind, log in enumerate(log_tables):
                Log[:, ind] = log.read()
                ics.append(list(Log[0, ind])[ind_x:])

            ## particle data
            particle = file_handle.root.particle[0]
            particle = pcl.Particle(particle['Mass0'], particle['KinEn0'], particle['G'])

        Log = Log.view(cls)
        Log.particle = particle
        Log.ics = StateList.from_list(ics)
        Log.n_ics = len(ics)

        return Log

    def write_file(self, file_handle, from_=0, to_=None):
        if to_ is None: to_ = self.nrow - 1
        pids = self[0]['PID']
        for pid in pids:
            tbl = getattr(file_handle.root.logs, 'P'+str(pid))
            tbl.append(self[from_: to_, pid])
            tbl.flush()

    def __setitem__(self, i, stamp_vector):
        """Scipy integrator returns a flat array of state vector values
        [x_i^0, y_i^0, ..., Sz_i^0, x_i^1, y_i^1, ...], where the superscript
        denotes the particle number.
        We want to append each particle's state variable vector with some
        metadata, and turn the resulting vector into a tuple.
        """
        stamp, vector = stamp_vector
        vector = vector.reshape(self.n_ics, -1)

        result = np.empty(self.n_ics, self.dtype)

        for ind, vec in enumerate(vector):
            result[ind] = stamp + (ind,) + tuple(vec)

        super(PLog, self).__setitem__(i, result)

    def __bundle_up(self, pid):
        log = self[:, pid]
        current_state = log[len(log)-1]
        return Bundle(pid=pid,
                      ics=self.ics[pid],
                      current_state=current_state,
                      log=log)

    def set_reference(self, pid=0):
        self._reference_particle = self.__bundle_up(pid)

    def get_reference(self):
        return self._reference_particle

    def plot(self, Ylab='-D dK', Xlab='-D Theta', pids='all',
             mark_special=None, new_plot=True, **kwargs):
        """Mark special may be redundant in the light of
        Lattice::plot_segment.
        """

        ## reading how to plot data: diff variable with reference value, or otherwise
        import re
        x_flags = re.findall('-.', Xlab)
        y_flags = re.findall('-.', Ylab)
        Xlab = re.subn('-.* ', '', Xlab)[0]
        Ylab = re.subn('-.* ', '', Ylab)[0]

        ## reading the reference data from host ensemble
        try:
            p_ref = self.get_reference()
            if p_ref.log is None:
                raise ValueError
        except (AttributeError, ValueError):
            self.set_reference()
            p_ref = self.get_reference()
            print('Reference not set; using default (pid: {})'.format(p_ref.pid))

        ## selecting only the required pids for plot
        names = list(range(self.n_ics))
        names.insert(0, names.pop(p_ref.pid)) # move reference particle pid to front
        if pids != 'all':
            pids = set(pids)
            pids.add(p_ref.pid)
            names = set(names)
            not_found = names - pids
            names = names - not_found
            print("Discarded PIDs: " + ','.join([str(e) for e in not_found]))

        ## creating plot
        from matplotlib import pyplot as PLT

        if new_plot:
            PLT.figure()

        if mark_special is not None:
            elems = [re.sub('_.*', '', str(e)).replace("b'", '') for e in self[:, 0]['Element']]

            def color_map(elist):
                d = lambda e: 'red' if e == mark_special else 'black'
                return [d(e) for e in elist]

            plot = lambda X, Y, lab, **kwargs: PLT.scatter(X, Y, label=mark_special,
                                                           c=color_map(elems), **kwargs)
            legend = lambda lab: PLT.legend([mark_special])
        else:
            plot = lambda X, Y, lab, **kwargs: PLT.plot(X, Y, label=lab, **kwargs)
            legend = lambda lab: PLT.legend(lab)

        nc = len(self[:, 0][Ylab])
        dX = np.empty([nc])
        dY = np.empty([nc])

        not_nan = [not e for e in np.isnan(p_ref.log['Turn'])]

        for i in names:
            dY = copy.deepcopy(self[:, i][Ylab][not_nan])
            dX = copy.deepcopy(self[:, i][Xlab][not_nan])
            if '-D' in x_flags: dX -= p_ref.log[Xlab][not_nan]
            if '-D' in y_flags: dY -= p_ref.log[Ylab][not_nan]
            plot(dX, dY, i, **kwargs)
        legend(names)

        ## creating pretty labels
        sub_map = {'Theta': '\Theta', 'dK': '\Delta K'}

        def sub(*labels):
            labels = list(labels) #otherwise labels is a tuple, so can't do assignment
            for i in range(len(labels)):
                for k, v in sub_map.items(): labels[i] = labels[i].replace(k, v)
            return labels

        Xlab, Ylab = sub(Xlab, Ylab)
        if '-D' in x_flags:
            Xlab = '${0} - {0}_0$'.format(Xlab)
        else:
            Xlab = '${0}$'.format(Xlab)
        if '-D' in y_flags:
            Ylab = '${0} - {0}_0$'.format(Ylab)
        else:
            Ylab = '${0}$'.format(Ylab)

        PLT.xlabel(Xlab)
        PLT.ylabel(Ylab)
        PLT.grid()

        return PLT.gcf()

#%%

if __name__ == '__main__':
    state_list = StateList(Sz=1, dK=(0, 1e-4, 1), x=(-1e-3, 1e-3, 1))
    log = PLog.from_file('test')
