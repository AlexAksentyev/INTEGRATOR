#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 15:04:54 2017

@author: alexa

TODO:
    * remove duplicate reference particle
    * vectorize tilted/untilted lattice:
        this'll require that all elements have to be tilted the same number of times,
        so as to have an equal number of tilt matrices,
        and hence length of state vector in the RHS,
        in all elements

"""
import copy
import numpy as np
import particle as pcl
import rhs

from utilities import Bundle
from particle_log import PLog

#%%

class StateList:
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

class Ensemble:

    def __init__(self, state_list, particle=pcl.Particle()):
        self.particle = particle

        self.log = None # Bundle()

        self.n_ics = len(state_list)
        self.n_var = len(state_list[0])

        self.ics = dict(zip(range(len(state_list)), state_list))
        
        self.log = None # a PLog created by tracker in track

    @classmethod
    def populate(cls, particle, **kwargs):
        state_list = StateList(**kwargs)
        return cls(state_list, particle)

    @classmethod
    def from_file(cls, filename, directory='./data/'):
        import tables as TBL

        filename = '{}{}.h5'.format(directory, filename)
        with TBL.open_file(filename) as f:
            particle = f.root.Particle[0]
            particle = pcl.Particle(particle['Mass0'], particle['KinEn0'], particle['G'])
            log, ic_dict = PLog.from_file(f)

            ens = cls(ic_dict, particle)
            ens.log = log

        return ens

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, copy.deepcopy(v, memo))

        if isinstance(self.log, PLog):   ## PLog keeps a reference to host ensemble,
            self.log.update_host(result) ## hence when deepcopying we want that reference
        return result                    ## to update to the current ensemble copy

    def __bundle_up(self, pid):
        log = self.log[:, pid]
        try:
            current_state = log[len(log)-1]
        except TypeError:
            current_state = self.ics[pid]
        return Bundle(pid=pid,
                      ics=self.ics[pid],
                      current_state=current_state,
                      log=log)

    def count(self):
        return self.n_ics

    def set_reference(self, name=0): #by default, a particle with all s.v. = 0 except Sz is put in index 0
        self.__reference_particle = self.__bundle_up(name)

    def get_reference(self):
        return self.__reference_particle

    def list_names(self):
        return list(self.ics.keys())

    def __getitem__(self, pid):
        return self.__bundle_up(pid)

    def __iter__(self):
        self.__current_pid = 0
        return self

    def __next__(self):
        last_pid = self.count() - 1
        if self.__current_pid <= last_pid:
            result = self[self.__current_pid]
            self.__current_pid += 1
            return result
        else:
            raise StopIteration

    def __repr__(self):
        from pandas import DataFrame
        return str(DataFrame(self.ics).T)
    
    def plot(self, Ylab='-D dK', Xlab='-D Theta', pids='all', mark_special=None, new_plot=True, **kwargs):
        """Wrapper for PLog::plot.
        """
        if self.log is None: 
            print('Nothing to plot.')
            return
        self.log.plot(Ylab='-D dK', Xlab='-D Theta', pids='all', mark_special=None, new_plot=True, **kwargs)
