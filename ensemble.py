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

        self.log = Bundle()

        self.n_ics = len(state_list)
        self.n_var = len(state_list[0])

        self.ics = dict(zip(range(len(state_list)), state_list))

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
            log = lambda: None
            ics = list()
            for tbl in f.root.Logs:
                setattr(log, tbl.name, tbl.read())
                ic = list(tbl[0])[3:]
                ics.append(dict(zip(rhs.VAR_NAME, ic)))

            ens = cls(ics, particle)
            ens.log = log

        return ens

    def __bundle_up(self, pid):
        log = getattr(self.log, 'P'+str(pid), None)
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

        ## reading how to plot data: diff variable with reference value, or otherwise
        import re
        x_flags = re.findall('-.', Xlab)
        y_flags = re.findall('-.', Ylab)
        Xlab = re.subn('-.* ', '', Xlab)[0]
        Ylab = re.subn('-.* ', '', Ylab)[0]

        try:
            p_ref = self.get_reference()
        except AttributeError:
            self.set_reference()
            p_ref = self.get_reference()
            print('Reference not set; using default (pid: {})'.format(p_ref.pid))


        ## selecting only the required pids for plot
        names = self.list_names()
        names.insert(0, names.pop(p_ref.pid))
        if pids != 'all':
            pids = set(pids)
            pids.add(p_ref.pid)
            names = set(names)
            not_found = names - pids
            names = names - not_found
            print("Discarded PIDs: " + ','.join([str(e) for e in not_found]))

        ## creating plot
        from matplotlib import pyplot as PLT

        if new_plot: PLT.figure()

        if mark_special is not None:
            elems = [re.sub('_.*', '', str(e)).replace("b'", '') for e in self[0].log['Element']]

            def cm(elist):
                d = lambda e: 'red' if e == mark_special else 'black'
                return [d(e) for e in elist]

            plot = lambda X, Y, lab, **kwargs: PLT.scatter(X, Y, label=mark_special, c=cm(elems), **kwargs)
            legend = lambda lab: PLT.legend([mark_special])
        else:
            plot = lambda X, Y, lab, **kwargs: PLT.plot(X, Y, label=lab, **kwargs)
            legend = lambda lab: PLT.legend(lab)

        nc = len(self[0].log[Ylab])
        dX = np.empty([nc])
        dY = np.empty([nc])

        not_nan = [not e for e in np.isnan(p_ref.log['Turn'])]

        for i in names:
            dY = copy.deepcopy(self[i].log[Ylab][not_nan])
            dX = copy.deepcopy(self[i].log[Xlab][not_nan])
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
