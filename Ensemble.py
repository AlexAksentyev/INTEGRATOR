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
from scipy.integrate import odeint
import numpy as NP
import Particle as PCL
import RHS
import copy
from time import clock

from Utilities import Bundle
    
#%%
        
class StateList:
    def __init__(self, **kwargs):
        
        keys = kwargs.keys()
        
        # create defined variables
        ntot = 1
        argDict = dict()
        for key, val in kwargs.items():
            try: 
                lb, ub, num = val
                val = NP.linspace(lb,ub,num)
            except TypeError: # if key = value, set value for all pcls
                num = 1
            ntot *= num
            argDict.update({RHS.imap[key]: val})
            
        # make mesh
        mesh = dict(zip(keys, NP.meshgrid(*list(argDict.values()))))
            
        vartype = list(zip(RHS.varname, NP.repeat(float, RHS.varnum)))
        self.SL = NP.zeros(ntot+1, dtype=vartype) # +1 for genuine refrence particle
        
        #write data
        for key, value in mesh.items():
            self.SL[key] = NP.array([0]+value.reshape(ntot).tolist())
            
        self.SL[0]['Sz'] = 1
        
        #here i should remove any duplicate reference particles
#        self.SL = NP.unique(self.SL) # this works but messes up the order
            
        # convert to list of dicts for use with ensemble
        self.SL = [dict(zip(self.SL.dtype.names, x)) for x in self.SL]
            
    def __len__(self):
        return len(self.SL)
    
    def __getitem__(self, pid):
        return self.SL[pid]
    
    def __repr__(self):
        from pandas import DataFrame
        return str(DataFrame(self.SL))
    
    def as_list(self):
        states = list()
        for d in self.SL: states.append(list(d.values()))
        return states        
        
#%%

class Ensemble:
    
    def __init__(self, state_list, Particle=PCL.Particle()):
        self.Particle = Particle
        
        self.Log = Bundle()
        
        self.n_ics = len(state_list)
        self.n_var = len(state_list[0])
        
        self.ics = dict(zip(range(len(state_list)), state_list))
        
        self.TrackControl = Bundle(FWD=True, inner = True, breaks=101, ncut=0) 
                                # ncut should be an int representing the number of turns to keep
                                # 0, or None, means don't cut, a.k.a. keep full log in memory
        
    @classmethod
    def populate(cls, Particle, **kwargs):
        sl = StateList(**kwargs)
        return cls(sl, Particle)
    
    @classmethod
    def from_file(cls, filename, directory = './data/'):
        import tables as TBL
        
        filename = '{}{}.h5'.format(directory, filename)
        with TBL.open_file(filename) as f:
            pcl = f.root.Particle[0]
            pcl = PCL.Particle(pcl['Mass0'], pcl['KinEn0'], pcl['G'])
            Log = lambda: None
            ics = list()
            for tbl in f.root.Logs:
                setattr(Log, tbl.name, tbl.read())
                ic = list(tbl[0])[3:]
                ics.append(dict(zip(RHS.varname, ic)))
                
            ens = cls(ics, pcl)
            ens.Log = Log
            
        return ens
        
    def __bundle_up(self, pid):
        log = getattr(self.Log, 'P'+str(pid), None)
        try: current_state = log[len(log)-1]
        except TypeError: 
            current_state = self.ics[pid]
        return Bundle(PID = pid, 
                      ics = self.ics[pid], 
                      current_state = current_state,
                      Log = log)
        
    def count(self):
        return self.n_ics
        
    def setReference(self, name = 0): #by default, a particle with all s.v. = 0 except Sz is put in index 0
        self.__fRefPart = self.__bundle_up(name)
        
    def getReference(self):
        return self.__fRefPart
    
    def listNames(self):
        return list(self.ics.keys())
    
    def __getitem__(self, pid):
        return self.__bundle_up(pid)

    def __iter__(self):
        self.__current_pid = 0
        return self
    
    def __next__(self):
        last_pid = self.count()-1
        if self.__current_pid <= last_pid:
            result = self[self.__current_pid]
            self.__current_pid +=1
            return result
        else:
            raise StopIteration
        
    def __repr__(self):
        from pandas import DataFrame
        return str(DataFrame(self.ics).T)
        
    def plot(self, Ylab='-D dK', Xlab='-D Theta', pids='all', mark_special=None, new_plot = True, **kwargs):

        ## reading how to plot data: diff variable with reference value, or otherwise
        import re
        x_flags = re.findall('-.', Xlab)
        y_flags = re.findall('-.', Ylab)
        Xlab = re.subn('-.* ', '', Xlab)[0]
        Ylab = re.subn('-.* ', '', Ylab)[0]
        
        try:
            pr = self.getReference()
        except AttributeError:
            self.setReference()
            pr = self.getReference()
            print('Reference not set; using default (pid: {})'.format(pr.PID))
            
        
        ## selecting only the required pids for plot
        names = self.listNames()
        names.insert(0, names.pop(pr.PID))
        if pids != 'all':
            pids = set(pids)
            pids.add(pr.PID)
            names = set(names)
            not_found = names - pids
            names = names - not_found
            print("Discarded PIDs: " + ','.join([str(e) for e in not_found]))
        
        ## creating plot
        from matplotlib import pyplot as PLT
        
        if new_plot: PLT.figure()     
        
        if mark_special is not None:
            Elt = [re.sub('_.*','',str(e)).replace("b'",'') for e in self[0].Log['Element']]
        
            def cm(elist):
                d = lambda e: 'red' if e == mark_special else 'black'
                return [d(e) for e in elist]
            
            plot = lambda X, Y, lab, **kwargs: PLT.scatter(X,Y, label=mark_special, c=cm(Elt), **kwargs)
            legend = lambda lab: PLT.legend([mark_special])
        else:
            plot = lambda X, Y, lab, **kwargs: PLT.plot(X,Y,label=lab, **kwargs)
            legend = lambda lab: PLT.legend(lab)
                
        nc = len(self[0].Log[Ylab])
        dX = NP.empty([nc])
        dY = NP.empty([nc])
        
        not_nan = [not e for e in NP.isnan(pr.Log['Turn'])]
        
        for i in names:
            dY = copy.deepcopy(self[i].Log[Ylab][not_nan])
            dX = copy.deepcopy(self[i].Log[Xlab][not_nan])
            if '-D' in x_flags: dX -= pr.Log[Xlab][not_nan]
            if '-D' in y_flags: dY -= pr.Log[Ylab][not_nan]
            plot(dX, dY, i, **kwargs)
            None
        legend(names)
        
        ## creating pretty labels
        sub_map = {'Theta':'\Theta','dK':'\Delta K'}
        
        def sub(*labels):
            labels = list(labels)
            for i in range(len(labels)):
                for k,v in sub_map.items(): labels[i] = labels[i].replace(k,v)
            return labels
        
        Xlab,Ylab = sub(Xlab,Ylab)
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
if __name__ is '__main__':

    pass
