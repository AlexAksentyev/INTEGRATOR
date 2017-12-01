#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 15:04:54 2017

@author: alexa
"""
from scipy.integrate import odeint
import numpy as NP
import Particle as PCL
import RHS
import copy

class Bundle(dict):
    """Bunch serves as an interface for easy access 
    to a bundle of data of a particle from ensemble
    """
    def __init__(self,**kw):
        dict.__init__(self,kw)
        self.__dict__ = self
        

class Ensemble:
    
    def __init__(self, state_list, Particle=PCL.Particle()):
        self.Particle = Particle
        self.__fRefPart = None
        
        self.Log = lambda: None 
        
        self.n_ics = len(state_list)
        self.n_var = len(state_list[0])
        
        self.ics = dict(zip(range(len(state_list)), state_list))
        
    def __bundle_up(self, pid):
        log = getattr(self.Log, 'P'+str(pid), None)
        lst_i = len(log)-1
        return Bundle(PID = pid, 
                      ics = self.ics[pid], 
                      current_state = log[lst_i],
                      Log = log)
        
    def count(self):
        return self.n_ics
        
    def setReference(self, name):
        self.__fRefPart = self.__bundle_up(name)
        
    def getReference(self):
        return self.__fRefPart
    
    def listNames(self):
        return list(self.ics.keys())
    
    def __getitem__(self, pid):
        return self.__bundle_up(pid)
    
    def __iter__(self):
        self.current_pid = 0
        return self
    
    def __next__(self):
        last_pid = self.count()-1
        if self.current_pid <= last_pid:
            result = self[self.current_pid]
            self.current_pid +=1
            return result
        else:
            raise StopIteration
        
    def __repr__(self):
        from pandas import DataFrame
        
        return str(DataFrame(self.ics).T)
    
    def saveData(self, filename):
        import tables
        
        filename = './data/{}.h5'.format(filename)
        with tables.open_file(filename, mode='w') as f:
            for p in self: f.create_table(f.root, 'P'+str(p.PID), p.Log)
        
    def plot(self, Ylab='-D dK', Xlab='-D Theta', pids='all', mark_special=None, new_plot = True, **kwargs):

        ## reading how to plot data: diff variable with reference value, or otherwise
        import re
        x_flags = re.findall('-.', Xlab)
        y_flags = re.findall('-.', Ylab)
        Xlab = re.subn('-.* ', '', Xlab)[0]
        Ylab = re.subn('-.* ', '', Ylab)[0]
        
        try:
            pid0 = self.getReference().PID
        except AttributeError:
            print('Reference not set')
            return
        
        ## selecting only the required pids for plot
        names = self.listNames()
        names.insert(0, names.pop(pid0))
        if pids != 'all':
            pids = set(pids)
            pids.add(pid0)
            names = set(names)
            not_found = names - pids
            names = names - not_found
            print("Discarded PIDs: " + ','.join([str(e) for e in not_found]))
        
        ## creating plot
        from matplotlib import pyplot as PLT
        
        if new_plot: PLT.figure()     
        
        if mark_special is not None:
            Elt = [re.sub('_.*','',e) for e in self[0].Log['Element']]
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
        
        for i in names:
            dY = copy.deepcopy(self[i].Log[Ylab])
            dX = copy.deepcopy(self[i].Log[Xlab])
            if '-D' in x_flags: dX -= self[pid0].Log[Xlab]
            if '-D' in y_flags: dY -= self[pid0].Log[Ylab]
            plot(dX, dY, i, **kwargs)
            
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

    def track(self, ElementSeq , ntimes, FWD=True, inner = True, breaks=101):
        from Element import Lattice
        if type(ElementSeq) == Lattice:
            filename = ElementSeq.Name # for saving data
            ElementSeq = ElementSeq.fSequence
        else: filename = 'Unnanmed_sequence'
        
        brks = breaks
        
        self.fIntBrks = brks
        
        names = ['START']+[e.fName for e in ElementSeq]
        n = str(len(names[NP.argmax(names)]))
        EType = 'U'+n
        vartype = [('Turn',int),('Element',EType),('Point', int)]
        vartype += list(zip(RHS.varname, NP.repeat(float, RHS.varnum)))
        
        self.__fLastPnt = -1
        
        ics = list()
        if inner: 
            nrow = ntimes*len(ElementSeq)*self.fIntBrks
            for pid, ic in self.ics.items():
                setattr(self.Log, 'P'+str(pid), NP.recarray(nrow,dtype=vartype))
                ics.append(list(ic.values()))
            ind = 0
        else: 
            nrow = ntimes*len(ElementSeq) + 1
            for pid, ic in self.ics.items():
                setattr(self.Log, 'P'+str(pid), NP.recarray(nrow,dtype=vartype))
                ic = list(ic.values())
                self[pid].Log[0] = 0,names[0],self.__fLastPnt, *ic
                ics.append(ic)
            ind = 1
        
        state = NP.array(ics)
        
        n_ics = self.n_ics
        n_var = self.n_var
        
        rhs = RHS.RHS(self)

        old_percent = -1
        ## tracking
        for n in range(1,ntimes+1):
            for i in range(len(ElementSeq)):
                percent = int(ind/nrow*100)
                if percent%10 == 0 and percent != old_percent:
                    print('Complete {} %'.format(percent))
                    old_percent = percent
                    
                if FWD: element = ElementSeq[i]
                else: element = ElementSeq[len(ElementSeq)-1-i]
                
                bERF = element.bSkip
                
                at = NP.linspace(0, element.fLength, brks)
                
                try:
                    element.frontKick(state)
                    state = state.reshape(n_ics*n_var)
                    if not bERF:
                        vals = odeint(rhs, state, at, args=(element,))
                        state = vals[brks-1]
                    else:
                        element.advance(state)
                    state = state.reshape(n_ics,n_var)
                    element.rearKick(state)
                    if not bERF and inner:
                        for k in range(brks-1):
                            valsk = vals[k].reshape(n_ics, n_var, order='C')
                            for pid in self.ics.keys():
                                self[pid].Log[ind] = n,element.fName, k, *valsk[pid]
                            ind += 1
                    for pid in self.ics.keys():
                        self[pid].Log[ind] = n,element.fName, self.__fLastPnt, *state[pid]
                    ind += 1
                except ValueError:
                    print('NAN error at: Element {}, turn {}'.format(element.fName, n))
                    for m in range(ind,len(self.Log.P0)):
                        for pid in self.ics.keys():
                            self[pid].Log[ind] = n, element.fName, self.__fLastPnt, *([NP.NaN]*(len(vartype)-3))
                        ind += 1
                    return
                
        print('Complete 100 %')
        
#%%
if __name__ is '__main__':
    import utilFunc as U
    import Element as ENT
    from matplotlib import pyplot as PLT
    
    states = U.StateList(dK=(0e-3,3e-4,2), x=(-1e-3,-1e-3,2), y=(-1e-3,1e-3,2), Sz=1)
    
    E = Ensemble(states)
    R3 = ENT.Wien(361.55403e-2,5e-2,PCL.Particle(),-120e5,.082439761)
    OD1 = ENT.Drift(.25, 'OD1')
    QD1 = ENT.MQuad(5e-2,-.82,"QD")
    QF1 = ENT.MQuad(5e-2,.736,"QF")
    
    E.track([QF1, OD1, OD1, OD1],100)
    
    #%%
    E.setReference(0)
    pids = [1,4,5]
    n=1
    for pid in pids:
        PLT.subplot(3,1,n); n += 1
        E.plot('-D y','s', pids=[E.getReference().PID, pid], new_plot=False)
        
