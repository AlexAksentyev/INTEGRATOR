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

class Ensemble:
    
    def __init__(self, state_list, Particle=PCL.Particle()):
        self.Particle = Particle
        self.__fRefPart = None
        
        self.Log = lambda: None 
        
        self.n_ics = len(state_list)
        self.n_var = len(state_list[0])
        
        self.ics = dict(zip(range(len(state_list)), state_list))
        
    def count(self):
        return self.n_ics
    
    def set(self, pid, **kwargs):
        self.ics[pid].update(**kwargs)
        
    def setReference(self, name):
        self.__fRefPart = getattr(self.Log, 'P'+str(name))
        self.__fPID = name
        
    def getReference(self):
        return self.__fRefPart
    
    def listNames(self):
        return list(self.ics.keys())
    
    def __getitem__(self, pid):
        return getattr(self.Log, 'P'+str(pid))
    
    def __setitem__(self, pid, name_value):
        from numpy.lib.recfunctions import append_fields
        
        name = name_value[0]
        value = name_value[1]
        
        log = self[pid]
        setattr(self.Log, 'P'+str(pid), append_fields(log, name, value,
                            dtypes=float, usemask=False, asrecarray=True))
        
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
        
    def plot(self, Ylab='-D dK', Xlab='-D Theta', pids='all', mark_special=None, **kwargs):

        ## reading how to plot data: diff variable with reference value, or otherwise
        import re
        x_flags = re.findall('-.', Xlab)
        y_flags = re.findall('-.', Ylab)
        Xlab = re.subn('-.* ', '', Xlab)[0]
        Ylab = re.subn('-.* ', '', Ylab)[0]
        
        try:
            pid0 = self.__fPID
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
        
        
        Elt = [re.sub('_.*','',e) for e in self[0]['Element']]
        
        def cm(elist):
            d = lambda e: 'red' if e == mark_special else 'black'
            return [d(e) for e in elist]
        
        
        nr = self.count()
        nc = len(self[0][Ylab])
                
        Y = NP.empty([nr,nc])
        X = NP.empty([nr,nc])
        dX = NP.empty([nc])
        dY = NP.empty([nc])
        
        from matplotlib import pyplot as PLT
        
        PLT.figure()     
        
        if mark_special is not None:
            plot = lambda X, Y, lab, **kwargs: PLT.scatter(X,Y, label=mark_special, c=cm(Elt), **kwargs)
            legend = lambda lab: PLT.legend([mark_special])
        else:
            plot = lambda X, Y, lab, **kwargs: PLT.plot(X,Y,label=lab, **kwargs)
            legend = lambda lab: PLT.legend(lab)
        
        for i in names:
            Y[i] = self[i][Ylab]
            dY = copy.deepcopy(Y[i])
            X[i] = self[i][Xlab]
            dX = copy.deepcopy(X[i])
            if '-D' in x_flags: dX -= X[pid0]
            if '-D' in y_flags: dY -= Y[pid0]
            plot(dX, dY, i, **kwargs)
            
        legend(names)
        
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
            
        return (X, Y, PLT.gcf())

    def track(self, ElementSeq , ntimes, FWD=True, inner = True, breaks=101):
        from Element import Lattice
        if type(ElementSeq == Lattice):
            ElementSeq = ElementSeq.fSequence
        
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
                self[pid][0] = 0,names[0],self.__fLastPnt, *ic
                ics.append(ic)
            ind = 1
        
        state = NP.array(ics)
        
        n_ics = self.n_ics
        n_var = self.n_var
        
        rhs = RHS.RHS(self)
        
        old_percent = -1
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
                                self[pid][ind] = n,element.fName, k, *valsk[pid]
                            ind += 1
                    for pid in self.ics.keys():
                        self[pid][ind] = n,element.fName, self.__fLastPnt, *state[pid]
                    ind += 1
                except ValueError:
                    print('NAN error at: Element {}, turn {}'.format(element.fName, n))
                    for m in range(ind,len(self.Log.P0)):
                        for pid in self.ics.keys():
                            self[pid][ind] = n, element.fName, self.__fLastPnt, *([NP.NaN]*(len(vartype)-3))
                        ind += 1
                    return
                
        print('Complete 100 %')
        
        pr = self.Particle
        
        try:
            check = pr.fRF
        except AttributeError:
            print('\n \t \t System w/o RF')
            check = {'Freq':0, 'Phase':0}
        
        th = lambda t: 2*NP.pi*check['Freq']*t + check['Phase']
        for pid in self.ics.keys():
            self[pid] = ('Theta', th(self[pid].t))
#%%
if __name__ is '__main__':
    import utilFunc as U
    import Element as ENT
    
    states=[list(e.values()) for e in U.form_state_list(xint=(1e-3,1e-3),yint=(-1e-3,-1e-3))]
    
    E = Ensemble(states)
    R3 = ENT.Wien(361.55403e-2,5e-2,PCL.Particle(),-120e5,.082439761)
    OD1 = ENT.Drift(.25, 'OD1')
    QD1 = ENT.MQuad(5e-2,-.82,"QD")
    QF1 = ENT.MQuad(5e-2,.736,"QF")
    
    E.track([QF1, OD1, OD1, OD1],100)
    
    E.plot_min('x','s')
