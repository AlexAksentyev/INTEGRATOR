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
        
    @classmethod
    def from_file(cls, filename):
        import tables as TBL
        
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
            pr = self.getReference()
        except AttributeError:
            print('Reference not set')
            return
        
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

    def track(self, ElementSeq , ntimes, FWD=True, inner = True, breaks=101, cut=False):
        """ if cut is true, particle logs keep only the values relevant for the current turn,
            else keep all values in RAM
        """
        from Element import Lattice
        if type(ElementSeq) == Lattice:
            latname = ElementSeq.Name
            ElementSeq = ElementSeq.fSequence
        else:
            latname = 'Unnamed_sequence'
        
        brks = breaks
        self.fIntBrks = brks
        
        import tables as TBL
        
        names = ['START']+[e.fName for e in ElementSeq]
        n = len(names[NP.argmax(names)])
        EType = TBL.StringCol(n)
        vartype = [('Turn',int),('Element',EType),('Point', int)]
        vartype += list(zip(RHS.varname, NP.repeat(float, RHS.varnum)))
        
        self.__fLastPnt = -1
        
        ics = list()
        n_elem = len(ElementSeq)
        if inner: 
            nrow = n_elem*self.fIntBrks
            if not cut:  nrow *= ntimes
            ind = 0 # odeint will return injection values
        else: 
            nrow = n_elem
            if not cut: nrow *= ntimes
            nrow += 1 # +1 for injection values
            ind = 1 # odeint won't return injection values; set them manually
        
        # check for memory error
        # if so, split log into chunks
        try: NP.recarray(nrow*self.n_ics,dtype=vartype)
        except MemoryError:
            cut = True
            print('Too much memory required; cutting logs')
            if inner: nrow /= ntimes
            else: nrow -=1; nrow /= ntimes; nrow += 1
        
        nrow = int(nrow)
        
        for pid, ic in self.ics.items():
                setattr(self.Log, 'P'+str(pid), NP.recarray(nrow,dtype=vartype))
                self[pid].Log.fill(NP.nan) # in case we later get a NaN error in tracking, 
                                        # pre-fill the log with nans
                ic = list(ic.values())
                self[pid].Log[0] = 0,names[0],self.__fLastPnt, *ic # saving injection values
                ics.append(ic)
        
        state = NP.array(ics)
        
        n_ics = self.n_ics
        n_var = self.n_var
        
        rhs = RHS.RHS(self)
        
        filename = './data/{}.h5'.format(latname)
        with TBL.open_file(filename, 'w', latname) as f: 
            ppar = self.Particle.getParams()
            f.create_table('/','Particle',ppar)
            f.create_group('/','Logs', 'Particle logs')
            for p in self: # creating data tables
                f.create_table(f.root.Logs, 'P'+str(p.PID), p.Log.dtype)
            
            old_percent = -1
            cnt = 0
            cnt_tot = n_elem*ntimes
            old_ind=0
            for n in range(1,ntimes+1):
                for i in range(len(ElementSeq)):       
                    # pick element
                    if FWD: element = ElementSeq[i]
                    else: element = ElementSeq[len(ElementSeq)-1-i]
                    
                    # print computation progress
                    percent = int(cnt/cnt_tot*100)
                    if percent%10 == 0 and percent != old_percent:
                        print('Complete {} %'.format(percent))
                        old_percent = percent
                    cnt += 1
                    
                    # choose if integrate or otherwise advance the state vector
                    bERF = element.bSkip
                    
                    #integrate at these points
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
                                valsk = vals[k].reshape(n_ics, n_var)
                                for pid in self.ics.keys():
                                    self[pid].Log[ind] = n,element.fName, k, *valsk[pid]
                                ind += 1
                        for pid in self.ics.keys():
                            self[pid].Log[ind] = n,element.fName, self.__fLastPnt, *state[pid]
                        ind += 1
                    except ValueError:
                        print('NAN error at: Element {}, turn {}, log index {}'.format(element.fName, n, ind))
                        return
                    
                # write data to file
                for p in self:
                    tbl = getattr(f.root.Logs, 'P'+str(p.PID))
                    tbl.append(p.Log[old_ind:ind])
                    tbl.flush()
                if cut: ind = 0 # if cut, logs can only keep data for one turn => overwrite
                                # otherwise, keep writing log
                old_ind = ind
                    
        print('Complete 100 %')
        
#%%
if __name__ is '__main__':
    import utilFunc as U
    import Element as ENT
    from matplotlib import pyplot as PLT
    

    states = U.StateList(dK=(0e-3,3e-4,5), x=(-1e-3,1e-3,2), Sz=1)
    
    E = Ensemble(states)
    R3 = ENT.Wien(361.55403e-2,5e-2,PCL.Particle(),-120e5,.082439761)
    OD1 = ENT.Drift(.25, 'OD1')
    QD1 = ENT.MQuad(5e-2,-.82,"QD")
    QF1 = ENT.MQuad(5e-2,.736,"QF")
    
    FODO = ENT.Lattice([QF1, OD1, QD1, OD1], 'FODO')
    FODO.insertRF(0,0,E,EField=15e4)
    
#%%
    E.track(FODO,int(1e3),cut=False)
    
#%%
    E.setReference(0)
    pids = [1,4,5]
    n=1
    for pid in pids:
        PLT.subplot(3,1,n); n += 1
        E.plot('-D y','s', pids=[E.getReference().PID, pid], new_plot=False)
