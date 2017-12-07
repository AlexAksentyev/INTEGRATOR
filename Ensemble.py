#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 15:04:54 2017

@author: alexa

TODO:        
    * vectorize tilted/untilted lattice

"""
from scipy.integrate import odeint
import numpy as NP
import Particle as PCL
import RHS
import copy
from time import clock

class Bundle(dict):
    """ Bundle serves as an interface for easy access 
        to a bundle of ensemble particle data
    """
    def __init__(self,**kw):
        dict.__init__(self,kw)
        self.__dict__ = self
        
    def __deepcopy__(self, memo):
        res = self.__class__(**self.__dict__)
        return res
        
    def __repr__(self):
        return str(list(self.keys()))
    
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
    
    def __setup_file(self, file_handle):
        # write particle parameters (Mass0, KinEn0, G)
        ppar = self.Particle.getParams() 
        file_handle.create_table('/','Particle',ppar)
        # separate group to write p-logs in
        file_handle.create_group('/','Logs', 'Particle logs')
    
    def saveData(self, filename, directory='./data/'):
        import tables
        
        filename = '{}{}.h5'.format(directory, filename)
        with tables.open_file(filename, mode='w') as f:
            self.__setup_file(f)
            for p in self: f.create_table(f.root.Logs, 'P'+str(p.PID), p.Log)
        
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
    
    def setTrackArgs(self, **kwargs):
        self.TrackControl.update(**kwargs)

    def track(self, ElementSeq , ntimes):
        """ if cut (see TrackControl) is true, particle logs keep only the values 
            relevant for the current turn, else keep all values in RAM
        """
        # get lattice name for saving data into file
        from Element import Lattice, ERF
        if isinstance(ElementSeq, Lattice):
            latname = ElementSeq.Name
            cnt = ElementSeq.RFCount
            RF = ElementSeq.getRF()
            ElementSeq = ElementSeq.Sequence
        else:
            latname = 'Unnamed_sequence'
            RF = None; cnt = 0
            for e in ElementSeq:
                if not isinstance(e, ERF):
                    continue
                RF = e; cnt += 1
        if cnt > 1: 
            print('\t\t More than one ({}) RF;\n aborting.'.format(cnt))
            return
        
        # number integration period subdivisions
        brks = self.TrackControl.breaks
        self.IntBrks = brks # Ensemble.IntBrks will be used in RHS
        
        ## creation of particle logs
        # recarray type
        import tables as TBL
        names = ['START']+[e.fName for e in ElementSeq]
        n = len(names[NP.argmax(names)])
        EType = TBL.StringCol(n) #otherwise problems writing into hdf5 file
        vartype = [('Turn',int),('Element',EType),('Point', int)]
        vartype += list(zip(RHS.varname, NP.repeat(float, RHS.varnum)))
        
        self.__fLastPnt = -1 # marker of the state vector upon exiting an element
        
        # counting the number of records in a p-log
        n_elem = len(ElementSeq)
        inner = self.TrackControl.inner
        ncut = self.TrackControl.ncut
        if inner: 
            nrow = n_elem*self.IntBrks
            if ncut == 0:  nrow *= ntimes 
            else: nrow *= ncut
            ind = 0 # odeint will return injection values
        else: 
            nrow = n_elem
            if ncut == 0: nrow *= ntimes
            else: nrow *= ncut
            nrow += 1 # +1 for injection values
            ind = 1 # odeint won't return injection values; set them manually
            
        cut = True # reset log index after writing to file
        if ncut == 0: 
            ncut = NP.floor(ntimes/10) # ncut is used below to decide whether to write data
                        # if we keep all data in RAM, still backup every 10% of turns
                        # 10% is arbitrary
            cut = False 
                                    
        print('Saving data to file {} every {} turns'.format(latname, ncut))
        
        # creating particle logs
        ics = list()
        for pid, ic in self.ics.items():
            setattr(self.Log, 'P'+str(pid), NP.recarray(nrow,dtype=vartype))
            self[pid].Log.fill(NP.nan) # in case we later get a NaN error in tracking, 
                                       # pre-fill the log with nans
                                       # Turn,Point fill fill with a big random integer
            ic = list(ic.values())
            self[pid].Log[0] = 0,names[0],self.__fLastPnt, *ic # saving injection values (will be overwritten if inner is true)
            ics.append(ic)
        
        # current state vector
        state = NP.array(ics) # [[x0,y0,...], [x1,y1,...], [x2,y2,...]]
        n_ics = self.n_ics
        n_var = self.n_var
        
        rhs = RHS.RHS(self, RF) # setting up the RHS
        
        # opening hdf5 file to output data
        # write the used particle parameters (Mass0, KinEn0, G)
        # write particle logs
        filename = './data/{}.h5'.format(latname)
        with TBL.open_file(filename, 'w', latname) as f: 
            self.__setup_file(f)
            for p in self: # creating data tables to fill
                f.create_table(f.root.Logs, 'P'+str(p.PID), p.Log.dtype)
            
            old_ind=0 # log is written into file from old_ind:ind
            old_turn = 0 # every ncut turns
            old_percent = -1 # to print 0 %
            cnt = 0; cnt_tot = n_elem*ntimes # for progress bar
            print('\t\t LATTICE: {} '.format(latname))
            FWD = self.TrackControl.FWD
            for n in range(1,ntimes+1): # turn
                for i in range(len(ElementSeq)): # element
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
                        state = state.reshape(n_ics*n_var) # flat [x0,y0,...,x1,y1,...,x2,y2]
                        if not bERF:
                            vals = odeint(rhs, state, at, args=(element,))
                            state = vals[brks-1] # [x0,y0,...,x1,y1,...]
                        else:
                            element.advance(state)
                        state = state.reshape(n_ics,n_var) # [[x0,y0,...], [x1,y1,...], [x2,y2,...]]
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
                # end element loop
                    
                if (n-old_turn)%ncut == 0:
                    print('turn {}, writing data ...'.format(n))
                    start = clock()
                    # write data to file
                    for p in self:
                        tbl = getattr(f.root.Logs, 'P'+str(p.PID))
                        tbl.append(p.Log[old_ind:ind])
                        tbl.flush()
                    if cut: ind = 0 # if cut, logs can only keep data for one turn => overwrite
                                    # otherwise, keep writing log
                    old_ind = ind
                    old_turn = n
                    print('old index {}'.format(old_ind))
                    print('writing took: {:04.2f} secs'.format(clock()-start))
                
            # end turn loop
            
            ## write any remainig data
            print('writing remaining ({}) data ...'.format(ind-old_ind))
            start = clock()
            # write data to file
            for p in self:
                tbl = getattr(f.root.Logs, 'P'+str(p.PID))
                tbl.append(p.Log[old_ind:ind])
                tbl.flush()
            print('writing took: {:04.2f} secs'.format(clock()-start))
                    
        print('Complete 100 %')
        
#%%
if __name__ is '__main__':
    import Element as ENT
    import Particle as PCL
    from matplotlib import pyplot as PLT
    
#    s = StateList(dK=(0e-3,3e-4,5), x=(-1e-3,1e-3,2), Sz=1)
    
    E = Ensemble.populate(PCL.Particle(), dK=(0e-3,3e-4,5), x=(-1e-3,1e-3,2), Sz=1)
    R3 = ENT.Wien(361.55403e-2,5e-2,PCL.Particle(),-120e5,.082439761)
    OD1 = ENT.Drift(.25, 'OD1')
    QD1 = ENT.MQuad(5e-2,-.82,"QD")
    QF1 = ENT.MQuad(5e-2,.736,"QF")
    
    FODO = ENT.Lattice([QF1, OD1, QD1, OD1], 'FODO')
    FODO.insertRF(0,0,E,EField=15e7)
#    
##%%
    E.track(FODO,int(3e2))
#    R3.tilt('x',5)
#    tE = copy.deepcopy(E)
#    tE.track([R3,OD1,OD1,OD1,OD1],int(5e1),cut=False)
#    
##%%
    E.setReference()
    E.plot('y','s')
#    tE.setReference(0)
#    pids = [1,4,5]
#    n=1
#    for pid in pids:
#        PLT.subplot(3,1,n); n += 1
#        tE.plot('-D x','s', pids=[E.getReference().PID, pid], new_plot=False)
