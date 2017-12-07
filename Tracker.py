#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:21:59 2017

@author: alexa

TODO:
    * Tracker should have Particle Logs, a pointer copy to which is given to Ensemble 
    at the end of tracking, and deleted

"""

from time import clock
import RHS
import numpy as NP
from scipy.integrate import odeint

from Utilities import Bundle, Writer

class Tracker:
    
    def __init__(self, Lattice, Ensemble):
        self.Lattice = Lattice
        self.Ensemble = Ensemble     
        
        
    
    def setControls(self, inner=True, FWD = True, breaks=101, ncut=0):
        self.Controls = Bundle(inner=inner,FWD=FWD,breaks=breaks,ncut=ncut)  

    def __form_log_type(self, injection_name):
        import tables as TBL
        names = [injection_name]+[e.fName for e in self.Lattice]
        n = len(names[NP.argmax(names)])
        EType = TBL.StringCol(n) #otherwise problems writing into hdf5 file
        vartype = [('Turn',int),('Element',EType),('Point', int)]
        vartype += list(zip(RHS.varname, NP.repeat(float, RHS.varnum)))
        
        return vartype
        
    def __getitem__(self, pid):
        return getattr(self.Log, 'P'+str(pid), None)
    
    def __create_logs(self, ntimes):
        """ returns first index from which to fill log
        """
        ElementSeq = self.Lattice.Sequence
        ## creating p-logs
        # counting the number of records in a p-log
        brks = self.Controls.breaks
        n_elem = len(ElementSeq)
        inner = self.Controls.inner
        ncut = self.Controls.ncut
        if inner: 
            nrow = n_elem*brks
            if ncut == 0:  nrow *= ntimes 
            else: nrow *= ncut
            ind = 0 # odeint will return injection values
        else: 
            nrow = n_elem
            if ncut == 0: nrow *= ntimes
            else: nrow *= ncut
            nrow += 1 # +1 for injection values
            ind = 1 # odeint won't return injection values; set them manually
            
        inj_name = 'START'
        self.__fLastPnt = -1 # marker of the state vector upon exiting an element
        vartype = self.__form_log_type(inj_name)
        Ensemble = self.Ensemble
        # creating particle logs
        ics = list()
        for pid, ic in Ensemble.ics.items():
            setattr(self.Log, 'P'+str(pid), NP.recarray(nrow,dtype=vartype))
            self[pid].fill(NP.nan) # in case we later get a NaN error in tracking, 
                                       # pre-fill the log with nans
                                       # Turn,Point fill fill with a big random integer
            ic = list(ic.values())
            self[pid][0] = 0,inj_name,self.__fLastPnt, *ic # saving injection values (will be overwritten if inner is true)
            ics.append(ic)
            
        self.ics = ics
            
        return ind
    
    def write_log(self, from_to):
        old_ind, ind = from_to
        # write data to file
        for p in self.Ensemble:
            tbl = getattr(self.file_handle.root.Logs, 'P'+str(p.PID))
            tbl.append(p.Log[old_ind:ind])
            tbl.flush()
    
    def track(self, ntimes):
        Lattice = self.Lattice
        latname = Lattice.Name
        
        ## check if there's more than one RF element
        if Lattice.RFCount > 1:
            print('\t\t More than one ({}) RF;\n aborting.'.format(Lattice.RFCount))
            return
        
        ind = self.__create_logs(ntimes)
        
        ncut = self.Controls.ncut    
        cut = True # reset log index after writing to file
        if ncut == 0: 
            ncut = NP.floor(ntimes/10) # ncut is used below to decide whether to write data
                        # if we keep all data in RAM, still backup every 10% of turns
                        # 10% is arbitrary
            cut = False 
                                   
        print('Saving data to file {} every {} turns'.format(latname, ncut))
        
         # current state vector
        state = NP.array(self.ics) # [[x0,y0,...], [x1,y1,...], [x2,y2,...]]
        n_ics = len(state)
        n_var = len(state[0])
        
        # create the RHS
        RF = Lattice.getRF()
        Ensemble = self.Ensemble
        rhs = RHS.RHS(Ensemble, RF) # setting up the RHS
        
        # opening hdf5 file to output data
        # write the used particle parameters (Mass0, KinEn0, G)
        # write particle logs
        filename = './data/{}.h5'.format(latname)
        ElementSeq = Lattice.Sequence
        n_elem = len(ElementSeq)
        import tables as TBL
        with TBL.open_file(filename, 'w', latname) as f: 
            writer = Writer(Ensemble, Lattice, f)
            old_ind=0 # log is written into file from old_ind:ind
            old_turn = 0 # every ncut turns
            old_percent = -1 # to print 0 %
            cnt = 0; cnt_tot = n_elem*ntimes # for progress bar
            print('\t\t LATTICE: {} '.format(latname))
            FWD = self.Controls.FWD
            brks = self.Controls.breaks
            inner = self.Controls.inner
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
                                    self[pid][ind] = n,element.fName, k, *valsk[pid]
                                ind += 1
                        for pid in self.ics.keys():
                            self[pid][ind] = n,element.fName, self.__fLastPnt, *state[pid]
                        ind += 1
                    except ValueError:
                        print('NAN error at: Element {}, turn {}, log index {}'.format(element.fName, n, ind))
                        return
                # end element loop
                    
                if (n-old_turn)%ncut == 0:
                    print('turn {}, writing data ...'.format(n))
                    start = clock()
                    writer.write_log((old_ind,ind))
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
            writer.write_log((old_ind,ind))
            print('writing took: {:04.2f} secs'.format(clock()-start))
                    
        print('Complete 100 %')
    
