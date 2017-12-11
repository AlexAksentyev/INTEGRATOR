#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 09:01:42 2017

@author: alexa
"""
import numpy as np
from tables import StringCol
import rhs


class PLog(np.recarray):
    """Particle log is a 3D recarray, where the 1st index refers to the particle (ics),
    the 2nd to the row (record), the 3rd to the variable (like x,y,s, &c).
    """
    
    def __new__(cls, ics_dict, n_records):
        
        max_len_name = 10
        el_field_type = StringCol(max_len_name) #otherwise problems writing into hdf5 file
        metadata_type = [('Turn', int), ('Element', el_field_type), 
                   ('EID', int) , ('Point', int)]
        variable_type = list(zip(rhs.VAR_NAME, np.repeat(float, rhs.VAR_NUM)))
        record_type = metadata_type + variable_type
        
        n_ics = len(ics_dict)
        
        obj = np.recarray((n_records, n_ics), dtype=record_type, order='C').view(cls)
        super(PLog, obj).fill(np.nan)
        obj._last_pnt_marker = -1
        obj.n_ics = n_ics
        
        ics = list()
        for ind, ic in enumerate(ics_dict.values()):
            ics.append(list(ic.values()))
        
        ics = np.array(ics).flatten()
            
        obj[0] = ((0, 'START', 0, obj._last_pnt_marker), ics)
        
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        
        self._last_pnt_marker = getattr(obj, '_last_pnt_marker', None)
        self.n_ics = getattr(obj, 'n_ics', None)
        
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
            result[ind] = stamp + tuple(vec)
        
        super(PLog, self).__setitem__(i, result)
    
#%%
    
if __name__ == '__main__':
    from ensemble import Ensemble
    from particle import Particle
    
    NRec = 5
    
    bunch = Ensemble.populate(Particle(), Sz=1, dK=(0, 1e-4, 1), x=(-1e-3, 1e-3, 1))    
    log = PLog(bunch.ics, NRec)
    
    #%%
    l12 = [1,2,3,4,5,6,7,8,9,10,11,12]
    statevec = np.array(l12+l12)
    metadata = (1,'RF',0,-1)
    
    for i in range(1,NRec):
        metadata = (i, 'RF',0,-1)
        log[i] = (metadata, statevec*i)