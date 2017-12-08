#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:53:58 2017

@author: alexa
"""

import numpy as NP

#import RHS

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
    
#class Writer:
#    """ Unfinished
#    """
#    def __init__(self, Ensemble, Lattice, file_handle, Log_type):
#        """ Lattice argument b/c in any case I will probably
#            want to write some lattice data as well,
#            like element tilts, parameters
#        """
#        self.Log_type = Log_type
#        self.file_handle = file_handle
#        self.Ensemble = Ensemble
#        self.Lattice = Lattice
#        
#        self.__setup_file()
#        self.__create_tables()
#        
#        self.__fLastPnt = -1 # marker of the state vector upon exiting an element
#        
#    def __setup_file(self):
#        ## setting up file
#        # write particle parameters (Mass0, KinEn0, G)
#        ppar = self.Ensemble.Particle.getParams() 
#        self.file_handle.create_table('/','Particle',ppar)
#        # separate group to write p-logs in
#        self.file_handle.create_group('/','Logs', 'Particle logs')
#        
#    def __create_tables(self):
#        ## creating data tables to fill
#        for p in self.Ensemble: 
#            self.file_handle.create_table(self.file_handle.root.Logs, 
#                                     'P'+str(p.PID), self.Log_type)
#    
#    def write_log(self, from_to):
#        old_ind, ind = from_to
#        # write data to file
#        for p in self.Ensemble:
#            tbl = getattr(self.file_handle.root.Logs, 'P'+str(p.PID))
#            tbl.append(p.Log[old_ind:ind])
#            tbl.flush()