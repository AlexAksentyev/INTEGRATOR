#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 09:19:48 2017

@author: alexa
"""
import PyDSTool as DST

class Lattice:
    
    fDSModel = None
    
    def __init__(self, ElementSeq, RefPart):
        
        #%% logging in the particle data
        pardict = {'Mass0':RefPart.fMass0, 'Kin0':RefPart.fKinEn0, 'P0c':RefPart.Pc(RefPart.fKinEn0)}
        fndict = {'KinEn':(['dK'],'Kin0*(1+dK)'), 
                  'Lgamma':(['dK'],'KinEn(dK)/Mass0 + 1'),
                  'Lbeta':(['dK'],'sqrt(pow(Lgamma(dK),2)-1)/Lgamma(dK)'),
                  'Pc':(['dK'],'sqrt(pow(Mass0 + KinEn(dK),2) - pow(Mass0,2))')}
        
        #%% logging in the element data
        size = len(ElementSeq)
        
        
        at=0
        for element in ElementSeq:
            at += element.fLength
            pardict.update({'L'+element.fName:at}) # log in the element position along the optical axis
            pardict.update({'kappa'+element.fName:element.fCurve})
        
        pardict.update({'Ltot':at}) # the total lattice length
        
        DSargs = DST.args(tdata=[0, pardict['Ltot']]) #in PyDSTool t is the independent variable; t is equiv to s
        
        #%% global event definitions
        event_args = {'name':'NaN_event','eventtol':1e-4,'eventdelay':0,'term':True, 'active':True}
        
        ## the NaN error handling event definition
        pardict.update({'offset':10000}) # require Ps**2 > 100**2
        NaN_event = DST.makeZeroCrossEvent('pow(Pc(dK),2)-pow(P0c,2)*(pow(px,2) + pow(py,2)) - offset',-1,
                                           event_args,varnames=['dK','px','py'],
                                           fnspecs=fndict,parnames=['Mass0','P0c','Kin0','offset'])
       
        #%% constructing the differential systems for each element
        # terminal element events are defined here
        
        DSargs.pars = pardict
        DSargs.fnspecs = fndict
        
        DSargs.ignorespecial = ['state','xp','yp','tp','pxp','pyp','dKp','Sxp','Syp','Ssp','Hp']
        DSargs.vfcodeinsert_start = """state = [x,y,ts,px,py,dK,Sx,Sy,Ss,H]
    xp,yp,tp,pxp,pyp,dKp,Sxp,Syp,Ssp,Hp = ds.Particle.RHS(state, [], ds.Element)
        """
        
        ## the derivatives
        DSargs.varspecs = {'x': 'xp', 'y': 'yp', 's':'1',
                           'ts':'tp', 'H':'Hp', 'start':'0',
                           'dK':'dKp', 'px':'pxp', 'py':'pyp', 
                           'Sx':'Sxp', 'Sy':'Syp', 'Ss':'Ssp'}
        
        ModList = list()
        MI_list = list()
        _id=0
        for element in ElementSeq:
            DSargs.update({'name':element.fName})
            DSargs.update({'xdomain':{'start':_id}}) #this is actually a very important string here, for initial model selection!
            DSargs.xtype={'start':DST.int}
            DSargs.varspecs.update({'start': str(_id)})
            DSargs.fnspecs.update({'hs':(['x'],'1 + x*kappa'+element.fName)})
            _id +=1
            event_args.update({'name':'passto'+str(_id%size)})
            if _id%size != 0:
                pass_event = DST.makeZeroCrossEvent('s%Ltot-L'+element.fName,1,event_args,varnames=['s'],parnames=list(pardict.keys()))
            else: 
                pass_event = DST.makeZeroCrossEvent('s-Ltot*ceil(s/(Ltot*1.001))',1,event_args,varnames=['s'],parnames=list(pardict.keys())) #the factor at Ltot must be > eventtol
            DSargs.events = [pass_event, NaN_event]
            DS = DST.Vode_ODEsystem(DSargs)
            DS.Element = element
            DS.Particle = RefPart
            DS = DST.embed(DS,name=element.fName)
            ModList.append(DS)
            MI_list.append(DST.intModelInterface(DS))
        
        all_names = [e.fName for e in ElementSeq]
        info = list()
        
        for i in range(len(MI_list)):
            info.append(DST.makeModelInfoEntry(MI_list[i],all_names,[('passto'+str((i+1)%size),MI_list[(i+1)%size].model.name)]))    
        
        modelInfoDict = DST.makeModelInfo(info)
        
        mod_args = {'name':'lattice','modelInfo':modelInfoDict}
        self.fDSModel = DST.Model.HybridModel(mod_args)
