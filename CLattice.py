#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 09:19:48 2017

@author: alexa
"""
import PyDSTool as DST
import collections as CLN
import CParticle as PCL
#import pathos.multiprocessing as MLP
from utilFunc import phi, sadd, smult

from subprocess import call
import re
import copy

class Lattice:
    
    def __init__(self, ElSeq, RefPart = PCL.Particle(), Options=dict()):

        try: Gen = Options['Generator']
        except KeyError: Gen = 'dopri'
        
        Gen = Gen.upper()
        if Gen == 'VODE': tlang = 'python'
        else: tlang = 'c'
        
        call('find ./dop853_temp/ -type f -exec rm {} + ', shell=True) # clear compiled code
        call('find ./radau5_temp/ -type f -exec rm {} + ', shell=True)
        
        self.fCount = len(ElSeq)
        self.fPardict = dict()
        
        #%% global event definitions
        try: passed_events = Options['Events']
        except KeyError:
            passed_events = []
        else: 
            if type(passed_events) is not list:  
                passed_events = list(passed_events)
        
        event_args = {'name':'NaN_event','eventtol':1e-4,'eventdelay':0,'term':True, 'active':True, 'precise':True}
        
        ## the NaN error handling event definition
        self.fPardict.update({'offset':10000}) # require Ps**2 > 100**2
        NaN_event = DST.makeZeroCrossEvent('pow(Pc(dK),2)-pow(Pc(0),2)*(pow(px,2) + pow(py,2)) - offset',-1,
                                           event_args,varnames=['dK','px','py'], targetlang=tlang,
                                           fnspecs=RefPart.fFndict,parnames=['Mass0','KinEn0','offset'])

        
        #%% constructing a DS for each element
        DS_list=list()
        MI_list = list()
        _id=0
        self.__fLength = 0 #lattice length
        for e in ElSeq:
            self.__fLength += e.fPardict['Length']
            
            ## RHS for DS 
            DSargs = Lattice.setup_element(e,RefPart)
            
            try: passed_algparams = Options['Algparams']
            except KeyError:
                passed_algparams = dict()
            else: 
                if not type(passed_algparams) == dict:
                    raise TypeError('Expected a dictionary')
            
            DSargs.algparams = passed_algparams
            
            ## model selection
            DSargs.update({'xdomain':{'start':_id}}) #can't select the initial model w/o this
            DSargs.xtype={'start':DST.int} # ensure integer type
            DSargs.varspecs.update({'start': '0'}) 
            
             ## events
            _id +=1
            event_args.update({'name':'passto'+str(_id%self.fCount)})
            pass_event = DST.makeZeroCrossEvent('s-'+str(e.fPardict['Length']),1,event_args,varnames=['s'],parnames=list(self.fPardict.keys()),targetlang=tlang)
            DSargs.events = [pass_event, NaN_event]+passed_events
            
            DSargs.pars.update(self.fPardict)
            if Gen == 'DOPRI': DS = DST.Generator.Dopri_ODEsystem(DSargs)
            elif Gen == 'RADAU': DS = DST.Generator.Radau_ODEsystem(DSargs)
            else: DS = DST.Generator.Vode_ODEsystem(DSargs)
            DS.Particle = RefPart
            DS.Element = e
            DS_list.append(DS)
            DS = DST.embed(DS,name=e.fName)
            MI_list.append(DST.intModelInterface(DS))
            
        self.fMods = DS_list
        
        #%% construction of the hybrid model
        try: passed_tmaps = Options['T-maps']
        except KeyError:
            passed_tmaps = dict()
        else:
            if not type(passed_tmaps) == dict:
                raise TypeError('Expected a dictionaty of dicts: element index: {x:map, y:map, ...}')
        
        all_names = [e.fName for e in ElSeq]
        info = list()
        for i in range(len(MI_list)):
            ## transition event mapping dictionary
            transdict = {'dK':"self.outin([x,y,ts,px,py,dK])"} # this is frontkick_n+1(backkick_n(state))
            transdict.update({'s':'0','at':'(at+1)%'+str(self.fCount)}) # then reset s for the next element
            # from Options
            try: tmap = passed_tmaps[i]
            except KeyError: pass
            else: transdict.update(tmap)
            # transition event state map
            tem = DST.EvMapping(transdict, model=DS_list[i]) 
            # transition from element i to element i+1
            tem.outin = lambda state: DS_list[(i+1)%self.fCount].Element.frontKick(DS_list[i%self.fCount].Element.rearKick(state))[5] # dK is state[5]
            info.append(DST.makeModelInfoEntry(MI_list[i],all_names,[('passto'+str((i+1)%self.fCount),(MI_list[(i+1)%self.fCount].model.name, tem))]))
            
        modelInfoDict = DST.makeModelInfo(info)
        
        mod_args = {'name':'lattice','modelInfo':modelInfoDict}
        self.fDSModel = DST.Model.HybridModel(mod_args)
    
    @classmethod
    def setup_element(cls, Element, RefPart):
        
        ## definitions
        arg = Element.fArgStr
        
        # fields
        sExA = 'Ex'+arg
        sEyA = 'Ey'+arg
        sEsA = 'Es'+arg
        sBxA = 'Bx'+arg
        sByA = 'By'+arg
        sBsA = 'Bs'+arg

        # spin-related
        t6 =  'v4_Tp* (q / (v0_Lgamma * m0 * m0* clight * clight)) * (G + 1/(1 + v0_Lgamma))'
        sp1 = 'v4_Tp*(-q / (v0_Lgamma*m0))*(1 + G * v0_Lgamma)'
        sp2 = 'v4_Tp*( q / (v0_Lgamma* m0*m0*m0*clight*clight)) * (G/(1 + v0_Lgamma))*(v2_Px*v0_Bx+v2_Py*v0_By+v3_Ps*v0_Bs)'
        
        ## derivative definitions
        sXpA = 'v0_hs*v0_P0c*px/v2_Ps'
        sYpA = 'v0_hs*v0_P0c*py/v2_Ps'    
        sHpA = 'v0_hs*v0_Pc/v2_Ps'
        sTpA = 'v3_Hp/v1_V'    
        
        det = lambda a,b,c,d: phi('-',smult(a,b),smult(c,d))
        sFxATp = '1e-6*clight*'+ sadd(smult(sExA,sTpA), det(sYpA,sBsA,'1',sByA)) 
        sFyATp = '1e-6*clight*'+ sadd(smult(sEyA,sTpA), det('1',sBxA,sXpA,sBsA)) 
        
        # these are probably from TBMT
        Sxp =      'Curve * Ss + v5_t6 * ((v3_Ps * v0_Ex - v2_Px * v0_Es) * Ss - (v2_Px * v0_Ey - v2_Py * v0_Ex) * Sy) + (v5_sp1*v0_By+v5_sp2*v2_Py)*Ss-(v5_sp1*v0_Bs+v5_sp2*v3_Ps)*Sy'
        Syp =                   'v5_t6 * ((v2_Px * v0_Ey - v2_Py * v0_Ex) * Sx - (v2_Py * v0_Es - v3_Ps * v0_Ey) * Ss) + (v5_sp1*v0_Bs+v5_sp2*v3_Ps)*Sx-(v5_sp1*v0_Bx+v5_sp2*v2_Px)*Ss'
        Ssp = '(-1)*Curve * Sx + v5_t6 * ((v2_Py * v0_Es - v3_Ps * v0_Ey) * Sy - (v3_Ps * v0_Ex - v2_Px * v0_Es) * Sx) + (v5_sp1*v0_Bx+v5_sp2*v2_Px)*Sy-(v5_sp1*v0_By+v5_sp2*v2_Py)*Sx'
        
        ## 
        reuse = RefPart.fReuse
        reuse.update({
                    '(1+Curve*x)':'v0_hs',
                    sExA:'v0_Ex', sEyA:'v0_Ey', sEsA:'v0_Es',
                    sBxA:'v0_Bx', sByA:'v0_By', sBsA:'v0_Bs',
                    sXpA:'v3_Xp',sYpA:'v3_Yp', sHpA:'v3_Hp',
                    sTpA:'v4_Tp',
                    t6:'v5_t6', sp1:'v5_sp1', sp2:'v5_sp2',
                    Sxp:'v6_Sxp', Syp:'v6_Syp', Ssp:'v6_Ssp'
                 })
        
        RHS = {
                'x'  : 'v3_Xp',
                'y'  : 'v3_Yp',  
                'H'  : 'v3_Hp',
                'ts' : 'v4_Tp',
                'px' : sadd(sFxATp, 'Curve*v2_Ps')+'/v0_P0c',
                'py' : sFyATp+'/v0_P0c',
                'dK' : sadd(smult(sExA,sXpA), smult(sEyA,sYpA), sEsA) + '*1e-6/KinEn0',
                'Sx' : 'v6_Sxp',
                'Sy' : 'v6_Syp',
                'Ss' : 'v6_Ssp',
                's'  : '1',
                'start':'0',
                'at':'0'
                }
        
        DSargs = DST.args(name=Element.fName)       
        DSargs.pars = {**RefPart.fPardict, **Element.fPardict}
        DSargs.fnspecs = {**RefPart.fFndict, **Element.fFndict}
        DSargs.reuseterms=reuse
        
        DSargs.varspecs = RHS      
        
        return DSargs

    def getLength(self):
        return self.__fLength
    
    def listModelNames(self):
        namelist = list()
        for model in self.fDSModel.sub_models():
            namelist.append(model.name)
            
        return namelist
    
    def __compute(self, ArgDict): # for parallel code
        tdata = ArgDict['tdata']
        inistate = ArgDict['inistate']
        name = ArgDict['name']
        print(name)
        self.fDSModel.compute(name,ics=inistate,tdata=tdata)
        return self.fDSModel
        
    
    def track(self, Ensemble, NTurns, StartID = '0'):
        tstp = NTurns * self.__fLength
        self.fDSModel.tdata = [0,tstp]
        
        ## for recurring computation
        # somehow i can't do self.fDSmodel.cleanupMemory(), because the objects
        # involved don't have some required attributes. Could be a python version-related issue
        trajnames = self.fDSModel.trajectories.keys()
        for name in Ensemble.fIniStateDict.keys():
            if name in trajnames: 
                try:
                    key, ind = re.split('_',name)
                    ind = str(int(ind) + 1)
                except ValueError:
                    key = re.split('_',name)[0]
                    ind = str(1)
                    
                key += '_'+ind
                Ensemble.rename(name, key)
        
        
        #%% parallel computation 
        # would be nice to fix it
        # however  unlikely, due to lack of knowledge of PyDSTool's internals
        
#        p = MLP.Pool(3)
#        arg = list()
#        for name,inistate in Ensemble.fStateDict.items():
#            inistate.update({'start':StartID})
#            arg.append({'name':name, 'inistate':inistate,'tdata':[0,tstp]})
#        
#        val = p.starmap(self.__compute, zip(arg))
#        p.close()
#        
#        self.fDSModel = val
    
        #%% sequential computation 
        # this works as expected
            
        IniStateDictCopy = copy.deepcopy(Ensemble.fIniStateDict)
        for name,inistate in IniStateDictCopy.items():
            inistate.update({'start':StartID, 'at':StartID})
            names = list(inistate.keys())
            linistate = [inistate[e] for e in names]
            dK = self.fMods[int(StartID)].Element.frontKick(linistate)[5] # state[5] = dK
            inistate.update({'dK':dK})
            self.fDSModel.compute(name,ics=inistate)
            
        Ensemble.fTrajectories = self.fDSModel.trajectories