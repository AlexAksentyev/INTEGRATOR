#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 09:19:48 2017

@author: alexa
"""
import PyDSTool as DST
#import pathos.multiprocessing as MLP
from utilFunc import phi

class Lattice:
    
    def __init__(self, ElSeq, RefPart):
        
        self.fCount = len(ElSeq)
        self.pardict = dict()
        
        #%% global event definitions
        event_args = {'name':'NaN_event','eventtol':1e-4,'eventdelay':0,'term':True, 'active':True, 'precise':True}
        
        ## the NaN error handling event definition
        self.pardict.update({'offset':10000}) # require Ps**2 > 100**2
        NaN_event = DST.makeZeroCrossEvent('pow(Pc(dK),2)-pow(Pc(0),2)*(pow(px,2) + pow(py,2)) - offset',-1,
                                           event_args,varnames=['dK','px','py'], targetlang='c',
                                           fnspecs=RefPart.fndict,parnames=['Mass0','KinEn0','offset'])

        
        #%% constructing a DS for each element
        DS_list=list()
        MI_list = list()
        _id=0
        self.__fLength = 0 #lattice length
        for e in ElSeq:
            self.__fLength += e.pardict['Length']
            self.pardict.update({'L'+e.fName:e.pardict['Length']}) # log in the element position along the optical axis
            self.pardict.update({'kappa'+e.fName:e.pardict['Curve']}) # and its curvature
            
            ## RHS for DS 
            DSargs = self.__setup_element(e,RefPart)
            
            ## model selection
            DSargs.update({'xdomain':{'start':_id}}) #can't select the initial model w/o this
            DSargs.xtype={'start':DST.int} # ensure integer type
            DSargs.varspecs.update({'start': '0'}) 
            
             ## events
            _id +=1
            event_args.update({'name':'passto'+str(_id%self.fCount)})
            pass_event = DST.makeZeroCrossEvent('s-'+str(e.pardict['Length']),1,event_args,varnames=['s'],parnames=list(self.pardict.keys()),targetlang='c')
            DSargs.events = [pass_event, NaN_event]
            
            DSargs.pars.update(self.pardict)
            DS = DST.Generator.Dopri_ODEsystem(DSargs)
            DS.Particle = RefPart
            DS.Element = e
            DS_list.append(DS)
            DS = DST.embed(DS,name=e.fName)
            MI_list.append(DST.intModelInterface(DS))
            
        self.fMods = DS_list
        
        #%% construction of the hybrid model
        all_names = [e.fName for e in ElSeq]
        info = list()
        for i in range(len(MI_list)):
            transdict = {'dK':"self.outin([x,y,ts,px,py,dK],self.RefPart)"} # this is frontkick_n+1(backkick_n(state))
            transdict.update({'s':'0'}) # then reset s for the next element
            epmapping = DST.EvMapping(transdict, model=MI_list[i].model) # transition event state map
            epmapping.outin = lambda state, part: DS_list[(i+1)%self.fCount].Element.frontKick(DS_list[i%self.fCount].Element.rearKick(state,part),part)[5] # dK is state[5]
            epmapping.RefPart = RefPart
            info.append(DST.makeModelInfoEntry(MI_list[i],all_names,[('passto'+str((i+1)%self.fCount),(MI_list[(i+1)%self.fCount].model.name, epmapping))]))
        
        modelInfoDict = DST.makeModelInfo(info)
        
        mod_args = {'name':'lattice','modelInfo':modelInfoDict}
        self.fDSModel = DST.Model.HybridModel(mod_args)
    
    def __setup_element(self, Element, RefPart):

        sadd = lambda a,b: phi('+',a,b)
        smult = lambda a,b: phi('*',a,b)
        
        ## definitions
        arg = Element.fArgStr
        defs = RefPart.defs
        
        # fields
        sExA = 'Ex'+arg
        sEyA = 'Ey'+arg
        sEsA = 'Es'+arg
        sBxA = 'Bx'+arg
        sByA = 'By'+arg
        sBsA = 'Bs'+arg
        
        # v cross B
        det = lambda a,b,c,d: phi('-',smult(a,b),smult(c,d))
        sVxBx = det(defs['Vy'],sBsA,sByA,defs['Vs'])
        sVxBy = det(defs['Vs'],sBxA,sBsA,defs['Vx'])
        
        # Lorentz forces
        
        sFxA = 'q*'+ sadd(sExA,sVxBx)
        sFyA = 'q*'+ sadd(sEyA,sVxBy)
        
        # spin-related
        t6 =  'v4_Tp* (q / (v0_Lgamma * m0 * m0* clight * clight)) * (G + 1/(1 + v0_Lgamma))'
        sp1 = 'v4_Tp*(-q / (v0_Lgamma*m0))*(1 + G * v0_Lgamma)'
        sp2 = 'v4_Tp*( q / (v0_Lgamma* m0*m0*m0*clight*clight)) * (G/(1 + v0_Lgamma))*(v1_Px*v0_Bx+v1_Py*v0_By+v2_Ps*v0_Bs)'
        
        ## derivative definitions
        sXpA = 'v0_hs*v0_P0c*px/v2_Ps'
        sYpA = 'v0_hs*v0_P0c*py/v2_Ps'    
        sHpA = 'v0_hs*v0_Pc/v2_Ps'
        sTpA = 'v3_Hp/v1_V'    
        
        # these are probably from TBMT
        Sxp =      'Curve * Ss + v5_t6 * ((v2_Ps * v0_Ex - v1_Px * v0_Es) * Ss - (v1_Px * v0_Ey - v1_Py * v0_Ex) * Sy) + (v5_sp1*v0_By+v5_sp2*v1_Py)*Ss-(v5_sp1*v0_Bs+v5_sp2*v2_Ps)*Sy'
        Syp =                   'v5_t6 * ((v1_Px * v0_Ey - v1_Py * v0_Ex) * Sx - (v1_Py * v0_Es - v2_Ps * v0_Ey) * Ss) + (v5_sp1*v0_Bs+v5_sp2*v2_Ps)*Sx-(v5_sp1*v0_Bx+v5_sp2*v1_Px)*Ss'
        Ssp = '(-1)*Curve * Sx + v5_t6 * ((v1_Py * v0_Es - v2_Ps * v0_Ey) * Sy - (v2_Ps * v0_Ex - v1_Px * v0_Es) * Sx) + (v5_sp1*v0_Bx+v5_sp2*v1_Px)*Sy-(v5_sp1*v0_By+v5_sp2*v1_Py)*Sx'
        
        ## 
        reuse = RefPart.reuse
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
                'px' : sadd(smult(sFxA,'v4_Tp'), 'Curve*v2_Ps')+'/v0_P0c',
                'py' : smult(sFyA,'v4_Tp')+'/v0_P0c',
                'dK' : sadd(smult(sExA,sXpA), smult(sEyA,sYpA), sEsA) + '*1e-6/KinEn0',
                'Sx' : 'v6_Sxp',
                'Sy' : 'v6_Syp',
                'Ss' : 'v6_Ssp',
                's'  : '1',
                'start':'0'
                }
        
        DSargs = DST.args(name='ODEs')       
        DSargs.pars = {**RefPart.pardict, **Element.pardict}
        DSargs.fnspecs = {**RefPart.fndict, **Element.fndict}
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
            
        for name,inistate in Ensemble.fIniStateDict.items():
            inistate.update({'start':StartID})
            self.fDSModel.compute(name,ics=inistate)
            
        Ensemble.fTrajectories = self.fDSModel.trajectories
            
