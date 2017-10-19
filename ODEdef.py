#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 14:44:21 2017

@author: alexa

Element and Particle classes set up DS arguments upon construction;
then inside the Lattice class' constructor, ODEs are initialized with 
those DS-args combined, and put into a hybrid model of the lattice.

"""

import PyDSTool as DST
from matplotlib import pyplot as PLT
import re

class Particle:
    
    def __init__(self,Mass0 = 1876.5592, KinEn0 = 270.11275, G = -.142987):
        q = 1.602176462e-19
        clight = 2.99792458e8
        self.pardict = {'Mass0':Mass0, 'KinEn0':KinEn0, 'G':G,
                        'q':q, 'clight':clight,
                        'm0': q*1e6*Mass0/clight**2}
        
        self.fndict = {'KinEn':(['dK'],'KinEn0*(1+dK)'), 
                       'Lgamma':(['dK'],'KinEn(dK)/Mass0 + 1'),
                       'Lbeta':(['dK'],'sqrt(pow(Lgamma(dK),2)-1)/Lgamma(dK)'),
                       'Pc':(['dK'],'sqrt(pow(Mass0 + KinEn(dK),2) - pow(Mass0,2))')
                       }
        
        self.reuse = {'Pc(0)':'v0_P0c','Pc(dK)':'v0_Pc', 
                      'Lbeta(dK)':'v0_Lbeta', 'Lgamma(dK)':'v0_Lgamma', 'KinEn(dK)':'v0_KinEn',
                      'v0_Lbeta*clight':'v1_V',
                      'v0_P0c*px':'v1_Px', 'v0_P0c*py':'v1_Py',
                      'sqrt(pow(v0_Pc,2)- pow(v1_Px,2) - pow(v1_Py,2))':'v2_Ps',
                      'v1_V*v1_Px/v0_Pc':'v2_Vx', 'v1_V*v1_Py/v0_Pc':'v2_Vy',
                      'v1_V*v2_Ps/v0_Pc':'v3_Vs'
                    }
        
        self.defs = dict()
        ell = list(self.reuse.items())
        for key,value in ell:
            self.defs.update({re.sub('.*_','',value):key})
            
            

class Element:
    
    fArgList = ['x','y','ts','px','py','dK','H','s','start','Sx','Sy','Ss']
    
    fArgStr = None
    
    fCount = 0
    
    def __init__(self, Curve, Length, Name = "Element"): 
        self.fName = Name+"_"+str(Element.fCount)
        
        Element.fCount += 1
        
        self.pardict = {'Curve':Curve, 'Length':Length}
        
        self.fndict = { # defines the element EM field
                'Ex':(self.fArgList, '0'),
                'Ey':(self.fArgList, '0'),
                'Es':(self.fArgList, '0'),
                'Bx':(self.fArgList, '0'),
                'By':(self.fArgList, '.46'),
                'Bs':(self.fArgList, '0')
                }
        
        arg = ''
        for e in self.fArgList: arg = arg+','+e
        arg = '('+arg[1:len(arg)]+')'
        self.fArgStr = arg # argument string '(x,y,...)' for RHS definition
        
    def frontKick(self, state, particle):
        return list(state)
    
    def rearKick(self, state, particle):
        return list(state)
        

class Lattice:
# inside Lattice class constructor,
# give a list of elements and a particle

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
    
    def __setup_element(self, Element,RefPart):

        ## definitions
        arg = Element.fArgStr
#        larg = Element.fArgList
        defs = RefPart.defs
        
        # fields
        sExA = 'Ex'+arg
        sEyA = 'Ey'+arg
        sEsA = 'Es'+arg
        sBxA = 'Bx'+arg
        sByA = 'By'+arg
        sBsA = 'Bs'+arg
        
        # v cross B
        sVxBx = defs['Vy']+'*'+sBsA+'-'+sByA+'*'+defs['Vs']
        sVxBy = defs['Vs']+'*'+sBxA+'-'+sBsA+'*'+defs['Vx']
        
        # Lorentz forces
        sFxA = 'q*('+sExA+'+'+sVxBx+')'
        sFyA = 'q*('+sEyA+'+'+sVxBy+')'
        
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
                'px' : '('+sFxA+'*v4_Tp'+ '+ Curve*v2_Ps)/v0_P0c',
                'py' : '('+sFyA+'*v4_Tp)/v0_P0c',
                'dK' : '('+sExA+'*'+sXpA+'+'+sEyA+'*'+sYpA+'+'+sEsA+')*1e-6/KinEn0',
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


#%%

if __name__ == '__main__':
    
    p = Particle()
    e1 = Element(0,5)
    e2 = Element(0,2.5)
    
    state = [0,0,0,
             1e-3,0,0,
             0,0,0,
             0,0,1]
    names = e1.fArgList
    icdict = dict(zip(names,state))
#%%    

Lat = Lattice([e1],p)
    
Lat.fDSModel.compute('test',ics=icdict,tdata=[0,100])

#%%
pts = Lat.fDSModel.sample('test')
PLT.plot(pts['t'],pts['x'],label='x')
PLT.plot(pts['t'],pts['Sx'],label='Sx')
PLT.legend()
