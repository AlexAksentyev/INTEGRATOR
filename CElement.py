import numpy as NP
import collections as CLN
from utilFunc import phi
import copy
import re


class Element:
    
    fArgList = ['x','y','ts','px','py','dK','H','s','Sx','Sy','Ss']
    
    fArgStr = None
    
    def __init__(self, Curve, Length, Name = "Element"):
        self.fName = Name
        
        self.fPardict = {'Curve':Curve, 'Length':Length}
        
        self.fGeomdict = copy.deepcopy(self.fPardict)
        
        self.fFndict = { # defines the element EM field
                'Ex':(self.fArgList, '0'),
                'Ey':(self.fArgList, '0'),
                'Es':(self.fArgList, '0'),
                'Bx':(self.fArgList, '0'),
                'By':(self.fArgList, '0'),
                'Bs':(self.fArgList, '0')
                }
        
        self.fArgStr = phi(',',*self.fArgList) # argument string '(x,y,...)' for RHS definition
        
    def __repr__(self):
        return str(self.getField())
    
    def setField(self, **kwargs):
        self.__setField(dict(**kwargs))
        
    def __setField(self, FldDict): # field definitions given in the form name:definition, w/o signature
        inFldDict = {key:(self.fArgList, value) for key, value in FldDict.items()}
        self.fFndict.update(inFldDict)
        
    def getField(self, fulldef = False):
        if fulldef: return self.fFndict
        else: return {key:value[1] for key,value in self.fFndict.items()}

    def getGeometry(self):
        return self.fGeomdict

    def frontKick(self, state, particle):
        return list(state)
    
    def rearKick(self, state, particle):
        return list(state)
    
class HasCounter:
    fCount = 0
    
    __fSep = "_"
    
    def __init__(self):
        if 'fName' not in self.__dict__.keys():
            self.fName = 'NoName'    
        self.fName += self.__fSep+str(self.__class__.fCount)
        self.__class__.fCount += 1

    def copy(self, Name = None):
        self.__class__.fCount += 1
        res = copy.deepcopy(self)
        if Name is None: res.fName = re.sub(self.__fSep+'.*',self.__fSep+str(res.__class__.fCount-1),res.fName)
        else: res.fName = Name
        return res

class Drift(Element, HasCounter):
    """ drift space
    """
    
    def __init__(self, Length, Name = "Drift"):
        Element.__init__(self, 0, Length, Name)
        HasCounter.__init__(self)

class MQuad(Element, HasCounter):
    """ magnetic quadrupole
    """

    def __init__(self, Length, Grad, Name = "MQuad"):
        Element.__init__(self, 0, Length, Name)
        HasCounter.__init__(self)
        self.setGrad(Grad)   
        
    def setGrad(self, value):
        self._Element__setField({'Bx':str(value)+'*(-y)','By':str(value)+'*(-x)'})
        self.__fGrad = value
        self.fGeomdict.update({'Grad':value})
        
    def getGrad(self):
        return self.__fGrad
        

class MDipole(Element, HasCounter):
    """ bending magnetic dipole (horizontally bending);
    define _BField as a tuple;
    could also be used as a solenoid, if _BField = (0,0,Bs)
    and fCurve = 0
    """
    
    def __init__(self, Length, R, BField, Name = "MDipole"):
        Element.__init__(self, 1/R, Length, Name)
        HasCounter.__init__(self)
        self.fGeomdict.update({'R':R})
        self.setBField(BField)
        
    def setBField(self,BField):
        if not isinstance(BField, CLN.Sequence): self.__fBField = (0,BField,0) # by default B = By
        else: self.__fBField = BField[:]
        BField = dict(zip(['Bx','By','Bs'],[str(b) for b in self.__fBField]))
        self._Element__setField(BField)

    def getBField(self):
        return self.__fBField
    
    @classmethod    
    def computeBStrength(cls,particle, R):
        return particle.Pc(particle.fKinEn0)*1e6/(R*particle.CLIGHT())
    
    @classmethod
    def computeRadius(cls,particle, BField): 
        return particle.Pc(particle.fKinEn0)*1e6/(BField*particle.CLIGHT())
  

class Solenoid(MDipole, HasCounter):
    
    def __init__(self, Length, Bs, Name = "Solenoid"):
        Element.__init__(self, 0, Length, Name)
        HasCounter.__init__(self)
        self.setBField((0,0,Bs))
        
    @classmethod
    def computeBStrength(cls, *args):
        print('Not defined')
        
    @classmethod
    def computeRadius(cls, *args):
        print('Infinite radius = zero curvature')
        

class MSext(Element, HasCounter):
    """ magnetic sextupole
    """
    
    def __init__(self, Length, Grad, Name = "MSext"):
        Element.__init__(self, 0, Length, Name)
        HasCounter.__init__(self)
        self.setGrad(Grad)
        
    def setGrad(self, value):
        self._Element__setField({'Bx':str(value)+'*(x*y)','By':str(value)+'*(x*x - y*y)','Bs':'0'})
        self.__fGrad = value
        self.fGeomdict.update({'Grad':value})
        
    def getGrad(self):
        return self.__fGrad
        
class Wien(Element, HasCounter):
    """ wien filter ?!?! 
    The magnetic field definition is weird
    """
    
    def __init__(self, Length, Hgap, RefPart, EField=None, BField=None, R = None, Name = "Wien?"):
        
        assert any([e is not None for e in [EField, R]]), "Either the E-field, or Radius must be defined"
        
        if R is None: Crv = None
        else: 
            Crv = 1/R
            self.__fR = [R, R - Hgap, R**2/(R-Hgap)]
            
        self.__fHgap = Hgap
        
        Element.__init__(self, Crv, Length, Name)
        HasCounter.__init__(self)
        self.fGeomdict.update({'Hgap':Hgap})
        self.fPardict.update(RefPart.fPardict)
        
        self.setEField(EField)
        self.setBField(BField)
        
    def __GammaBeta(self):
        Mass0 = self.fPardict['Mass0']
        K0 = self.fPardict['KinEn0']
        gamma = K0 / Mass0 + 1
        beta = NP.sqrt(gamma**2-1)/gamma
        return (gamma, beta)
    
    def __Pc(self, KNRG):
        return NP.sqrt((self.fPardict['Mass0'] + KNRG)**2 - self.fPardict['Mass0']**2)
    
    def setEField(self, EField=None):
        P0c = self.__Pc(self.fPardict['KinEn0'])
        gamma, beta = self.__GammaBeta()
        if EField is None:
            R = self.__fR
            EField = - P0c*beta/R[0] * 1e6
        else:
            R = - P0c*beta/EField * 1e6
            self.__fR = [R, R - self.__fHgap, R**2/(R-self.__fHgap)]
            R = self.__fR
        
        self.__fEField = (EField,0,0)
        self.__fVolt = (EField * R[0] * NP.log(R[2] / R[1])) / (-2)
        self._Element__setField({'Ex':str(EField)})
        self.fPardict.update({'Curve':1/R[0]})
        self.fGeomdict.update({'R':R[0], 'Curve':1/R[0]})
        
    def setBField(self, BField=None):        
        e0 = self.fPardict['q']
        m0 = self.fPardict['m0']
        clight = self.fPardict['clight']
        qm = e0/(m0*clight**2)
        
        gamma, beta = self.__GammaBeta()
        k = 1.18 * qm * ((2 - 3*beta**2 - .75*beta**4)/beta**2 + 1/gamma)
        v = beta*clight
        
        if BField is None:
            if self.__fEField is None: self.setEField()
            BField = -self.__fEField[0]/v
            
        h = 1/self.__fR[0]
        
        B = phi('*', 0.018935,-BField,phi('+',-h, phi('*',k,BField,v)), 'x')
        B = B[1:len(B)-1]
        
        self._Element__setField({'By':B})
    
    def frontKick(self, state, particle):
        x=state[0]
        Xk = list(state)
        R = self.__fR[0]
        R1 = self.__fR[1]
        R2 = self.__fR[2]
        V = self.__fVolt
        KinEn0 = particle.fPardict['KinEn0']
        u = -V + 2*V*NP.log((R+x)/R1)/NP.log(R2/R1)
        Xk[5] -= u*1e-6/KinEn0
        return Xk
        
    def rearKick(self, state, particle):
        x=state[0]
        Xk = list(state)
        R = self.__fR[0]
        R1 = self.__fR[1]
        R2 = self.__fR[2]
        V = self.__fVolt
        KinEn0 = particle.fPardict['KinEn0']
        u = -V + 2*V*NP.log((R+x)/R1)/NP.log(R2/R1)
        Xk[5] += u*1e-6/KinEn0
        return Xk
