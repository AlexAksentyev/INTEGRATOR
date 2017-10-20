import numpy as NP
import collections as CLN
from utilFunc import phi



class Element:
    
    fArgList = ['x','y','ts','px','py','dK','H','s','start','Sx','Sy','Ss']
    
    fArgStr = None
    
    def __init__(self, Curve, Length, Name = "Element"):
        self.fName = Name
        
        self.pardict = {'Curve':Curve, 'Length':Length}
        
        self.fndict = { # defines the element EM field
                'Ex':(self.fArgList, '0'),
                'Ey':(self.fArgList, '0'),
                'Es':(self.fArgList, '0'),
                'Bx':(self.fArgList, '0'),
                'By':(self.fArgList, '0'),
                'Bs':(self.fArgList, '0')
                }
        
        self.fArgStr = phi(',',*self.fArgList) # argument string '(x,y,...)' for RHS definition
        
    def __setField(self, FldDict): # field definitions given in the form name:definition, w/o signature
        inFldDict = {key:(self.fArgList, value) for key, value in FldDict.items()}
        self.fndict.update(inFldDict)
        
    def getField(self, fulldef = False):
        if fulldef: return self.fndict
        else: return {key:value[1] for key,value in self.fndict.items()}

    def frontKick(self, state, particle):
        return list(state)
    
    def rearKick(self, state, particle):
        return list(state)

class Drift(Element):
    """ drift space
    """
    fCount = 0 # counts the number of existing drifts
    
    def __init__(self, Length, Name = "Drift"):
        Element.__init__(self, 0, Length, Name+"_"+str(Drift.fCount))
        Drift.fCount += 1

class MQuad(Element):
    """ magnetic quadrupole
    """
    
    fCount = 0 # counts the # of existing mquads

    def __init__(self, Length, Grad, Name = "MQuad"):
        Element.__init__(self, 0, Length, Name+"_"+str(MQuad.fCount))
        self.setGrad(Grad)   
        MQuad.fCount += 1
        
    def setGrad(self, value):
        self._Element__setField({'Bx':str(value)+'*(-y)','By':str(value)+'*(-x)'})
        self.__fGrad = value
        
    def getGrad(self):
        return self.__fGrad
        

class MDipole(Element):
    """ bending magnetic dipole (horizontally bending);
    define _BField as a tuple;
    could also be used as a solenoid, if _BField = (0,0,Bs)
    and fCurve = 0
    """
    
    fCount = 0
    
    def __init__(self, Length, R, BField, Name = "MDipole"):
        Element.__init__(self, 1/R, Length, Name+"_"+str(MDipole.fCount))
        self.setBField(BField)
        MDipole.fCount += 1
        
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
  

class Solenoid(MDipole):
    
    fCount = 0
    
    def __init__(self, Length, Bs, Name = "Solenoid"):
        Element.__init__(self, 0, Length, Name+"_"+str(Solenoid.fCount))
        self.setBField((0,0,Bs))
        Solenoid.fCount += 1
        
    @classmethod
    def computeBStrength(cls, *args):
        print('Not defined')
        
    @classmethod
    def computeRadius(cls, *args):
        print('Infinite radius = zero curvature')
        

class MSext(Element):
    """ magnetic sextupole
    """
    
    fCount = 0
    
    def __init__(self, Length, Grad, Name = "MSext"):
        Element.__init__(self, 0, Length, Name+"_"+str(MSext.fCount))
        self.setGrad(Grad)
        MSext.fCount += 1
        
    def setGrad(self, value):
        self._Element__setField({'Bx':str(value)+'*(x*y)','By':str(value)+'*(x*x - y*y)','Bs':'0'})
        self.__fGrad = value
        
    def getGrad(self):
        return self.__fGrad
        
class Wien(Element):
    """ wien filter ?!?! 
    The magnetic field definition is weird
    """
    
    fCount = 0
    
    def __init__(self, Length, R, Hgap, RefPart, BField=None, EField=None, Name = "Wien?"):
        Element.__init__(self, 1/R, Length, Name+"_"+str(Wien.fCount))
        self.pardict.update(RefPart.pardict)
        
        self.__fR = [R, R - Hgap, R**2/(R-Hgap)]
        
        self.setEField(EField)
        self.setBField(BField)
        
        Wien.fCount += 1
        
    def __GammaBeta(self):
        Mass0 = self.pardict['Mass0']
        K0 = self.pardict['KinEn0']
        gamma = K0 / Mass0 + 1
        beta = NP.sqrt(gamma**2-1)/gamma
        return (gamma, beta)
    
    def __Pc(self, KNRG):
        return NP.sqrt((self.pardict['Mass0'] + KNRG)**2 - self.pardict['Mass0']**2)
    
    def setEField(self, EField=None):
        R = self.__fR
        if EField is None:
            P0c = self.__Pc(self.pardict['KinEn0'])
            gamma, beta = self.__GammaBeta()
            EField = - P0c*beta/R[0] * 1e6
        
        self.__fEField = (EField,0,0)
        self.__fVolt = (EField * R[0] * NP.log(R[2] / R[1])) / (-2)
        self._Element__setField({'Ex':str(EField)})
        
    def setBField(self, BField=None):        
        e0 = self.pardict['q']
        m0 = self.pardict['m0']
        clight = self.pardict['clight']
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
        u = -V + 2*V*NP.log((R+x)/R1)/NP.log(R2/R1)
        Xk[5] -= u*1e-6/particle.fKinEn0
        return Xk
        
    def rearKick(self, state, particle):
        x=state[0]
        Xk = list(state)
        R = self.__fR[0]
        R1 = self.__fR[1]
        R2 = self.__fR[2]
        V = self.__fVolt
        u = -V + 2*V*NP.log((R+x)/R1)/NP.log(R2/R1)
        Xk[5] += u*1e-6/particle.fKinEn0
        return Xk
