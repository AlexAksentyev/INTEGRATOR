import numpy as NP
import collections as CLN
import copy
import re
import math


class Element:
    
    def __init__(self, Curve, Length, Name = "Element"):
        self.fCurve = Curve
        self.fLength = Length
        self.fName = Name
        
        self.__fEField = (0,0,0)
        self.__fBFIeld = (0,0,0)
        
        super().__init__()
    
    def EField(self,arg):
        return self.__fEField
    
    def BField(self,arg):
        return self.__fBFIeld

    def frontKick(self,particle):
        pass # do nothing
    
    def rearKick(self,particle):
        pass # do nothing
        
        
class HasCounter:
    fCount = 0
    
    __fSep = "_"
    
    def __init__(self):
        if 'fName' not in self.__dict__.keys(): self.fName = 'NoName'    
        self.fName += self.__fSep+str(self.__class__.fCount)
        self.__class__.fCount += 1
        super().__init__()

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
        super().__init__(Curve=0, Length=Length, Name=Name)
        

class MQuad(Element, HasCounter):
    """ magnetic quadrupole
    """
    
    def __init__(self, Length, Grad, Name = "MQuad"):
        super().__init__(Curve=0, Length=Length, Name=Name)
        self.__fGrad = Grad
        
    def BField(self, arg):
        x = arg['x']
        y = arg['y']
        return (self.__fGrad*y, self.__fGrad*x,0)
        

class MDipole(Element, HasCounter):
    """ bending magnetic dipole (horizontally bending);
    define _BField as a tuple;
    could also be used as a solenoid, if _BField = (0,0,Bz)
    and fCurve = 0
    """
    
    def __init__(self, Length, R, BField, Name = "MDipole"):
        super().__init__(Curve=1/R, Length=Length, Name=Name)
        self.setBField(BField)
        
    def setBField(self,BField):
        if not isinstance(BField, CLN.Sequence): BField = (0,BField,0) # by default B = By
        self.__fBField = BField[:]
        
    def getBField(self):
        return self.__fBField[:]
    
    @classmethod    
    def computeBStrength(cls,particle, R):
        return particle.Pc(particle.fKinEn0)*1e6/(R*particle.CLIGHT())
    
    @classmethod
    def computeRadius(cls,particle, BField): 
        return particle.Pc(particle.fKinEn0)*1e6/(BField*particle.CLIGHT())
  

class Solenoid(Element, HasCounter):
    
    def __init__(self, Length, Bs, Name = "Solenoid"):
        super().__init__(self, Curve=0, Length=Length, Name=Name)
        self.__fField = (0,0,Bs)
        

class MSext(Element, HasCounter):
    """ magnetic sextupole
    """
    
    def __init__(self, Length, Grad, Name = "MSext"):
        super().__init__(Curve=0, Length=Length, Name=Name)
        self.__fGrad = Grad
        
    def BField(self, arg):
        x = arg['x']
        y = arg['y']
        return (self.__fGrad*x*y,.5*self.__fGrad*(x**2 - y**2), 0)
        
class Wien(Element, HasCounter):
    """ wien filter ?!?! 
    The magnetic field definition is weird
    """
    
    def __init__(self, Length, Hgap, RefPart, EField=None, BField=None, R=None, Name = "Wien?"):
        
        if all([e is None for e in [EField, R]]): raise Exception("Either the E-field, or Radius must be defined")
        
        if R is None: Crv = None
        else: 
            Crv = 1/R
            self.__fR = [R, R - Hgap, R**2/(R-Hgap)]
            
        self.__fHgap = Hgap
        
        q = RefPart.EZERO()
        clight = RefPart.CLIGHT()
        self.fPardict = {'KinEn0':RefPart.fKinEn0, 'Mass0':RefPart.fMass0,
                         'q':q, 'clight':clight,
                         'm0':RefPart.fMass0/clight**2*1e6*q}
        
        super().__init__(Curve=Crv,Length=Length,Name=Name)
        
        self.setEField(EField)
        self.setBField(BField)
        
    def setEField(self, EField=None):
        P0c = self.__Pc(self.fPardict['KinEn0'])
        gamma, beta = self.__GammaBeta()
        if EField is None:
            R = self.__fR
            EField = - P0c*beta/R[0] * 1e6
        else:
            assert EField < 0, "Incorrect field value ({} >= 0)".format(EField)
            R = - P0c*beta/EField * 1e6
            self.__fR = [R, R - self.__fHgap, R**2/(R-self.__fHgap)]
            R = self.__fR
        
        self.fCurve = 1/R[0]
        self.__fEField = (EField,0,0)
        self.__fVolt = (EField * R[0] * NP.log(R[2] / R[1])) / (-2)
        
        #define a subfunction for use in kicks
        R0 = float(R[0])
        R1 = float(self.__fR[1])
        R2 = float(self.__fR[2])
        V = float(self.__fVolt)
        
        self.__f0 = lambda x: (x/R0 - .5*(x/R0)**2 + 1/3*(x/R0)**3 - .25*(x/R0)**4 +
                     .2*(x/R0)**5 - 1/6*(x/R0)**6 + 1/7*(x/R0)**7 -
                     .125*(x/R0)**8 + 1/9*(x/R0)**9 - .1*(x/R0)**10) #log(1+x/R)
        
        self.__U = lambda x: (-V + 2*V*math.log((R0+x)/R1)/math.log(R2/R1)) # DK = q*U
    
    
    def kickVolts(self, x):
        return (self.__fVolt, self.__U(x))
    
    def setBField(self, BField=None):        
        
        clight = self.fPardict['clight']
        gamma, beta = self.__GammaBeta()        
        
        v = beta*clight
        
        if BField is None:
            if self.__fEField is None: self.setEField()
            BField = -self.__fEField[0]/v
        
        self.__fBField = (0, BField, 0)
        
    def EField(self, arg):
        x = arg['x']
        Ex = self.__fEField[0]/(1+self.fCurve*x)
        return (Ex, 0, 0)
    
    def BField(self, arg):
        x =  arg['x']
        
        e0 = self.fPardict['q']
        m0 = self.fPardict['m0']
        clight = self.fPardict['clight']
        qm = e0/(m0*clight**2)
        
        gamma, beta = self.__GammaBeta()
        v = beta*clight
        
        k = 1.18 * qm * ((2 - 3*beta**2 - .75*beta**4)/beta**2 + 1/gamma)
        h = 1/self.__fR[0]
        
        B0 = self.__fBField[1]
        B1 = .018935*(-B0)*(-h+k*v*B0)*x
        return (0, B1, 0)
    
    def __GammaBeta(self):
        Mass0 = self.fPardict['Mass0']
        K0 = self.fPardict['KinEn0']
        gamma = K0 / Mass0 + 1
        beta = float(NP.sqrt(gamma**2-1)/gamma)
        return (gamma, beta)
    
    def __Pc(self, KNRG):
        return float(NP.sqrt((self.fPardict['Mass0'] + KNRG)**2 - self.fPardict['Mass0']**2))
    
    def frontKick(self, particle):
        x=particle.getState()['x']
        u = self.__U(x)
        Xk = particle.getState()
        Xk['dK'] -= u*1e-6/particle.fKinEn0
        print('Kick voltage {}'.format(u))
        particle.setState(Xk)
        
    def rearKick(self, particle):
        x=particle.getState()['x']
        u = self.__U(x)
        Xk = particle.getState()
        Xk['dK'] += u*1e-6/particle.fKinEn0
        print('Kick voltage {}'.format(u))
        particle.setState(Xk)
        
    
class ERF(Element, HasCounter):
    """ RF element
    """
    
    def __init__(self, Length, RF_params = (0,0,0), Name = "RF"):
        super().__init__(Curve=0,Length=Length,Name=Name)
        self.fAmplitude, self.fFreq, self.fPhase = RF_params
        
    def EField(self, arg):
        t = arg['t']
        A = self.fAmplitude
        w = self.fFreq*2*NP.pi
        phi = self.fPhase
        return (0, 0, A*NP.cos(w*t+phi))