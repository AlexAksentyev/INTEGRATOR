import CParticle as PCL
import pandas as PDS
import numpy as NP
import collections as CLN
import copy
import re
import math

class Element:
    
    def __init__(self, Curve, Length, Name = "Element",**kwargs):
        self.fCurve = Curve
        self.fLength = Length
        self.fName = Name
        
        self.__fEField = (0,0,0)
        self.__fBField = (0,0,0)
        
        self.bSkip = False # for testing ERF.advance()
        
        self.__fChars = PDS.DataFrame({'Curve':self.fCurve, 'Length':self.fLength}, 
                                            index=[self.fName])
        
        super().__init__(**kwargs)
    
    def EField(self,arg):
        return self._Element__fEField
    
    def Eprime_tp(self, arg): # added for testing with ERF
        return self._Element__fEField
    
    def BField(self,arg):
        return self._Element__fBField

    def frontKick(self,particle):
        pass # do nothing
    
    def rearKick(self,particle):
        pass # do nothing
        
    def printFields(self):
        print('Ex {}, Ey {}, Es {}'.format(*self._Element__fEField))
        print('Bx {}, By {}, Bs {}'.format(*self._Element__fBField))
    
    def __repr__(self):
        return str(self._Element__fChars)
        
class Bend:
    def __init__(self,RefPart,**kwargs):
        q = PCL.ezero
        clight = PCL.clight
        self.fPardict = {'KinEn0':RefPart.fKinEn0, 'Mass0':RefPart.fMass0,
                         'q':q, 'clight':clight,
                         'm0':RefPart.fMass0/clight**2*1e6*q}
        
        super().__init__(**kwargs)
        
    def __GammaBeta(self):
        Mass0 = self.fPardict['Mass0']
        K0 = self.fPardict['KinEn0']
        gamma = K0 / Mass0 + 1
        beta = math.sqrt(gamma**2-1)/gamma
        return (gamma, beta)
    
    def __Pc(self, KNRG):
        return math.sqrt((self.fPardict['Mass0'] + KNRG)**2 - self.fPardict['Mass0']**2)
        
class HasCounter:
    fCount = 0
    
    __fSep = "_"
    
    def __init__(self, **kwargs):
        if 'fName' not in self.__dict__.keys(): self.fName = 'NoName' 
#        self.fName = Name   
        self.fName += self.__fSep+str(self.__class__.fCount)
        self.__class__.fCount += 1
        super().__init__(**kwargs)

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
        self._Element__fChars['Grad'] = self.__fGrad
        
    def BField(self, arg):
        x = arg['x']
        y = arg['y']
        return (self.__fGrad*y, self.__fGrad*x,0)
        

class MDipole(Element, HasCounter, Bend):
    """ bending magnetic dipole (horizontally bending);
    define _BField as a tuple
    """
    
    def __init__(self, Length, RefPart, R=None, BField=None, Name = "MDipole"):
        if all([e is None for e in [BField, R]]): raise Exception("Either the B-field, or Radius must be defined")
        
        if R is None: # if R isn't given, BField is
            if isinstance(BField, CLN.Sequence): By = BField[1] # we want the By component
            else: By = BField
            if By == 0: 
                Crv = 0 # if By == 0, dipole is turned off == drift space; skip computeRadius
                self.__fR = NP.Inf
            else: Crv = None # otherwise, computeRadius
        else:
            Crv = 1/R
            self.__fR = R
        
        super().__init__(Curve=Crv, Length=Length, Name=Name, RefPart=RefPart)
        
        if Crv is None:
             R = self.computeRadius(RefPart.fKinEn0, BField)
             self.__fR = R
             self.fCurve = 1/self.__fR
             self._Element__fChars['Curve'] = self.fCurve
             
        self.setBField(BField)
        
    def setBField(self,BField):
        if BField is None:
            R = self.__fR
            BField = (0, self.computeBStrength(self.__fPardict['KinEn0'], R), 0)
        elif not isinstance(BField, CLN.Sequence): BField = (0,BField,0) # by default B = By
        
        self._Element__fBField = BField[:]
        self._Element__fChars['Bx'] = BField[0]
        self._Element__fChars['By'] = BField[1]
        self._Element__fChars['Bs'] = BField[2]
        
        
    def getBField(self):
        return self._Element__fBField[:]
    
    def computeBStrength(self,KinEn, R):
        return self._Bend__Pc(KinEn)*1e6/(R*PCL.clight)
    
    def computeRadius(self,KinEn, BField): 
        return self._Bend__Pc(KinEn)*1e6/(BField*PCL.clight)
  

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
        self._Element__fChars['Grad'] = self.__fGrad
        
    def BField(self, arg):
        x = arg['x']
        y = arg['y']
        return (self.__fGrad*x*y,.5*self.__fGrad*(x**2 - y**2), 0)
        
class Wien(Element, HasCounter, Bend):
    """ wien filter
    """
    
    def __init__(self, Length, Hgap, RefPart, EField=None, BField=None, R=None, Name = "Wien"):
        
        if all([e is None for e in [EField, R]]): raise Exception("Either the E-field, or Radius must be defined")
        
        if R is None: Crv = None
        else: 
            Crv = 1/R
            self.__fR = [R, R - Hgap, R**2/(R-Hgap)]
            
        self.__fHgap = Hgap
        
        super().__init__(Curve=Crv,Length=Length,Name=Name, RefPart=RefPart)
        
        self.setEField(EField)
        self.setBField(BField)
        
    def setEField(self, EField=None):
        P0c = self._Bend__Pc(self.fPardict['KinEn0'])
        gamma, beta = self._Bend__GammaBeta()
        if EField is None:
            R = self.__fR
            EField = - P0c*beta/R[0] * 1e6
        else:
            assert EField < 0, "Incorrect field value ({} >= 0)".format(EField)
            R = - P0c*beta/EField * 1e6
            self.__fR = [R, R - self.__fHgap, R**2/(R-self.__fHgap)]
            R = self.__fR
        
        self.fCurve = 1/R[0]
        self._Element__fEField = (EField,0,0)
        self.__fVolt = (EField * R[0] * math.log(R[2] / R[1])) / (-2)
        self._Element__fChars['Curve'] = self.fCurve
        
        #define a subfunction for use in kicks
        R0 = float(R[0])
        R1 = float(self.__fR[1])
        R2 = float(self.__fR[2])
        V = float(self.__fVolt)
        
        self.__U = lambda x: (-V + 2*V*NP.log((R0+x)/R1)/NP.log(R2/R1)) # DK = q*U
    
    def setBField(self, BField=None):        
        
        clight = self.fPardict['clight']
        gamma, beta = self._Bend__GammaBeta()        
        
        v = beta*clight
        
        if BField is None:
            if self._Element__fEField is None: self.setEField()
            BField = -self._Element__fEField[0]/v
        
        self._Element__fBField = (0, BField, 0)
        
    def EField(self, arg):
        x = arg['x']
        Ex = self._Element__fEField[0]/(1+self.fCurve*x)
        return (Ex, 0, 0)
    
    def BField(self, arg):
        x =  arg['x']
        
        e0 = self.fPardict['q']
        m0 = self.fPardict['m0']
        clight = self.fPardict['clight']
        qm = e0/(m0*clight**2)
        
        gamma, beta = self._Bend__GammaBeta()
        v = beta*clight
        
        k = 1.18 * qm * ((2 - 3*beta**2 - .75*beta**4)/beta**2 + 1/gamma)
        h = 1/self.__fR[0]
        
        B0 = self._Element__fBField[1]
        B1 = .018935*(-B0)*(-h+k*v*B0)*x
        return (0, B1, 0)
    
    def frontKick(self, particle):
        x=particle.getState()['x']
        u = self.__U(x)
        Xk = particle.getState()
        Xk['dK'] -= u*1e-6/particle.fKinEn0
#        print('Kick voltage {}'.format(u))
        particle.setState(Xk)
        
    def rearKick(self, particle):
        x=particle.getState()['x']
        u = self.__U(x)
        Xk = particle.getState()
        Xk['dK'] += u*1e-6/particle.fKinEn0
#        print('Kick voltage {}'.format(u))
        particle.setState(Xk)
        
    def kickVolts(self, x):
        return (self.__fVolt, self.__U(x))
        
    
class ERF(Element, HasCounter):
    """ RF element
    """
    
    def __init__(self, Length, RefPart, Acc_length, EField = 1500, Phase = 1.5*NP.pi, H_number = 50, Name = "RF"):
        super().__init__(Curve=0,Length=Length,Name=Name)
        if type(RefPart) is PCL.Ensemble: RefPart = RefPart.getReference()
        elif type(RefPart) is PCL.Particle: pass
        else: raise ValueError('Wrong type Reference Particle')
        
        self.fAmplitude = EField
        self.fPhase = Phase
        self.fFreq = RefPart.revFreq(Acc_length) * H_number
        self.__fH_number = H_number
        
        self._Element__fChars = PDS.DataFrame({'Amplitude':self.fAmplitude, 
                       'Frequency':self.fFreq,'h-number': self.__fH_number, 
                       'Phase':self.fPhase},index=[self.fName]).T
            
        self.__fU = self.fAmplitude*self.fLength
        
        RefPart.fRF = {'Amplitude':self.fAmplitude,'Freq':self.fFreq, 'Phase':self.fPhase}
        
    def EField(self, arg):
        t = arg['t']
        A = self.fAmplitude
        w = self.fFreq*2*NP.pi
        phi = self.fPhase
        return (0, 0, A*NP.cos(w*t+phi))
    
    def Eprime_tp(self, arg): # Es prime divided by time prime
        t = arg['t']
        A = self.fAmplitude
        w = self.fFreq*2*NP.pi
        phi = self.fPhase
        return (0, 0, -A*w*NP.sin(w*t+phi))
    
    def EField_vec(self, time_vec):
        time_vec = NP.array(time_vec)
        z = NP.zeros(len(time_vec))
        A = self.fAmplitude
        w = self.fFreq*2*NP.pi
        phi = self.fPhase
        return list(zip(z,z, A*NP.cos(w*time_vec+phi)))
    
    def advance(self, particle):
        arg = particle.getState()
        
        w = self.fFreq*2*NP.pi
        u = self.__fU
        K = particle.fKinEn0 * (1 + arg['dK'])
        
        arg['dK'] += u*NP.cos(w*arg['t']+self.fPhase)*1e-6/particle.fKinEn0
        arg['s'] += self.fLength
        gamma,beta = particle.GammaBeta(K)
        arg['t'] += self.fLength/beta/particle.CLIGHT()
        
        particle.setState(arg)
    
#    def frontKick(self, particle):
#        u = self.__fU
#        Xk = particle.getState()
#        Xk['dK'] -= u*1e-6/particle.fKinEn0
##        print('Kick voltage {}'.format(u))
#        particle.setState(Xk)
#        
#    def rearKick(self, particle):
#        u = self.__fU
#        Xk = particle.getState()
#        Xk['dK'] += u*1e-6/particle.fKinEn0
##        print('Kick voltage {}'.format(u))
#        particle.setState(Xk)
        
class Lattice:
    def __init__(self, ElSeq, RefPart):
        
        super().__init__()
        
        if type(RefPart) is PCL.Ensemble: self.fRefPart = RefPart.getReference()
        elif type(RefPart) is PCL.Particle: self.fRefPart = RefPart
        else: raise ValueError('Wrong type Reference Particle')
        
        self.fSequence = ElSeq[:]
        
        self.fCount = len(ElSeq)
        self.fLength = 0
        for e in ElSeq: self.fLength += e.fLength
        
    def insertRF(self, position, length, **ERF_pars):
        full_acc_len = self.fLength + length
        rf = ERF(length,self.fRefPart, full_acc_len, **ERF_pars)
        self.fSequence.insert(position, rf)
        
    def __repr__(self):
        pass