import numpy as NP
import pandas as PDS
import collections as CLN
import CParticle as PCL
import utilFunc as U
import copy
import re
import PyDSTool as DST


class Element:
    
    fArgList = ['x','y','ts','px','py','dK','H','s','Sx','Sy','Ss']
    
    fArgStr = None
    
    def __init__(self, Curve, Length, Name = "Element"):
        self.fName = Name
        
        self.fGeomdict = {'Curve':Curve, 'Length':Length}
        
        self.fPardict = {key:value for key,value in self.fGeomdict.items()}
        
        self.fFndict = { # defines the element EM field
                'Ex':(self.fArgList, '0'),
                'Ey':(self.fArgList, '0'),
                'Es':(self.fArgList, '0'),
                'Bx':(self.fArgList, '0'),
                'By':(self.fArgList, '0'),
                'Bs':(self.fArgList, '0')
                }
        
        self.fArgStr = U.phi(',',*self.fArgList) # argument string '(x,y,...)' for RHS definition
        super().__init__()
        
    def __repr__(self):
        return str(self.getField())
    
    def setField(self, **kwargs):
        self.__setField(dict(**kwargs))
        
    def __setField(self, FldDict): # field definitions given in the form name:definition, w/o signature
        inFldDict = {key:(self.fArgList, value) for key, value in FldDict.items()}
        self.fFndict.update(inFldDict)
        
    def getField(self, *args):
        if len(args)==0 : rval = {key:value[1] for key,value in self.fFndict.items()}    
        elif len(args) == 1 : rval = self.fFndict[args[0]][1]        
        else:rval = {key:self.fFndict[key][1] for key in args}
        
        return rval

    def getGeometry(self):
        return self.fGeomdict
    
    def frontKick(self):
        dK = DST.Var('dK')
        print('front kick, {}'.format(self.fName))
        return DST.Fun(dK, self.fArgList,'Front')
    
    def rearKick(self):   
        dK = DST.Var('dK')
        print('rear kick, {}'.format(self.fName))
        return DST.Fun(dK, self.fArgList,'Rear')
    
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
        super().__init__(Curve=0,Length=Length,Name=Name)
        

class MQuad(Element, HasCounter):
    """ magnetic quadrupole
    """

    def __init__(self, Length, Grad, Name = "MQuad"):
        super().__init__(Curve=0,Length=Length,Name=Name)
        self.setGrad(Grad)   
        
    def setGrad(self, value):
        self._Element__setField({'Bx':str(value)+'*(y)','By':str(value)+'*(x)'})
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
        super().__init__(Curve=1/R,Length=Length,Name=Name)
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
  

class MSol(Element, HasCounter):
    
    def __init__(self, Length, Bs, Name = "Solenoid"):
        super().__init__(Curve=0,Length=Length,Name=Name)
        self.setBField(Bs)
        
    def setBField(self, BField):
        self.__fBField = (0,0,BField)
        self._Element__setField({'Bs':BField})
        
    def getBField(self):
        return self.__fBField
        

class MSext(Element, HasCounter):
    """ magnetic sextupole
    """
    
    def __init__(self, Length, Grad, Name = "MSext"):
        super().__init__(Curve=0,Length=Length,Name=Name)
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
    
    __f0 = None
    
    def __init__(self, Length, Hgap, RefPart, EField=None, BField=None, R = None, Name = "Wien?"):
        
        if all([e is None for e in [EField, R]]): raise Exception("Either the E-field, or Radius must be defined")
        
        if R is None: Crv = None
        else: 
            Crv = 1/R
            self.__fR = [R, R - Hgap, R**2/(R-Hgap)]
            
        self.__fHgap = Hgap
        
        super().__init__(Curve=Crv,Length=Length,Name=Name)
        self.fGeomdict.update({'Hgap':Hgap})
        self.fPardict.update(RefPart.fPardict)
        
        self.setEField(EField)
        self.setBField(BField)
        
    def __GammaBeta(self):
        Mass0 = self.fPardict['Mass0']
        K0 = self.fPardict['KinEn0']
        gamma = K0 / Mass0 + 1
        beta = float(NP.sqrt(gamma**2-1)/gamma)
        return (gamma, beta)
    
    def __Pc(self, KNRG):
        return float(NP.sqrt((self.fPardict['Mass0'] + KNRG)**2 - self.fPardict['Mass0']**2))
    
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
        
        self.__fEField = (EField,0,0)
        self.__fVolt = (EField * R[0] * NP.log(R[2] / R[1])) / (-2)
        self.fPardict.update({'Curve':1/R[0]})
        self.fGeomdict.update({'R':R[0], 'Curve':1/R[0]})
        crv = str(self.fPardict['Curve'])
        self._Element__setField({'Ex':str(EField)+'/(1+'+crv+'*x)'})
        
        #define the voltage function for use in kicks
        R0 = float(R[0])
        R1 = float(self.__fR[1])
        R2 = float(self.__fR[2])
        V = float(self.__fVolt)
        x = DST.Var('x')
        
        self.__f0 = DST.Fun(x/R0 - .5*(x/R0)**2 + 1/3*(x/R0)**3 - .25*(x/R0)**4 +
                     .2*(x/R0)**5 - 1/6*(x/R0)**6 + 1/7*(x/R0)**7 -
                     .125*(x/R0)**8 + 1/9*(x/R0)**9 - .1*(x/R0)**10,['x'],'taylor') #log(1+x/R)
        
        self.__U = DST.Fun(-V + 2*V*(DST.Log(R0/R1) + self.__f0(x))/DST.Log(R2/R1), ['x'], 'volts') # DK = q*U
        
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
        
        B = U.phi('*', 0.018935,-BField,U.phi('+',-h, U.phi('*',k,BField,v)), 'x')
        B = B[1:len(B)-1]
        
        self._Element__setField({'By':B})
        
    def frontKick(self):
        KinEn0 = float(self.fPardict['KinEn0'])
        
        x = DST.Var('x')
        dK = DST.Var('dK')
        
        f = DST.Fun(dK - self.__U(x)*1e-6/KinEn0,self.fArgList,'Front')
        print('front kick, {}'.format(self.fName))
        return f
        
    def rearKick(self):

        KinEn0 = float(self.fPardict['KinEn0'])
            
        x = DST.Var('x')
        dK = DST.Var('dK')
        
        f = DST.Fun(dK + self.__U(x)*1e-6/KinEn0,self.fArgList,'Rear')
        print('rear kick, {}'.format(self.fName))
        return f
    
    def kickVolts(self, x):
        return (self.__fVolt, self.__U(x).tonumeric())
    
    
    
class ERF(Element, HasCounter):
    """ RF element
    """
    
    def __init__(self, Length, RefPart, Acc_length, EField = 1500, Phase = NP.pi/2, H_number = 50, Name = "RF"):
        super().__init__(Curve=0,Length=Length,Name=Name)
        if type(RefPart) is PCL.Ensemble: RefPart = RefPart.getReference()
        elif type(RefPart) is PCL.Particle: pass
        else: raise ValueError('Wrong type Reference Particle')
        
        self.fAmplitude = EField
        self.fPhase = Phase
        self.fFreq = RefPart.revFreq(Acc_length) * H_number
        self.__fH_number = H_number
        
        self.__fChars = PDS.DataFrame({'Amplitude':self.fAmplitude, 
                       'Frequency':self.fFreq,'h-number': self.__fH_number, 
                       'Phase':self.fPhase},index=[self.fName]).T
            
        self.__fU = self.fAmplitude*self.fLength
        
        RefPart.fRF = {'Amplitude':self.fAmplitude,'Freq':self.fFreq, 'Phase':self.fPhase}
        
    def __repr__(self):
        return str(self.__fChars)
    
    def __setEField(self):
        A = str(self.fAmplitude)
        w = str(self.fFreq*2*NP.pi)
        phi = str(self.fPhase)
        Es = A+'*cos('+w+'*ts+'+phi+')'
        self._Element__setField({'Es': Es})
    
    def EField_vec(self, time_vec):
        time_vec = NP.array(time_vec)
        z = NP.zeros(len(time_vec))
        A = self.fAmplitude
        w = self.fFreq*2*NP.pi
        phi = self.fPhase
        return list(zip(z,z, A*NP.cos(w*time_vec+phi)))
    
    def frontKick(self, particle):
        dK = DST.Var('dK')
        KinEn0 = float(self.fPardict['KinEn0'])
        
        f = DST.Fun(dK - self.__fU*1e-6/KinEn0,self.fArgList,'Rear')
        print('rear kick, {}'.format(self.fName))
        return f
        
    def rearKick(self, particle):
        dK = DST.Var('dK')
        KinEn0 = float(self.fPardict['KinEn0'])
        
        f = DST.Fun(dK + self.__fU*1e-6/KinEn0,self.fArgList,'Rear')
        print('rear kick, {}'.format(self.fName))
        return f
        
        