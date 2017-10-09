import numpy as NP
import collections as CLN



class Element:
    
    def __init__(self, Curve, Length, Name = "Element"):
        self.fCurve = Curve
        self.fLength = Length
        self.fName = Name
    
    def EField(self,arg):
        return (0,0,0)
    
    def BField(self,arg):
        return (0,0,0)

    def frontKick(self, state, particle):
        return list(state)
    
    def rearKick(self, state, particle):
        return list(state)

class Drift(Element):
    """ drift space
    """
    fCount = 0 # counts the number of existing drifts
    
    def __init__(self, Length, Name = "Driftspace"):
        Element.__init__(self, 0, Length, Name+str(Drift.fCount))
        Drift.fCount += 1

class MQuad(Element):
    """ magnetic quadrupole
    """
    
    fCount = 0 # counts the # of existing mquads
    
    def __init__(self, Length, Grad, Name = "MQuadrupole"):
        Element.__init__(self, 0, Length, Name+str(MQuad.fCount))
        self.__fGrad = Grad
        MQuad.fCount += 1
        
    def setGrad(self, value):
        self.__fGrad = value
        
    def getGrad(self):
        return self.__fGrad
        
    def BField(self, arg):
        x,y = arg[0:2]
        return (-self.__fGrad*x, -self.__fGrad*y,0)
        

class MDipole(Element):
    """ bending magnetic dipole (horizontally bending);
    define _BField as a tuple;
    could also be used as a solenoid, if _BField = (0,0,Bz)
    and fCurve = 0
    """
    
    fCount = 0
    
    def __init__(self, Length, R, BField, Name = "MDipole"):
        Element.__init__(self, 1/R, Length, Name+str(MDipole.fCount))
        self.setBField(BField)
        MDipole.fCount += 1
        
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
        
        
    def BField(self, arg):
        return self.__fBField
  

class Solenoid(MDipole):
    
    fCount = 0
    
    def __init__(self, Length, Bs, Name="Sol"):
        Element.__init__(self, 0, Length, Name+str(Solenoid.fCount))
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
    
    def __init__(self, Length, Grad, Name="MSext"):
        Element.__init__(self, 0, Length, Name+str(MSext.fCount))
        self.__fGrad = Grad
        MSext.fCount += 1
        
    def BField(self, arg):
        x,y=arg[0:2]
        return (self.__fGrad*x*y,.5*self.__fGrad*(x**2 - y**2), 0)
        
class Wien(Element):
    """ wien filter ?!?! 
    The magnetic field definition is weird
    """
    
    fCount = 0
    
    def __init__(self, Length, R, Hgap, Voltage, BField, Name = "WF"):
        Element.__init__(self, 1/R, Length, Name+str(Wien.fCount))
        self.__fR = [R, R - Hgap, R**2/(R-Hgap)]
        self.__fVolt = Voltage
        self.__fBField = BField
        Wien.fCount += 1
        
    @staticmethod
    def computeVoltage(particle, R, Hgap):
        gamma,beta = particle.GammaBeta(particle.fKinEn0)
        R = [R, R-Hgap, R**2/(R-Hgap)]
        E0 = - particle.Pc(particle.fKinEn0)*beta/R[0] * 1e6
        return (E0 * R[0] * NP.log(R[2] / R[1])) / (-2)
    
    @staticmethod
    def computeBStrength(particle, x, R, Hgap): # no idea about the end-formula here
        gamma,beta = particle.GammaBeta(particle.fKinEn0)
        R = [R, R-Hgap, R**2/(R-Hgap)]
        E0 = - particle.Pc(particle.fKinEn0)*beta/R[0] * 1e6
        qm = particle.EZERO()/particle.fMass0
        k = qm * ((2 - 3*beta**2 - .75*beta**4)/beta**2 + 1/gamma)
        
        v=beta*particle.CLIGHT()
        B0 = -E0/v # this yields .46 for the deuteron at 270 MeV, as expected 
        k = k*1.18
        return 0.018935*(-B0)*(-1/R[0] + k*B0*v)*x
        
    def EField(self, arg):
        x = arg[0]
        Ex = -2*self.__fVolt/(NP.log(self.__fR[2]/self.__fR[1])*(self.__fR[0]+x))
        return (Ex, 0, 0)
    
    def BField(self, arg):        
        return (0, self.__fBField, 0)
    
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