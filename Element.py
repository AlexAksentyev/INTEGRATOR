import Particle as PCL
import pandas as PDS
import numpy as NP
import collections as CLN
import copy
import re
import math

from RHS import imap, index

"""
1. scipy.linalg.block_diag for constructing diagonal block matrices:
    m0,m1,m2,... --- n_tilt, 3x3 tilt matrices
    f --- vectorized force (3x n_states matrix)
    M = scipy.linalg.block_diag(m0,m1,..)
    F = scipy.linalg.block_diag(f,f,f,...)
    MF = M.dot(F) --- [[m0*F, 0, 0,...],[0, m1*F, 0, ..], ...]
2. for picking out diagonal blocks (f'_i = m_i*F):
    x = [MF[i*3:(i+1)*3,i*n_state:(i+1)*n_state] for i in range(int(MF.shape[0]/3))]
3. for making an RHS compatible tilted force vector:
    x = NP.reshape(x,n_tilt*3*n_state) # for flat array
    x = x.reshape((3,-1)) # for array [[fx00,fx01,..., fxij],[fyij],[fsij]] #i -- tilt matrix, j -- state vector index
"""

#%%

class Field(NP.ndarray):

    def __new__(cls, array, Element, dtype=float):
        obj = NP.asarray(array, dtype=dtype, order='C').view(cls)  #force C-representation by default 
        obj.Host = Element
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.Host = getattr(obj, 'Host', None)
        
    def vectorize(self,arg):
        n = len(index(arg,'x')[0])
        return NP.repeat(self, n).reshape(3,n) # Field type is preserved
    
    def tilt(self):
        return Field(self.Host.Tilt.Matrix.dot(self).A, self.Host)
    
#    def updateHost(self, new_host):
#        """ In case i have to implement element::deepcopy, to use in there
#        """
#        self.Host = new_host

#%%
class Tilt:
    def __init__(self):
        self.Angle = {}
        self.Matrix = NP.matrix(NP.identity(3))
        self.__set_rep_str()
        
        # rotation matrices for tilting
        Rx = lambda c,s: NP.matrix([[1, 0,  0], 
                       [0, c, -s],
                       [0, s,  c]])
        Ry = lambda c,s: NP.matrix([[ c, 0, s],
                        [ 0,  1, 0],
                        [-s, 0, c]])
        Rs = lambda c,s: NP.matrix([[c, -s, 0],
                        [s,  c, 0],
                        [0,   0,  1]])
        
        self.Rotation = {'X':Rx,'Y':Ry,'S':Rs}
        
    def tilt(self, order, *args, append=False):
        i0 = 0
        if append:
            keys = list(self.Angle.keys())
            try:
                i0 = keys[len(keys)-1][1] # index of the last made rotation
                i0 += 1 # start writing from next one
            except IndexError:
                pass
            
        order = order.upper()
        i0 = NP.repeat(i0, len(args))
        i = list(range(len(args)))
        keys = list(zip(order, i0 + i))
        ang = dict(zip(keys, zip(args, NP.radians(args))))
        
        if append: self.Angle.update(ang)
        else: self.Angle = ang.copy()
        
        c = {key:NP.cos(x[1]) for key,x in ang.items()}
        s = {key:NP.sin(x[1]) for key,x in ang.items()}
        
        if append: res = self.Matrix
        else: res = NP.matrix(NP.identity(3))
        for key in ang.keys():
            res = self.Rotation[key[0]](c[key],s[key]).dot(res)
         
        self.Matrix = res
        
        self.__set_rep_str()
        
    def __set_rep_str(self):
        axis = [x[0] for x in self.Angle.keys()]
        adeg = [x[0] for x in self.Angle.values()]
        arad = [x[1] for x in self.Angle.values()]
        pd = PDS.DataFrame({'Axis':axis,'Deg':adeg,'Rad':arad}).T
        self.__RepStr = str(pd)
        
    def __repr__(self):
        return self.__RepStr
#%%

class Element:
    
    def __init__(self, Curve, Length, Name = "Element",**kwargs):
        self.fCurve = Curve
        self.fLength = Length
        self.fName = Name
        
        self.__fEField = (0,0,0)
        self.__fBField = (0,0,0)
        
        self.Tilt = Tilt()
        
        self.bSkip = False # for testing ERF.advance()
        
        self.__fChars = PDS.DataFrame({'Curve':self.fCurve, 'Length':self.fLength}, 
                                            index=[self.fName])
        
        super().__init__(**kwargs)
    
    def getFields(self):
        return self.__fEField, self.__fBField
    
    def EField(self,arg):   
        return Field(self.__fEField,self).vectorize(arg).tilt()
    
    def Eprime_tp(self, arg): # added for testing with ERF
        return Field((0,0,0),self).vectorize(arg).tilt()
    
    def BField(self,arg):
        return Field(self.__fBField,self).vectorize(arg).tilt()
                
    def frontKick(self,particle):
        pass # do nothing
    
    def rearKick(self,particle):
        pass # do nothing
        
    def printFields(self):
        print('Ex {}, Ey {}, Es {}'.format(*self.__fEField))
        print('Bx {}, By {}, Bs {}'.format(*self.__fBField))
    
    def __repr__(self):
        return str(self._Element__fChars)
    
    def tilt(self, order, *args, append=False):
        self.Tilt.tilt(order, *args, append=append)
        
class Bend:
    def __init__(self,RefPart,**kwargs):
        q = PCL.ezero
        clight = PCL.clight
        self.fPardict = {'KinEn0':RefPart.KinEn0, 'Mass0':RefPart.Mass0,
                         'q':q, 'clight':clight,
                         'm0':RefPart.Mass0/clight**2*1e6*q}
        
        super().__init__(**kwargs)
        
    def __GammaBeta(self):
        Mass0 = self.fPardict['Mass0']
        K0 = self.fPardict['KinEn0']
        gamma = K0 / Mass0 + 1
        beta = math.sqrt(gamma**2-1)/gamma
        return (gamma, beta)
    
    def __Pc(self, KNRG):
        return math.sqrt((self.fPardict['Mass0'] + KNRG)**2 - self.fPardict['Mass0']**2)
        
#%%
class Drift(Element):
    """ drift space
    """
    
    def __init__(self, Length, Name = "Drift"):
        super().__init__(Curve=0, Length=Length, Name=Name)
        

class MQuad(Element):
    """ magnetic quadrupole
    """
    
    def __init__(self, Length, Grad, Name = "MQuad"):
        super().__init__(Curve=0, Length=Length, Name=Name)
        self.__fGrad = Grad
        self._Element__fChars['Grad'] = self.__fGrad
        
    def BField(self, arg):
        i_x, i_y = index(arg, 'x','y')
        x = arg[i_x]
        y = arg[i_y]
        fld = (self.__fGrad*y, self.__fGrad*x, NP.zeros(len(x)))
        return  Field(fld, self).tilt()
        

class MDipole(Element, Bend):
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
             R = self.computeRadius(RefPart.KinEn0, BField)
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
  

class Solenoid(Element):
    
    def __init__(self, Length, Bs, Name = "Solenoid"):
        super().__init__(Curve=0, Length=Length, Name=Name)
        self._Element__fBField = (0,0,Bs)
        

class MSext(Element):
    """ magnetic sextupole
    """
    
    def __init__(self, Length, Grad, Name = "MSext"):
        super().__init__(Curve=0, Length=Length, Name=Name)
        self.__fGrad = Grad
        self._Element__fChars['Grad'] = self.__fGrad
        
    def BField(self, arg):
        i_x, i_y = index(arg, 'x','y')
        x = arg[i_x]
        y = arg[i_y]
        fld = (self.__fGrad*x*y,.5*self.__fGrad*(x**2 - y**2), NP.zeros(len(x)))
        return Field(fld, self).tilt()
        
class Wien(Element, Bend):
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
        i_x, = index(arg,'x')
        x = arg[i_x]
        Ex = self._Element__fEField[0]/(1+self.fCurve*x)
        z = NP.zeros(len(x))
        fld = (Ex, z, z)
        return Field(fld, self).tilt()
    
    def BField(self, arg):
        i_x, = index(arg,'x')
        x = arg[i_x]
        
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
        z = NP.zeros(len(x))
        fld = (z, B1, z)
        return Field(fld, self).tilt()
    
    def frontKick(self, state):
        u = self.__U(state[:,imap['x']])
        state[:,imap['dK']] -= u*1e-6/self.fPardict['KinEn0']
        
    def rearKick(self, state):
        u = self.__U(state[:,imap['x']])
        state[:,imap['dK']] += u*1e-6/self.fPardict['KinEn0']
        
    def kickVolts(self, x):
        return (self.__fVolt, self.__U(x))
        
    
class ERF(Element):
    """ RF element
    """
    
    def __init__(self, Length, Ensemble, Acc_length, EField = 15e5, Phase = 1.5*NP.pi, H_number = 50, Name = "RF"):
        super().__init__(Curve=0,Length=Length,Name=Name)
        
        if Length==0: 
            self.bSkip = True
            Length = 5e-4
        
        self.RefPart = Ensemble.Particle
        
        self.fAmplitude = EField
        self.fPhase = Phase
        self.fFreq = self.RefPart.revFreq(Acc_length) * H_number
        self.__fH_number = H_number
        
        self._Element__fChars = PDS.DataFrame({'Amplitude':self.fAmplitude, 
                       'Frequency':self.fFreq,'h-number': self.__fH_number, 
                       'Phase':self.fPhase},index=[self.fName]).T
            
        self.__fU = self.fAmplitude*Length # Length instead self.fLength for compatibility with Length 0
        
    def EField(self, arg):
        i_t, = index(arg, 't')
        t = arg[i_t]
        A = self.fAmplitude
        w = self.fFreq*2*NP.pi
        phi = self.fPhase
        z = NP.zeros(len(t))
        fld = (z, z, A*NP.cos(w*t+phi))
        return Field(fld, self).tilt()
    
    def Eprime_tp(self, arg): # Es prime divided by time prime
        i_t, = index(arg, 't')
        t = arg[i_t]
        A = self.fAmplitude
        w = self.fFreq*2*NP.pi
        phi = self.fPhase
        z = NP.zeros(len(t))
        fld = (z, z, -A*w*NP.sin(w*t+phi))
        return Field(fld, self).tilt()
    
    def advance(self, state):  
        """ alternative to integration,
            supposed to mimick a discrete kick
        """
        w = self.fFreq*2*NP.pi
        u = self.__fU
        i_dK, i_s, i_t = index(state, 'dK','s','t')
        K = self.RefPart.KinEn0 * (1 + state[i_dK])
        
        state[i_dK] += u*NP.cos(w*state[i_t]+self.fPhase)*1e-6/self.RefPart.KinEn0
        state[i_s] += self.fLength
        gamma,beta = self.RefPart.GammaBeta(K)
        state[i_t] += self.fLength/beta/PCL.clight
    
    def frontKick(self, state):
        u = self.__fU
        state[:,imap['dK']] -= u*1e-6/self.RefPart.KinEn0
        
    def rearKick(self, state):
        u = self.__fU
        state[:,imap['dK']] += u*1e-6/self.RefPart.KinEn0
        
    def kickVolts(self):
        return self.__fU
        
class Lattice:
    def __init__(self, ElSeq, Name):
        
        super().__init__()
        
        self.Sequence = list(ElSeq)
        self.ElCount = len(ElSeq)
        self.Name = Name
        
        self.RFCount = 0
        self.Length = 0
        for e in self:
            self.Length += e.fLength
            if isinstance(e, ERF): self.RFCount += 1
            
        if self.RFCount > 1: 
            print('More than one RF element in element sequence')
            return
        
    def copy(self):
        pass
        
    def __getitem__(self, idx):
        return self.Sequence[idx]
        
    def __repr__(self):
        return self.Sequence.__repr__()
    
    def __iter__(self):
        self.__current_id=0
        return self
    
    def __next__(self):
        last_id = self.ElCount - 1
        if self.__current_id <= last_id:
            result = self[self.__current_id]
            self.__current_id += 1
            return result
        else:
            raise StopIteration
            
    def insertRF(self, position, length, Ensemble, **ERF_pars):
        if isinstance(self[position], ERF):
            print('Replacing RF {}'.format(self.pop(position)))
            self.RFCount -= 1
        full_acc_len = self.Length + length
        rf = ERF(length,Ensemble, full_acc_len, **ERF_pars)
        self.IndRF = position
        self.RFCount += 1
        self.ElCount += 1
        self.Sequence.insert(position, rf)
        
    def getRF(self):
        try:
            return self[self.IndRF]
        except AttributeError:
            return None
        
    def pop(self, index):
        self.ElCount -= 1
        return self.Sequence.pop(index)
        
    def listNames(self, full=False):
        names = [e.fName for e in self]
        if full:
            return names
        else:
            return NP.unique([re.sub('_.*','',e) for e in names])        
        
    def tilt(self, order='S', mean_angle=(0,), sigma=(0,), append=False):
        n = len(order)
        if not isinstance(mean_angle, CLN.Sequence): mean_angle = (mean_angle,)
        if not isinstance(sigma, CLN.Sequence): sigma = (sigma,)
        nmean = len(mean_angle)
        nsig = len(sigma)           
        
        try:
            angle = NP.random.normal(mean_angle, sigma, size=(self.ElCount, n))
        except ValueError:
            print('Dimension mismatch: order {}, mean {}, sigma {}'.format(n, nmean, nsig))
            return
        
        nuids = len(NP.unique([id(e.Tilt) for e in self]))
        if nuids != self.ElCount:
            print('Non-unique elements ({}/{}) in lattice. Smart tilting.'.format(self.ElCount-nuids,self.ElCount))
            print('\t Not tilting:')
        
        i=0
        ids = set()
        cnt = 1
        for element in self:
            eid = id(element.Tilt)
            if eid in ids:
                print('\t\t {} element {} at lattice index {}'.format(cnt, element.fName, i))
                cnt += 1
                i += 1
                continue
            element.tilt(order,*angle[i], append=append)
            ids.add(eid)
            i +=1
    
#%%
if __name__ is '__main__':
    import Ensemble as ENS
    import Element as ENT
    from BNL import SSb1H2, BDA, QFS, QDA2
    
    E = ENS.Ensemble.populate(PCL.Particle(), Sz=1, x=(-1e-3,1e-3,5),dK=(0,3e-4,3))    
    tE = copy.deepcopy(E)
    state = NP.array(ENS.StateList(Sz=1, x=(-1e-3,1e-3,5),dK=(0,3e-4,3)).as_list()).flatten()
    
    el = BDA()    
    el.BField(state)
    #%%
    if False:
        mqf = MQuad(5e-2,86,'QF')
        mqd = MQuad(5e-2,-83,'QD')
        dft = Drift(25e-2)
        dft2 = Drift(25e-2)
        FODO = [mqf, dft, mqd, dft]
        
        lat = ENT.Lattice([BDA(),dft,dft,dft],'QFS')
        lat.insertRF(0,0,E,EField=15e7)
        
        tlat = copy.deepcopy(lat)
        
        tlat.tilt('ys',10,.003)
        
        #%%
        if True:
            E.track(lat, 500)
            tE.track(tlat, 30)
        
        #%%
            from matplotlib import pyplot as PLT
            PLT.figure()
            PLT.subplot(2,1,1)
            E.plot('-D x','s', pids=[3,6],new_plot=False)
            PLT.subplot(2,1,2)
            tE.plot('-D x','s', pids=[3,6],new_plot=False)
        

