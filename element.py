
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
    x = x.reshape((3,-1)) # for array [[fx00,fx01,..., fxij],[fyij],[fsij]]
    # i -- tilt matrix, j -- state vector index
    
TODO:
    * think about using slots in elements, tilt, field
"""

import pandas as pds
import numpy as np
import collections as cln
import math

import particle as pcl

from rhs import IMAP, index

class Field(np.ndarray):
    """Representation of an element's EM-field;
    makes more sense to call field.tilt, instead of
    element.tilt(field). Same for vectorization.
    """
    def __new__(cls, array, element, dtype=float):
        obj = np.asarray(array, dtype=dtype, order='C').view(cls)
        obj.host = element
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.host = getattr(obj, 'host', None)

    def vectorize(self, arg):
        vec_num = len(index(arg, 'x')[0])
        return np.repeat(self, vec_num).reshape(3, vec_num) # Field type is preserved

    def tilt(self):
        return Field(self.host.tilt_.matrix.dot(self).A, self.host)

#    def updateHost(self, new_host):
#        """In case i have to implement element::deepcopy, to use in there."""
#        self.Host = new_host

#%%
class Tilt:
    def __init__(self):
        self.angle = {}
        self.matrix = np.matrix(np.identity(3))
        self.__set_rep_str()

        # rotation matrices for tilting
        Rx = lambda c, s: np.matrix([[1, 0, 0], [0, c, -s], [0, s, c]])
        Ry = lambda c, s: np.matrix([[c, 0, s], [0, 1, 0], [-s, 0, c]])
        Rs = lambda c, s: np.matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]])

        self.rotation = {'X': Rx, 'Y': Ry, 'S': Rs}

    def __call__(self, order, *degree_angles, append=False):
        i0 = 0
        if append:
            keys = list(self.angle.keys())
            try:
                i0 = keys[len(keys)-1][1] # index of the last made rotation
                i0 += 1 # start writing from next one
            except IndexError:
                pass

        order = order.upper()
        i0 = np.repeat(i0, len(degree_angles))
        i = list(range(len(degree_angles)))
        keys = list(zip(order, i0 + i))
        ang = dict(zip(keys, zip(degree_angles, np.deg2rad(degree_angles))))

        if append:
            self.angle.update(ang)
        else:
            self.angle = ang.copy()

        c = {key:np.cos(x[1]) for key, x in ang.items()}
        s = {key:np.sin(x[1]) for key, x in ang.items()}

        if append:
            res = self.matrix
        else:
            res = np.matrix(np.identity(3))
        for key in ang.keys():
            res = self.rotation[key[0]](c[key], s[key]).dot(res)

        self.matrix = res

        self.__set_rep_str()

    def __set_rep_str(self):
        axis = [x[0] for x in self.angle.keys()]
        adeg = [x[0] for x in self.angle.values()]
        arad = [x[1] for x in self.angle.values()]
        df = pds.DataFrame({'Axis': axis, 'Deg': adeg, 'Rad': arad}).T
        self.__rep_str = str(df)

    def __repr__(self):
        return self.__rep_str
#%%

class Element:

    def __init__(self, curve, length, name="Element", **kwargs):
        self.curve = curve
        self.length = length
        self.name = name

        self.__E_field = (0, 0, 0)
        self.__B_field = (0, 0, 0)

        self.tilt_ = Tilt()

        self.__bool_skip = False # for testing ERF.advance()

        self.__chars = pds.DataFrame({'Curve':self.curve, 'Length':self.length}, index=[self.name])

        super().__init__(**kwargs)
        
    @property
    def skip(self):
        return self.__bool_skip
    
#    @skip.setter
#    def skip(self, value):
#        print("Element doesn't have the advance method.")

    def get_fields(self):
        return self.__E_field, self.__B_field

    def EField(self, arg):
        return Field(self.__E_field, self).vectorize(arg).tilt()

    def EField_prime_t(self, arg): # added for testing with ERF
        return Field((0, 0, 0), self).vectorize(arg).tilt()

    def BField(self, arg):
        return Field(self.__B_field, self).vectorize(arg).tilt()

    def front_kick(self, state):
        """If I want not to reshape the state vector before writing to
        PLog, all I need to do is change state[:,IMAP['dK']] = ....
        to  i_dK = index('dK'); state[i_dK] = ...
        """
        pass # do nothing

    def rear_kick(self, state):
        pass # do nothing

    def print_fields(self):
        print('Ex {}, Ey {}, Es {}'.format(*self.__E_field))
        print('Bx {}, By {}, Bs {}'.format(*self.__B_field))

    def __repr__(self):
        return str(self.__chars)

    def tilt(self, order, *args, append=False):
        self.tilt_(order, *args, append=append)

class Bend:
    def __init__(self, reference_particle, **kwargs):
        q = pcl.EZERO
        clight = pcl.CLIGHT
        self.pardict = {'KinEn0':reference_particle.kin_nrg_0, 'Mass0':reference_particle.mass0,
                        'q':q, 'clight':clight,
                        'm0':reference_particle.mass0/clight**2*1e6*q}

        super().__init__(**kwargs)

    def __GammaBeta(self):
        mass0 = self.pardict['Mass0']
        K0 = self.pardict['KinEn0']
        gamma = K0 / mass0 + 1
        beta = math.sqrt(gamma**2-1)/gamma
        return (gamma, beta)

    def __Pc(self, KNRG):
        return math.sqrt((self.pardict['Mass0'] + KNRG)**2 - self.pardict['Mass0']**2)

#%%
class Drift(Element):
    """Drift space."""

    def __init__(self, length, name="Drift"):
        super().__init__(curve=0, length=length, name=name)


class MQuad(Element):
    """Magnetic quadrupole."""

    def __init__(self, length, grad, name="MQuad"):
        super().__init__(curve=0, length=length, name=name)
        self.__grad = grad
        self._Element__chars['Grad'] = self.__grad

    def BField(self, arg):
        i_x, i_y = index(arg, 'x', 'y')
        x = arg[i_x]
        y = arg[i_y]
        fld = (self.__grad*y, self.__grad*x, np.zeros(len(x)))
        return  Field(fld, self).tilt()


class MDipole(Element, Bend):
    """(Horizontally) Bending magnetic dipole."""

    def __init__(self, length, reference_particle, R=None, B_field=None, name="MDipole"):
        if all([e is None for e in [B_field, R]]):
            raise Exception("Either the B-field, or Radius must be defined")

        if R is None: # if R isn't given, B_field is
            if isinstance(B_field, cln.Sequence):
                By = B_field[1] # we want the By component
            else:
                By = B_field
            if By == 0:
                crv = 0 # if By == 0, dipole is turned off == drift space; skip computeRadius
                self.__radius = np.Inf
            else: crv = None # otherwise, computeRadius
        else:
            crv = 1/R
            self.__radius = R

        super().__init__(curve=crv, length=length, name=name, reference_particle=reference_particle)

        if crv is None:
            R = self.compute_radius(reference_particle.kin_nrg_0, B_field)
            self.__radius = R
            self.curve = 1/self.__radius
            self._Element__chars['Curve'] = self.curve

        self.set_B_field(B_field)

    def set_B_field(self, B_field):
        if B_field is None:
            R = self.__radius
            B_field = (0, self.compute_B_strength(self.pardict['KinEn0'], R), 0)
        elif not isinstance(B_field, cln.Sequence): B_field = (0, B_field, 0) # by default B = By

        self._Element__B_field = B_field[:]
        self._Element__chars['Bx'] = B_field[0]
        self._Element__chars['By'] = B_field[1]
        self._Element__chars['Bs'] = B_field[2]


    def get_B_field(self):
        return self._Element__B_field[:]

    def compute_B_strength(self, KinEn, R):
        return self._Bend__Pc(KinEn)*1e6/(R*pcl.CLIGHT)

    def compute_radius(self, KinEn, B_field):
        return self._Bend__Pc(KinEn)*1e6/(B_field*pcl.CLIGHT)

class Solenoid(Element):

    def __init__(self, length, Bs, name="Solenoid"):
        super().__init__(curve=0, length=length, name=name)
        self._Element__B_field = (0, 0, Bs)


class MSext(Element):
    """Magnetic sextupole."""

    def __init__(self, length, grad, name="MSext"):
        super().__init__(curve=0, length=length, name=name)
        self.__grad = grad
        self._Element__chars['Grad'] = self.__grad

    def BField(self, arg):
        i_x, i_y = index(arg, 'x', 'y')
        x = arg[i_x]
        y = arg[i_y]
        fld = (self.__grad*x*y, .5*self.__grad*(x**2 - y**2), np.zeros(len(x)))
        return Field(fld, self).tilt()

class Wien(Element, Bend):
    """Wien filter."""

    def __init__(self, length, h_gap, reference_particle, E_field=None, B_field=None, R=None, name="Wien"):

        if all([e is None for e in [E_field, R]]):
            raise Exception("Either the E-field, or Radius must be defined")

        if R is None: crv = None
        else:
            crv = 1/R
            self.__radius = [R, R - h_gap, R**2/(R-h_gap)]

        self.__h_gap = h_gap

        super().__init__(curve=crv, length=length, name=name, reference_particle=reference_particle)

        self.set_E_field(E_field)
        self.set_B_field(B_field)

    def set_E_field(self, E_field=None):
        P0c = self._Bend__Pc(self.pardict['KinEn0'])
        _, beta = self._Bend__GammaBeta()
        if E_field is None:
            R = self.__radius
            E_field = - P0c*beta/R[0] * 1e6
        else:
            assert E_field < 0, "Incorrect field value ({} >= 0)".format(E_field)
            R = - P0c*beta/E_field * 1e6
            self.__radius = [R, R - self.__h_gap, R**2/(R-self.__h_gap)]
            R = self.__radius

        self.curve = 1/R[0]
        self._Element__E_field = (E_field, 0, 0)
        self.__volt = (E_field * R[0] * math.log(R[2] / R[1])) / (-2)
        self._Element__chars['Curve'] = self.curve

        #define a subfunction for use in kicks
        R0 = float(R[0])
        R1 = float(self.__radius[1])
        R2 = float(self.__radius[2])
        V = float(self.__volt)

        self.__U = lambda x: (-V + 2*V*np.log((R0+x)/R1)/np.log(R2/R1)) # DK = q*U

    def set_B_field(self, B_field=None):

        clight = self.pardict['clight']
        _, beta = self._Bend__GammaBeta()

        v = beta*clight

        if B_field is None:
            if self._Element__E_field is None:
                self.set_E_field()
            B_field = -self._Element__E_field[0]/v

        self._Element__B_field = (0, B_field, 0)

    def EField(self, arg):
        i_x, = index(arg, 'x')
        x = arg[i_x]
        Ex = self._Element__E_field[0]/(1+self.curve*x)
        z = np.zeros(len(x))
        fld = (Ex, z, z)
        return Field(fld, self).tilt()

    def BField(self, arg):
        i_x, = index(arg, 'x')
        x = arg[i_x]

        e0 = self.pardict['q']
        m0 = self.pardict['m0']
        clight = self.pardict['clight']
        qm = e0/(m0*clight**2)

        gamma, beta = self._Bend__GammaBeta()
        v = beta*clight

        k = 1.18 * qm * ((2 - 3*beta**2 - .75*beta**4)/beta**2 + 1/gamma)
        h = 1/self.__radius[0]

        B0 = self._Element__B_field[1]
        B1 = .018935*(-B0)*(-h+k*v*B0)*x
        z = np.zeros(len(x))
        fld = (z, B1, z)
        return Field(fld, self).tilt()

    def front_kick(self, state):
        u = self.__U(state[:, IMAP['x']])
        state[:, IMAP['dK']] -= u*1e-6/self.pardict['KinEn0']

    def rear_kick(self, state):
        u = self.__U(state[:, IMAP['x']])
        state[:, IMAP['dK']] += u*1e-6/self.pardict['KinEn0']

    def kickVolts(self, x):
        return (self.__volt, self.__U(x))


class ERF(Element):
    """RF element."""

    def __init__(self, length, reference_particle, acc_length, E_field=15e5, phase=1.5*np.pi, H_number=50, name="RF"):
        super().__init__(curve=0, length=length, name=name)

        if length == 0:
            self._Element__bool_skip = True
            length = 5e-4

        self.reference_particle = reference_particle

        self.amplitude = E_field
        self.phase = phase
        self.freq = self.reference_particle.revolution_freq(acc_length) * H_number
        self.__H_number = H_number

        self._Element__chars = pds.DataFrame({'Amplitude':self.amplitude,
                                              'Frequency':self.freq, 'h-number': self.__H_number,
                                              'Phase':self.phase}, index=[self.name]).T

        self.__U = self.amplitude*length # length instead self.length for compatibility with length 0
        
#    @Element.skip.setter
#    def skip(self, value):
#        if not isinstance(value, bool):
#            raise ValueError('Must be boolean.')
#        self._Element__bool_skip = value

    def EField(self, arg):
        i_t, = index(arg, 't')
        t = arg[i_t]
        A = self.amplitude
        w = self.freq*2*np.pi
        phi = self.phase
        z = np.zeros(len(t))
        fld = (z, z, A*np.cos(w*t+phi))
        return Field(fld, self).tilt()

    def EField_prime_t(self, arg): # Es prime divided by time prime
        i_t, = index(arg, 't')
        t = arg[i_t]
        A = self.amplitude
        w = self.freq*2*np.pi
        phi = self.phase
        z = np.zeros(len(t))
        fld = (z, z, -A*w*np.sin(w*t+phi))
        return Field(fld, self).tilt()

    def advance(self, state):
        """ alternative to integration,
            supposed to mimick a discrete kick
        """
        w = self.freq*2*np.pi
        u = self.__U
        i_dK, i_s, i_t = index(state, 'dK', 's', 't')
        K = self.reference_particle.kin_nrg_0 * (1 + state[i_dK])

        state[i_dK] += u*np.cos(w*state[i_t]+self.phase)*1e-6/self.reference_particle.kin_nrg_0
        state[i_s] += self.length
        _, beta = self.reference_particle.GammaBeta(K)
        state[i_t] += self.length/beta/pcl.CLIGHT

    def front_kick(self, state):
        u = self.__U
        state[:, IMAP['dK']] -= u*1e-6/self.reference_particle.kin_nrg_0

    def rear_kick(self, state):
        u = self.__U
        state[:, IMAP['dK']] += u*1e-6/self.reference_particle.kin_nrg_0

    def kickVolts(self):
        return self.__U

class Observer:
    """This element of zero length is put where we want
    to check the state vector.
    """
    
    def __init__(self):
        self.__bool_skip = True
        
    @property
    def skip(self):
        return self.__bool_skip
    
    def advance(self, state):
        """In Tracker::track, when an element's skip is true, it uses 
        the element's advance method. Hence it is defined here.
        """
        pass
        
        
#%%

if __name__ == '__main__':
    from particle import Particle
    el = Element(0,1)
    erf = ERF(3,Particle(), 5)
    erf1 = ERF(0,Particle(), 5)
    obs = Observer()