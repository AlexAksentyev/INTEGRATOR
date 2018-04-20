
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
from utilities import MutableNamedTuple
from scipy.integrate import odeint

import particle as pcl

from rhs import IMAP, index, select, VAR_NUM

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

#%%
class Tilt:
    def __init__(self):
        self.angle = {'X':0, 'S':0} # angles in rads
        self.__set_rep_str()

    def __call__(self, order, *degree_angles, append=False):
        order = order.upper()

        find_axis = lambda axis_name: np.array([i for i, axis in enumerate(order) if axis == axis_name])

        # find the positions of 'X' and 'S', but not 'Y', b/c rotations about 'Y'
        # don't have any effect on the dipole field
        x_index = find_axis("X")
        s_index = find_axis("S")
        # remove the corresponding anglesand the corresponding angles
        angle_x = math.radians(np.sum(np.array(degree_angles)[x_index])) if x_index.size else 0
        angle_s = math.radians(np.sum(np.array(degree_angles)[s_index])) if s_index.size else 0

        if append:
            angle_x += self.angle['X']
            angle_s += self.angle['S']

        self.angle.update({'X':angle_x, 'S':angle_s})

        self.__set_rep_str()

        return angle_x, angle_s

    def __set_rep_str(self):
        ## to update
        axis = ['X', 'S']
        adeg = [ math.degrees(x) for x in [self.angle['X'], self.angle['S']] ]
        arad = [self.angle['X'], self.angle['S']]
        df = pds.DataFrame({'Axis': axis, 'Deg': adeg, 'Rad': arad})
        self.__rep_str = str(df)

    def __repr__(self):
        return self.__rep_str


class Shift(MutableNamedTuple):
    __slots__ = ['x_shift', 'y_shift']

    def __init__(self, x=0, y=0):
        super().__init__(x,y)

    def __add__(self, pair):
        self.x_shift += pair[0]
        self.y_shift += pair[1]
        return self
#%%

class Element:

    def __init__(self, curve, length, name="Element", **kwargs):
        self.curve = curve
        self.length = length
        self.name = name

        self.__E_field = (0, 0, 0)
        self.__B_field = (0, 0, 0)
        self.__delta_B_field = (0, 0, 0)
        self.__delta_E_field = (0, 0, 0)

        self.tilt_ = Tilt()
        self.shift_ = Shift()

        self.__bool_skip = False # for ERF.advance()

        self.__chars = pds.DataFrame({'Curve':self.curve, 'Length':self.length}, index=[self.name])

        super().__init__(**kwargs)

    def V(self, x):
        return np.zeros_like(x)

    def get_fields(self):
        return self.__E_field, self.__B_field

    def EField(self, arg):
        return Field(self.__E_field, self).vectorize(arg)

    def EField_prime_s(self, arg): # added for ERF
        return Field((0, 0, 0), self).vectorize(arg)

    def BField(self, arg):
        return Field(self.__B_field, self).vectorize(arg) 

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
        return self.tilt_(order, *args, append=append)

    def shift(self, x_shift=0, y_shift=0, append=False):
        if append:
            self.shift_ += (x_shift, y_shift)
        else:
            self.shift_ = Shift(x_shift, y_shift)

    def track(self, rhs, state, brks, controls_dict={}):
        """ Default element tracking behavior: integrate;
        if there exists a map for the specific element, use that instead;
        or whatever"""
        at = np.linspace(0, self.length, brks)
        vals = odeint(rhs, state, at, args=(self, ),**controls_dict)
        state[:] = vals[brks-1]
        return vals

#%%

class TiltableElement(Element):

    def __init__(self, reference_particle, **kwargs):
        _, self.ref_beta = reference_particle.GammaBeta()

        self.__delta_B_field = (0,0,0)
        self.__delta_E_field = (0,0,0)

        super().__init__(**kwargs)

    def tilt(self, order, *deg_angles, append=False):
        angle_x, angle_s = super().tilt(order, *deg_angles, append)
        tan_angle_x = math.tan(angle_x)
        tan_angle_s = math.tan(angle_s)
        By = self._Element__B_field[1]
        dBx = By*tan_angle_s
        v = self.ref_beta * pcl.CLIGHT
        self.__delta_B_field = (dBx, 0, By*tan_angle_x)
        self.__delta_E_field = (0, -v*dBx, 0)

    def BField(self, arg):
        pure_fld = super().BField(arg) # calls Element.BField
        delta = Field(self.__delta_B_field, self).vectorize(arg) # field from tilt
        return pure_fld + delta

    def Ey_comp(self, arg):
        """Compensatory E-field that we get after tilting to
        keep the reference orbit"""
        return Field(self.__delta_E_field, self).vectorize(arg)

    def EField(self, arg):
        pure_fld = super().EField(arg) # calls Element.EField
        return pure_fld + self.Ey_comp(arg)


#%%
class Drift(Element):
    """Drift space."""

    def __init__(self, length, name="Drift"):
        super().__init__(curve=0, length=length, name=name)
        self.map = None


class MQuad(Element):
    """Magnetic quadrupole."""

    def __init__(self, length, grad, name="MQuad"):
        super().__init__(curve=0, length=length, name=name)
        self.__grad = grad
        self._Element__chars['Grad'] = self.__grad

    def BField(self, arg):
        i_x, i_y = index(arg, 'x', 'y')
        x = arg[i_x]
        x0 = self.shift_.x_shift
        y0 = self.shift_.y_shift
        y = arg[i_y]
        fld = (self.__grad*y, self.__grad*x, np.zeros(len(x)))
        return  Field(fld, self) - self.__grad*np.array([[y0], [x0], [0]])


class MDipole(TiltableElement):
    """(Horizontally) Bending magnetic dipole."""

    def __init__(self, reference_particle, length, R=None, B_field=None, name="MDipole"):
        if all([e is None for e in [B_field, R]]):
            raise Exception("Either the B-field, or Radius must be defined")

        if R is None: # if R isn't given, B_field is
            if isinstance(B_field, cln.Sequence):
                By = B_field[1] # we want the By component
            else:
                By = B_field
            if By == 0:
                crv = 0 # if By == 0, dipole is turned off == drift space; skip computeRadius
                R = np.Inf
            else: crv = None # otherwise, computeRadius
        else:
            crv = 1/R

        super().__init__(reference_particle=reference_particle, curve=crv, length=length, name=name)

        if crv is None:
            R = reference_particle.Pc()*1e6/(B_field*pcl.CLIGHT)
            self.curve = 1/R
            self._Element__chars['Curve'] = self.curve

        if B_field is None:
            B_field = (0, reference_particle.Pc()*1e6/(R*pcl.CLIGHT), 0)
        elif not isinstance(B_field, cln.Sequence): B_field = (0, B_field, 0) # by default B = By

        self._Element__B_field = B_field[:]
        self._Element__chars['Bx'] = B_field[0]
        self._Element__chars['By'] = B_field[1]
        self._Element__chars['Bs'] = B_field[2]

    def get_B_field(self):
        return self._Element__B_field[:]
       

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
        return Field(fld, self)

class ECylDeflector0(Element):
    """Electrostatic cylindrical deflector, as defined in Andrey's matlab code.
    """

    def __init__(self, reference_particle, length, h_gap, R, name="ECylDefl"):
        self._ref_kinetic_energy = reference_particle.kinetic_energy

        crv = 1/R
        R1 = R - h_gap
        R2 = R*R/R1

        # Lorentz force is not necessarily zero
        # Reference bend radius R0 is computed via the
        # centrifugal force formula
        gamma, beta = reference_particle.GammaBeta()
        v = pcl.CLIGHT * beta
        mv2 = gamma*reference_particle.mass0_kg *  v**2
        E0 = mv2 / (- pcl.EZERO * R)
        V = (E0 * R * np.log(R2 / R1)) / (-2)

        super().__init__(curve=crv, length=length, name=name)

        self._U = lambda x: -V + 2*V*np.log((R+x)/R1)/np.log(R2/R1)
        self._Element__E_field = (E0, 0, 0)

    def V(self, x):
        return self._U(x)

    def EField(self, arg):
        x, = select(arg, 'x')
        Ex = self._Element__E_field[0]/(1+self.curve*x)

        z = np.zeros(len(x))
        fld = (Ex, z, z)
        return Field(fld, self)

    def front_kick(self, state):
        u = self._U(select(state, 'x')[0])
        i_dK, = index(state, 'dK')
        state[i_dk] -= u*1e-6/self._ref_kinetic_energy

    def rear_kick(self, state):
        u = self._U(select(state, 'x')[0])
        i_dK, = index(state, 'dK')
        state[i_dk] += u*1e-6/self._ref_kinetic_energy

class ECylDeflector(Element):
    """Electrostatic cylindrical deflector, as defined in Andrey's matlab code.
    Changed: instead of passing curvature radius R in the constructor, pass
    the radial electric field magnitude (signed).
    """

    def __init__(self, reference_particle, length, h_gap, E_field, name="ECylDefl"):
        self._ref_kinetic_energy = reference_particle.kinetic_energy

        # Lorentz force is not necessarily zero
        # Reference bend radius R0 is computed via the
        # centrifugal force formula
        gamma, beta = reference_particle.GammaBeta()
        v = pcl.CLIGHT * beta
        mv2 = gamma*reference_particle.mass0_kg *  v**2
        R = mv2 / (- pcl.EZERO * E_field )
        crv = 1/R

        R1 = R - h_gap
        R2 = R*R/R1
        V = (E_field * R * np.log(R2 / R1)) / (-2)

        super().__init__(curve=crv, length=length, name=name)

        self._U = lambda x: -V + 2*V*np.log((R+x)/R1)/np.log(R2/R1)
        self._Element__E_field = (E_field, 0, 0)

    def V(self, x):
        return self._U(x)

    def EField(self, arg):
        x, = select(arg, 'x')
        Ex = self._Element__E_field[0]/(1+self.curve*x)

        z = np.zeros(len(x))
        fld = (Ex, z, z)
        return Field(fld, self)

    def front_kick(self, state):
        u = self._U(select(state, 'x')[0])
        i_dK, = index(state, 'dK')
        state[i_dK] -= u*1e-6/self._ref_kinetic_energy

    def rear_kick(self, state):
        u = self._U(select(state, 'x')[0])
        i_dK, = index(state, 'dK')
        state[i_dK] += u*1e-6/self._ref_kinetic_energy

class CylWien(TiltableElement):
    """Cylindtical Wien filter.
    For the field definitnion confer:
        http://home.fnal.gov/~ostiguy/OptiM/OptimHelp/electrostatic_combined_function.html
    """

    def __init__(self, reference_particle, length, h_gap, E_field, B_field, name="Wien"):
        self._ref_kinetic_energy = reference_particle.kinetic_energy

        # Lorentz force is not necessarily zero
        # Reference bend radius R0 is computed via the
        # centrifugal force formula
        gamma, beta = reference_particle.GammaBeta()
        v = pcl.CLIGHT * beta
        mv2 = gamma*reference_particle.mass0_kg *  v**2
        R = mv2 / (- pcl.EZERO * (E_field - v*B_field))
        crv = 1/R
        if crv < 0:
            print("Negative curvature, bending in the direction -x.")
        else:
            print("Positive curvature, bending in the direction +x.")

        super().__init__(reference_particle=reference_particle, curve=crv, length=length, name=name)

        R1 = R - h_gap
        R2 = R*R/R1
        V = (E_field * R * np.log(R2 / R1)) / (-2)

        self._U = lambda x: -V + 2*V*np.log((R+x)/R1)/np.log(R2/R1)
        self._Element__B_field = (0, B_field, 0)
        self._Element__E_field = (E_field, 0, 0)

    def V(self, x):
        return self._U(x)

    def EField(self, arg):
        # redefinition of Element.EField
        x, = select(arg, 'x')
        Ex = self._Element__E_field[0]/(1+self.curve*x)

        z = np.zeros(len(x))
        pure_fld = (Ex, z, z)
        return Field(pure_fld, self) + self.Ey_comp(arg) # redefines TiltableElement.EField()

    def front_kick(self, state):
        u = self._U(select(state, 'x')[0])
        i_dK, = index(state, 'dK')
        state[i_dK] -= u*1e-6/self._ref_kinetic_energy

    def rear_kick(self, state):
        u = self._U(select(state, 'x')[0])
        i_dK, = index(state, 'dK')
        state[i_dK] += u*1e-6/self._ref_kinetic_energy


class StraightWien(TiltableElement):
    """Straight Wien filter."""

    def __init__(self, reference_particle, length, h_gap, E_field, B_field, name="Wien"):
        self._ref_kinetic_energy = reference_particle.kinetic_energy

        super().__init__(reference_particle=reference_particle, curve=0, length=length, name=name)

        self._U = h_gap*E_field
        self._Element__E_field = (E_field, 0, 0)
        self._Element__B_field = (0, B_field, 0)

    def V(self, x):
        return np.repeat(self._U, len(x))

    def front_kick(self, state):
        i_dK, = index(state, 'dK')
        state[i_dK] -= self._U*1e-6/self._ref_kinetic_energy

    def rear_kick(self, state):
        i_dK, = index(state, 'dK')
        state[i_dK] += self._U*1e-6/self._ref_kinetic_energy



class ERF(Element):
    """RF element."""

    def __init__(self, reference_particle, length, acc_length, E_field=15e5, phase=-0.5*np.pi, H_number=50, name="RF"):
        super().__init__(curve=0, length=length, name=name)

        if length == 0:
            self._Element__bool_skip = True
            length = 5e-4

        self.reference_particle = reference_particle

        self.amplitude = E_field
        self.phase = phase
        self.freq = self.reference_particle.revolution_freq(acc_length) * H_number
        _, beta = reference_particle.GammaBeta()
        self.wave_num = self.freq*2*np.pi/pcl.CLIGHT/beta
        self.__H_number = H_number

        self._Element__chars = pds.DataFrame({'Amplitude':self.amplitude,
                                              'Frequency':self.freq, 'h-number': self.__H_number,
                                              'Phase':self.phase}, index=[self.name]).T

        self.__U = self.amplitude*length # length instead self.length for compatibility with length 0

    def V(self, x):
        return np.repeat(self._U, len(x))

    def EField(self, arg):
        z, = select(arg, 'z')
        A = self.amplitude
        beta = self.reference_particle.beta
        w = self.freq*2*np.pi
        phase = -z*w/beta/CLIGHT + self.phase # delta phi + phi_ref
        z = np.zeros(len(s))
        fld = (z, z, A*np.cos(phase))
        return Field(fld, self)

    def EField_prime_s(self, arg): # Es prime 
        z, = select(arg, 'z')
        A = self.amplitude
        beta = self.reference_particle.beta
        w = self.freq*2*np.pi
        phase = -z*w/beta/CLIGHT + self.phase # delta phi + phi_ref
        z = np.zeros(len(s))
        fld = (z, z, A*w/betaa/CLIGHT*np.sin(phase))
        return Field(fld, self)

    def track(self, rhs, state, brks, controls_dict={}):
        """ Alternative to integration, supposed to mimick a discrete kick.
        Only uses the **state** argument.
        """
        w = self.freq*2*np.pi
        u = self.__U
        i_dK, i_z = index(state, 'dK', 'z')
        K = self.reference_particle.kinetic_energy * (1 + state[i_dK])
        _, beta = self.reference_particle.GammaBeta(K)
        dphi = -w/beta/pcl.CLIGHT*state[i_z]

        state[i_dK] += u*np.cos(+self.phase)*1e-6/self.reference_particle.kinetic_energy
        state[i_z] += self.length
        _, beta = self.reference_particle.GammaBeta(K)
        return state

    def front_kick(self, state):
        u = self.__U
        i_dK, = index(state, 'dK')
        state[i_dK] -= u*1e-6/self.reference_particle.kinetic_energy

    def rear_kick(self, state):
        u = self.__U
        i_dK, = index(state, 'dK')
        state[i_dK] += u*1e-6/self.reference_particle.kinetic_energy

class Observer(Element):
    """This element of zero length is put where we want
    to check the state vector.
    """

    def __init__(self, name='Observer'):
        super().__init__(0, 0, name)
        self._Element__bool_skip = True

    def V(self, x):
        return np.zeros_like(x)

    def track(self, rhs, state, brks, controls_dict={}):
        """In Tracker::track, when an element's skip is true, it uses
        the element's advance method. Hence it is defined here.
        """
        pass

# class EVert(Observer):
#     """Adds a vertical field Ey."""
#     def __init__(self, reference_particle, B_hor, name="Ey"):
#         super().__init__(name)
#         v = pcl.CLIGHT*self.ref_beta # WRONG GET BETA FROM REFERENCE PARTICLE
        
#         self._Element__E_field = (0, v*B_hor, 0)


#%%
# setup
if __name__ == '__main__':
    """test of CylWien"""

    from particle import Particle
    from tracker import Tracker
    from particle_log import StateList
    from lattice import Lattice
    import matplotlib.pyplot as plt

    deu = Particle()
    trkr = Tracker()
    trkr.set_controls(inner=True)

    bunch = StateList(Sz=1, x=[-1e-3, 1e-3], y=[-1e-3, 1e-3])
    
    element0 = MDipole(deu, 180e-2, B_field=1, name="DIP")
    element1 = MQuad(25e-2, 8.6, name="QF")

    state = np.array(bunch.as_list()).flatten()

    B0 = element1.BField(state)
    
    element1.shift(1e-3, -1e-3)

    B1 = element1.BField(state)
    element1.shift(-1e-3, 1e-3, True)

    B2 = element1.BField(state)
