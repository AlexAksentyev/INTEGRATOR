import numpy as np
from particle import CLIGHT, EZERO
from copy import deepcopy

MeV2J = EZERO*1e6

class Element:

    def __init__(self, particle, length, curve, name="Element"):
        self._length = length
        self._curve = curve
        self._name=name
        self._particle = particle
        I = np.eye(6)
        self._matrix = I
        self._tilt_matrix = I
        self._tilt_matrix_inv = I

        self._beta0 = self._particle.beta
        self._gamma0 = self._particle.gamma
        self._b2g2 = self._beta0*self._beta0*self._gamma0*self._gamma0
        

    @property
    def TM(self):
        return self._tilt_matrix*self._matrix*self._tilt_matrix_inv

    @property
    def length(self):
        return self._length

    @property
    def name(self):
        return self._name

    def __repr__(self):
        return self.name

    def __call__(self, state):
        """Applies element map to the state matrix"""
        # by default __call__ uses a linear transformation of the state matrix;
        # other elements may overwrite this behavior to use non-linear mapping
        return self.TM*state
    
    def s_tilt(self, angle):
        """Give angle in radians."""
        c, s = np.cos(angle), np.sin(angle)
        I = np.eye(2)
        Z = np.zeros((2,2))
        self._tilt_matrix = np.bmat([[c*I, s*I, Z], [-s*I, c*I, Z], [Z, Z, I]])
        self._tilt_matrix_inv = np.linalg.inv(self._tilt_matrix)

    def clear_tilt(self):
        self._tilt_matrix = np.eye(6)
        self._tilt_matrix_inv = np.eye(6)

class Drift(Element):

    def __init__(self,particle, length, name="Drift"):
        super().__init__(particle, length, 0, name)

        b2g2 = self._b2g2
        Mx = np.array([[1, length], [0, 1]])
        Z = np.zeros((2,2))
        Mz = np.array([[1, length/b2g2], [0, 1]])

        self._matrix = np.bmat([[Mx, Z, Z], [Z, Mx, Z], [Z, Z, Mz]])

    def __call__(self, state):
        """This map is taken from the Cockroft lectures, NOT the MAD manual.
        It is different; even if we assume delta_s = 0 (reference Ps == design P0).
        The matrix is from the MAD manual, as for all elements.
        """
        state_c = np.array(deepcopy(state))
        px = state_c[1]
        py = state_c[3]
        delta = state_c[5]

        beta0 = self._beta0
        gamma0 = self._gamma0
        normed_E = 1/beta0 + delta

        root = np.sqrt(normed_E**2 - px**2 - py**2 - 1/self._b2g2)

        state_c[0] += px*self._length/root
        state_c[2] += py*self._length/root
        state_c[4] += 1/beta0 - normed_E/root

        return state_c

class MDipoleSect(Element):
    """Sector magnetic dipole. Implements only the body transfer map.
    NO FRINGE FIELDS.
    """
    def __init__(self, particle, length, B_field, name="MDS"):
        P0c = particle.Pc()
        P0 = P0c*MeV2J/CLIGHT
        k0 = particle.charge*B_field/P0
        super().__init__(particle, length, k0, name)
        # http://pcwww.liv.ac.uk/~awolski/Teaching/Cockcroft/LinearDynamics/LinearDynamics-Lecture3.pdf [p 32]:
        w = h = k0
        L = length
        c = np.cos(w*L)
        s = np.sin(w*L)
        d = (1 - c)/w
        f = (w*L - s)/w
        b0 = self._beta0
        b2g2 = self._b2g2

        Z = np.zeros((2,2))
        Mxx = np.array([[c, s/w], [-w*s, c]])
        Mxy = Z
        Mxz = np.array([[0, d/b0], [0, s/b0]])
        Myx = Z
        Myy = np.array([[1, L], [0, 1]])
        Myz = Z
        Mzx = np.array([[-s/b0, -d/b0], [0, 0]])
        Mzy = Z
        Mzz = np.array([[1, L/b2g2 - f/b0/b0], [0, 1]])
        
        self._matrix = np.bmat([[Mxx, Mxy, Mxz],
                                [Myx, Myy, Myz],
                                [Mzx, Mzy, Mzz]])
        
    
class MQuad(Element):

    def __init__(self, particle, length, grad, name="MQ"):
        super().__init__(particle, length, 0, name)
        self._grad = grad

        P0 = particle.Pc()*MeV2J/CLIGHT
        k1 = particle.charge*grad/P0
        w = np.sqrt(abs(k1))
        L = self._length
        wL = w*L
        
        c = np.cos(wL)
        s = np.sin(wL)
        ch = np.cosh(wL)
        sh = np.sinh(wL)

        Mt = np.array([[c, s/w], [-w*s, c]])
        Mh = np.array([[ch, sh/w], [w*sh, ch]])
        Mz = np.array([[1, L/self._b2g2], [0, 1]])
        Z = np.zeros((2,2))
        
        if k1>0:
            self._matrix = np.bmat([[Mt, Z, Z], [Z, Mh, Z], [Z, Z, Mz]])
        else:
            self._matrix = np.bmat([[Mh, Z, Z], [Z, Mt, Z], [Z, Z, Mz]])

class MSext(Element):

    def __init__(self, particle, length, grad, name="SXT"):
        super().__init__(particle, length, 0, name)
        self._grad = grad

        b2g2 = self._b2g2
        Mx = np.array([[1, length], [0, 1]])
        Z = np.zeros((2,2))
        Mz = np.array([[1, length/b2g2], [0, 1]])

        self._matrix = np.bmat([[Mx, Z, Z], [Z, Mx, Z], [Z, Z, Mz]])

    def __call__(self, state):
        state_c = np.array(deepcopy(state))
        x, px, y, py, t, pt = state_c
        beta_s = self._beta0
        gamma_s = self._gamma0
        b2g2 = self._b2g2
        L = self._length
        K2 = self._grad
        L2 = L*L
        L3 = L*L*L
        L4 = L*L*L*L/24
        px2_py2 = px*px - py*py
        xpx_ypy = x*px - y*py
        xpy_ypx = x*py + y*px
        x2_y2 = x*x - y*y
        ftr = L*(1 - pt/beta_s)

        state_c[0] += ftr*px - K2*(L2/4*x2_y2 + L3/12*xpx_ypy + L4/24*px2_py2) - .5*L/beta_s * px*pt
        state_c[1] += -K2*(L/2*x2_y2 + L2/4*xpx_ypy + L3/6*px2_py2)
        state_c[2] += ftr*py + K2*(L2/4*x*y + L3/12*xpy_ypx + L4/24*px*py) - .5*L/beta_s * py*pt
        state_c[3] += K2*(L/2*x*y + L2/4*xpy_ypx + L3/6*px*py)
        # in the equation for t2 from page 25 of the MAD manual there's a term L * eta*delta_s/beta_s,
        # which i set to 0 here (beta_s and gamma_s in my equations refer to the reference particle,
        # and delta_s is 0 there; besides, I don't know eta!) !!
        # this is equivalent to using the equations from Andy Wolski, i think...
        state_c[4] += L/b2g2 - .5*L/beta_s*(px*px + py*py + 3*pt*pt/b2g2)

        return state_c

class RF(Element):

    def __init__(self, particle, acc_length, voltage, ref_phase=np.pi/2, h_number=50, name="RF"):
        super().__init__(particle, 0, 0, name)
        Ts = acc_length/CLIGHT/self._beta0
        w = 2*np.pi/Ts*h_number
        P0c = particle.Pc()

        Z = particle.Z
        A = Z*voltage*1e-6/P0c # qV/Pc w/energy in MeVs
        self._w = w
        self._ref_phase = ref_phase
        self._voltage = voltage
        self._A = A

        I = np.eye(2)
        Z = np.zeros((2,2))
        Mz = np.array([[1, 0],[-w*A*np.cos(ref_phase), 1]])

        self._matrix = np.bmat([[I, Z, Z],[Z, I, Z],[Z, Z, Mz]])

    def __call__(self, state):
        state_c = deepcopy(state)
        t = state_c[4]/CLIGHT
        state_c[5] += self._A*np.sin(self._ref_phase - self._w*t)
        return state_c

class CylWien(Element):
    def __init__(self, particle, length, E_field, B_field, name="CWF"):
        gamma, beta = particle.GammaBeta()
        g2 = gamma*gamma
        v = CLIGHT * beta
        mv2 = gamma*particle.mass0_kg *  v*v
        ## the expression for R0 is different from Optim, where it is:
        ## 1/R0 = EZERO * E_field / (beta*CLIGHT) + EZERO * B_field/ CLIGHT
        R0 = mv2 / (-EZERO * (E_field - v*B_field) )
        super().__init__(particle, length, 1/R0, name)

        ## matrix from OptimX:
        ## http://home.fnal.gov/~ostiguy/OptiM/OptimHelp/electrostatic_combined_function.html
        term = EZERO*E_field/mv2/gamma**2
        kappa_p = 1 + (term*R0/gamma)**2
        # gradients G_B, G_E = 0
        kx = 1/R0**2 + term**2 
        ky = 0
        kx2 = kx*kx
        kpR0 = kappa_p/R0
        kpkx = kappa_p/kx

        L = length
        cx = np.cos(kx*L)
        sx = np.sin(kx*L)
        dx = 1 - cx
        cy = 1; sy = 0 # b/c ky == 0
        
        Z = np.zeros((2,2))
        Mxx = np.array([[cx, sx/kx], [-kx*sx, cx]])
        Mxy = Z
        Mxz = np.array([[0, kpR0*dx/kx2], [0, kpR0/kx*sx]])
        Myx = Z
        Myy = np.array([[1, 1], [0, 1]])
        Myz = Z
        Mzx = np.array([[kpR0/kx*sx, kpR0/kx2*dx], [0, 0]])
        Mzy = Z
        Mzz = np.array([[1, L/g2 - kpR0*kpkx/kx*(L*kx - sx)], [0, 1]])

        self._matrix = np.bmat([[Mxx, Mxy, Mxz], [Myx, Myy, Myz], [Mzx, Mzy, Mzz]])
        
class Observer(Element):
    def __init__(self, particle, name="OBS"):
        super().__init__(particle, 0, 0, name)
