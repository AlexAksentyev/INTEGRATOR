import numpy as np
import matplotlib.pyplot as plt

EZERO = 1.62E-19
CLIGHT = 2.9997E8

VAR_NAME=['x','px','y','py','z','d']

class Particle:
    def __init__(self, mass0=1876, gamma0=1.14):
        self._gamma0 = gamma0
        self._beta0 = np.sqrt(gamma0**2-1)/gamma0
        self._mass0 = mass0

        MeV2J = EZERO*1e6

class PLog(np.recarray):
    rec_type = list(zip(VAR_NAME, np.repeat(float, 6)))

    def __new__(cls, ini_states, n_rec):
        n_ics = ini_states.shape[1]
        obj = np.recarray((n_rec, n_ics), dtype=cls.rec_type, order='C').view(cls)
        super(PLog, obj).fill(np.nan)
        
        obj[0] = ini_states
        obj.n_ics = n_ics
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

        self.n_ics = getattr(obj, 'n_ics', None)

    def __setitem__(self, item):
        vector = np.empty(self.n_ics, self.dtype)
        for ind, vec in enumerate(vector):
            stamped[ind] = stamp + (ind,) + tuple(vec)

        super(PLog, self).__setitem__(i, stamped)

class Element:

    def __init__(self, particle, length, curve, name="Element"):
        self._length = length
        self._curve = curve
        self._name=name
        self._particle = particle

class Drift(Element):

    def __init__(self,particle, length, name="Drift"):
        super().__init__(particle, length, 0, name)

        beta0 = self._particle._beta0
        gamma0 = self._particle._gamma0
        b2g2 = beta0*beta0*gamma0*gamma0
        Mx = np.array([[1, length], [0, 1]])
        Z = np.zeros((2,2))
        Mz = np.array([[1, length/b2g2], [0, 1]])

        self._matrix = np.array([[Mx, Z, Z], [Z, Mx, Z], [Z, Z, Mz]])
        self._matrix.shape = (6,6)

    def map(self, state):
        px = state[:,1]
        py = state[:,3]
        delta = state[:,5]

        beta0 = self._particle._beta0
        gamma0 = self._particle._gamma0
        normed_E = 1/beta0 + delta

        root = np.sqrt(normed_E**2 - px**2 - py**2 - 1/(beta0*gamma0)**2)

        state[:,0] += px*self._length/root
        state[:,2] += py*self._length/root
        state[:,4] += 1/beta0 - normed_E/root

class MQuad(Element):

    def __init__(self, particle, length, grad, name="MQ"):
        super().__init__(particle, length, 0, name)
        self._grad = grad

        w = np.sqrt(abs(grad))
        L = self._length
        wL = w*L
        c = np.cos(wL)
        s = np.sin(wL)
        ch = np.cosh(wL)
        sh = np.sinh(wL)
        beta0 = self._particle._beta0
        gamma0 = self._particle._gamma0
        b2g2 = beta0*beta0*gamma0*gamma0
        Mf = np.array([[c, s/w], [-w*s, c]]) # focusing matrix
        Md = np.array([[ch, sh/w], [w*sh, ch]]) # defocusing matrix
        Mz = np.array([[1, L/b2g2], [0, 1]])
        Z = np.zeros((2,2))

        if self._grad > 0:
            self._matrix = np.array([[Mf, Z, Z], [Z, Md, Z], [Z, Z, Mz]])
        else:
            self._matrix = np.array([[Md, Z, Z], [Z, Mf, Z], [Z, Z, Mz]])
        self._matrix.shape = (6,6)


if __name__ == '__main__':
    p = Particle()
    O = Drift(p, 25e-2)
    F = MQuad(p, 25e-2, 8.6)
    D = MQuad(p, 25e-2, -8.11)

    n_ics = 5
    state = np.array([np.linspace(-1e-3, 1e-3, n_ics),
                      np.repeat(1e-3,n_ics),
                      np.zeros(n_ics),
                      np.random.normal(0,1e-3,n_ics),
                      np.zeros(n_ics),
                      np.random.normal(0, 1e-4,n_ics)])

    Om = O._matrix
    Fm = F._matrix
    Dm = D._matrix
    M = Om.dot(Dm.dot(Om.dot(Fm)))

    state_list = [state]
    for i in range(n_trn):
        state = M.dot(state)
        state_list.append(state)
        
    print("test success!")
