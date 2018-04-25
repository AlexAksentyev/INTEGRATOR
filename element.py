import numpy as np

class Element:

    def __init__(self, particle, length, curve, name="Element"):
        self._length = length
        self._curve = curve
        self._name=name
        self._particle = particle
        self._matrix = None

    @property
    def M(self):
        return self._matrix

class Drift(Element):

    def __init__(self,particle, length, name="Drift"):
        super().__init__(particle, length, 0, name)

        beta0 = self._particle.beta
        gamma0 = self._particle.gamma
        b2g2 = beta0*beta0*gamma0*gamma0
        Mx = np.array([[1, length], [0, 1]])
        Z = np.zeros((2,2))
        Mz = np.array([[1, length/b2g2], [0, 1]])

        self._matrix = np.bmat([[Mx, Z, Z], [Z, Mx, Z], [Z, Z, Mz]])

    def __call__(self, state):
        state_c = state[:]
        px = state_c[1]
        py = state_c[3]
        delta = state_c[5]

        beta0 = self._particle.beta
        gamma0 = self._particle.gamma
        normed_E = 1/beta0 + delta

        root = np.sqrt(normed_E**2 - px**2 - py**2 - 1/(beta0*gamma0)**2)

        state_c[0] += px*self._length/root
        state_c[2] += py*self._length/root
        state_c[4] += 1/beta0 - normed_E/root

        return state_c

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
        beta0 = self._particle.beta
        gamma0 = self._particle.gamma
        b2g2 = beta0*beta0*gamma0*gamma0
        Mf = np.array([[c, s/w], [-w*s, c]])
        Md = np.array([[ch, sh/w], [w*sh, ch]])
        Mz = np.array([[1, L/b2g2], [0, 1]])
        Z = np.zeros((2,2))

        if grad > 0:
            self._matrix = np.bmat([[Mf, Z, Z], [Z, Md, Z], [Z, Z, Mz]])
        else:
            self._matrix = np.bmat([[Md, Z, Z], [Z, Mf, Z], [Z, Z, Mz]])
