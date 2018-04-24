import numpy as np
import matplotlib.pyplot as plt

EZERO = 1.62E-19
CLIGHT = 2.9997E8

class Particle:
    def __init__(self, mass0=1876, gamma0=1.14):
        self._gamma0 = gamma0
        self._beta0 = np.sqrt(gamma0**2-1)/gamma0
        self._mass0 = mass0

        MeV2J = EZERO*1e6

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

    state = np.array([np.linspace(-1e-3, 1e-3, 5),
                      np.repeat(1e-3,5),
                      np.zeros(5),
                      np.random.normal(0,1e-3,5),
                      np.zeros(5),
                      np.random.normal(0, 1e-4,5)])

    Om = O._matrix
    Fm = F._matrix
    Dm = D._matrix
    M = Om.dot(Dm.dot(Om.dot(Fm)))

    state_list = [state]
    for i in range(10):
        state = M.dot(state)
        state_list.append(state)


    x = np.zeros(len(state_list))
    y = np.zeros_like(x)
    px = np.zeros_like(x)
    py = np.zeros_like(x)
    z = np.zeros_like(x)
    d = np.zeros_like(x)
    for n, s in enumerate(state_list):
        
        
    print("test success!")
