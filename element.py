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

        kx = np.sqrt(grad + 0j) #turn grad complex to take sqrt if negative
        ky = 1j*kx
        L = self._length

        cx = np.cos(kx*L)
        sx = np.sin(kx*L)/kx 
        cy = np.cos(ky*L) 
        sy = np.sin(ky*L)/ky
        
        beta0 = self._particle.beta
        gamma0 = self._particle.gamma
        b2g2 = beta0*beta0*gamma0*gamma0
        
        Mx = np.array([[cx, sx], [-grad*sx, cx]])
        My = np.array([[cy, sy], [ grad*sy, cy]])
        Mz = np.array([[1, L/b2g2], [0, 1]])
        Z = np.zeros((2,2))

        def check_complex(matrix):
            if np.all(matrix.imag == 0):
                matrix = matrix.real
            else:
                print("WARNING: Complex transfer matrix!")
            return matrix

        Mx = check_complex(Mx)
        My = check_complex(My)

        self._matrix = np.bmat([[Mx, Z, Z], [Z, My, Z], [Z, Z, Mz]])

# class RF(Element):

#     def __init__(self, particle, acc_length, ref_phase, h_number):
#         beta0 = self._particle.beta
#         gamma0 = self._particle.gamma
#         Ts = 
