from particle import Particle, EZERO, CLIGHT

import numpy as np
import matplotlib.pyplot as plt


VAR_NAME=['x','px','y','py','z','d']
VAR_NUM=len(VAR_NAME)

class PLog(np.recarray):
    """Particle log is a 3D recarray, where the 1st index refers to the particle (ics),
    the 2nd to the row (record), the 3rd to the variable (like x,y,s, &c).
    """

    max_len_name = 10
    el_field_type = object#tbl.StringCol(max_len_name) # otherwise problems writing into hdf5 file
    metadata_type = [('Turn', int)]
    variable_type = list(zip(VAR_NAME, np.repeat(float, VAR_NUM)))
    record_type = metadata_type + [('PID', int)] + variable_type
                ## PID field is weird like this to enable automatic
                ## pid setting in setitem (if it were metadata, i'd have to supply it,
                ## and that causes problems in __array_finalize__
                ## i could actually do without it, now that i know that subsetting
                ## a recarray makes it 1D, and I have to manually reshape it.

    _state_var_1st_ind = len(metadata_type) + 1 # len(metadata) + len(PID) - 1 + 1

    def __new__(cls, ini_states, n_records):

        if not ini_states.shape[0] == VAR_NUM:
            print('Wrong number of state variables')
            return

        n_ics = ini_states.shape[1]

        obj = np.recarray((n_records, n_ics), dtype=cls.record_type, order='C').view(cls)
        super(PLog, obj).fill(np.nan)
        obj.n_ics = n_ics
        obj.ics = ini_states

        obj[0] = ((0,), ini_states)

        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

        self.n_ics = getattr(obj, 'n_ics', None)
        self.ics = getattr(obj, 'ics', None)

    def trajectories(self, from_=0, to_=None):
        """Returns a generator of particle logs. That is,
        in the expression:
            for p in plog.trajectories():
                plot(p['s'], p['x'])
        p = plog[:, i] for the next i, and hence
        the code above will plot x vs s for all of the
        logged particles.
        """
        if to_ is None or to_ > self.n_ics:
            to_ = self.n_ics
        for pid in range(from_, to_):
            yield self[:, pid]

    def __setitem__(self, i, stamp_vector):
        """Scipy integrator returns a flat array of state vector values
        [x_i^0, y_i^0, ..., Sz_i^0, x_i^1, y_i^1, ...], where the superscript
        denotes the particle number.
        We want to append each particle's state variable vector with some
        metadata, and turn the resulting vector into a tuple.
        """
        stamp, vector = stamp_vector
        vector = vector.T

        stamped = np.empty(self.n_ics, self.dtype)
        for ind, vec in enumerate(vector):
            stamped[ind] = stamp + (ind,) + tuple(vec.A1)

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

        beta0 = self._particle.beta
        gamma0 = self._particle.gamma
        b2g2 = beta0*beta0*gamma0*gamma0
        Mx = np.array([[1, length], [0, 1]])
        Z = np.zeros((2,2))
        Mz = np.array([[1, length/b2g2], [0, 1]])

        self._matrix = np.bmat([[Mx, Z, Z], [Z, Mx, Z], [Z, Z, Mz]])

    def map(self, state):
        px = state[:,1]
        py = state[:,3]
        delta = state[:,5]

        beta0 = self._particle.beta
        gamma0 = self._particle.gamma
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
        beta0 = self._particle.beta
        gamma0 = self._particle.gamma
        b2g2 = beta0*beta0*gamma0*gamma0
        Mf = np.array([[c, s/w], [-w*s, c]]) # focusing matrix
        Md = np.array([[ch, sh/w], [w*sh, ch]]) # defocusing matrix
        Mz = np.array([[1, L/b2g2], [0, 1]])
        Z = np.zeros((2,2))

        if self._grad > 0:
            self._matrix = np.bmat([[Mf, Z, Z], [Z, Md, Z], [Z, Z, Mz]])
        else:
            self._matrix = np.bmat([[Md, Z, Z], [Z, Mf, Z], [Z, Z, Mz]])


if __name__ == '__main__':
    p = Particle()
    O = Drift(p, 25e-2)
    F = MQuad(p, 25e-2, 8.6)
    D = MQuad(p, 25e-2, -8.11)

    n_ics = 5
    state = np.matrix([np.linspace(-1e-3, 1e-3, n_ics),
                      np.repeat(1e-3,n_ics),
                      np.zeros(n_ics),
                      np.random.normal(0,1e-3,n_ics),
                      np.zeros(n_ics),
                      np.random.normal(0, 1e-4,n_ics)])

    Om = O._matrix
    Fm = F._matrix
    Dm = D._matrix
    M = Om.dot(Dm.dot(Om.dot(Fm)))

    n_trn = 10
    log = PLog(state, n_trn+1)
    for i in range(1, n_trn+1):
        state = M.dot(state)
        log[i] = ((i,), state)

    print("test success!")
