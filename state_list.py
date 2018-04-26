from plog import VAR_NAME, VAR_NUM, IMAP
from collections import Iterable
import numpy as np

class StateList:
    """Create an ensemble of initial conditions.
    Keyword arguments are the state variables:
        :x, y:       the particle position in the transverse plane
        :s:           position along the optical axis
        :t:           time of injection
        :Theta:       RF phase at injection
        :H:           sqrt(x^2 + y^2 + s^2)
        :px, py:      the reference-momentum-normalized momentum projections Px/P0, Py/P0
        :dK:          normalized energy deviation (K-K0)/K0
        :Sx, Sy, Sz:  spin projections.

    Use
    ----------------
    To define an ensemble of initial conditions varying in the variable **x**, write: \n
        StateList(x=[x0, x1, x2, ..., xn], ...)

    To set a variable **y** constant accross all states, write: \n
        StateList(..., y = shared_value)

    """
    def __init__(self, **kwargs):

        keys = kwargs.keys()

        # create defined variables
        ntot = 1
        arg_dict = dict()
        for key, val in kwargs.items():
            if isinstance(val, Iterable):
                num = len(val)
            else: # passed a single value, cannot use len()
                num = 1
            ntot *= num
            arg_dict.update({IMAP[key]: val})

        # make mesh
        mesh = dict(zip(keys, np.meshgrid(*list(arg_dict.values()))))

        vartype = list(zip(VAR_NAME, np.repeat(float, VAR_NUM)))
        self.state_list = np.zeros(ntot+1, dtype=vartype) # +1 for genuine refrence particle

        #write data
        for key, value in mesh.items():
            self.state_list[key] = np.array([0]+value.reshape(ntot).tolist())

        # self.state_list[0]['Sz'] = 1

        #here i should remove any duplicate reference particles
#        self.state_list = np.unique(self.state_list) # this works but messes up the order

        # convert to list of dicts for use with ensemble
        self.state_list = [dict(zip(self.state_list.dtype.names, x)) for x in self.state_list]

    @classmethod
    def from_list(cls, state_list):
        if not isinstance(state_list, list):
            print('Need a list!')
            return
        if len(state_list[0]) != VAR_NUM:
            print('Wrong state vector length')
            return
        if isinstance(state_list[0], list):
            state_list = [dict(zip(VAR_NAME, x)) for x in state_list]
            print('Converted to a list of dicts.')

        result = cls()
        result.state_list = state_list

        return result

    def __len__(self):
        return len(self.state_list)

    def __getitem__(self, pid):
        return self.state_list[pid]

    def pop(self, index):
        self.state_list.pop(index)

    def __repr__(self):
        from pandas import DataFrame
        return str(DataFrame(self.state_list))

    def as_list(self):
        states = list()
        for d in self.state_list:
            states.append(list(d.values()))
        return states

    @property
    def array(self):
        return np.array(self.as_list()).T

    def write_to_file(self, filename, directory):
        from pandas import DataFrame
        from rhs import VAR_NAME
        DataFrame(self.state_list).to_csv(directory+'/'+filename,
                 columns=VAR_NAME, index=False)
