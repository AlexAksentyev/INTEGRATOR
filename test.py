import numpy as NP


StateVars = ['x','y']
SVM = {'x':0,'y':1}


class Element:
    
    def __init__(self, Name, Length, force):
        self.name = Name
        self.length = Length
        self.force = force
        
    def frontKick(self, state):
        i_y = SVM['y']
        n = len(SVM)
        i = NP.arange(i_y, len(state), n)
        state[i] -= .1
        
    def rearKick(self, state):
        i_y = SVM['y']
        n = len(SVM)
        i = NP.arange(i_y, len(state), n)
        state[i] += .1



class Ensemble:
    
    def __init__(self, state_list):
        self.w = 3.
        
        self.n_ics = len(state_list)
        self.n_var = len(state_list[0])
        
        self.ics = dict(zip(range(len(state_list)), state_list))
        
    def __RHS(self, state, at, field):
        x,y = state.reshape(self.n_var,self.n_ics)
        
        n = len(SVM)
        
        i_x = NP.arange(SVM['x'], len(state), n)
        i_y = NP.arange(SVM['y'], len(state), n)
        state = {'x':state[i_x], 'y':state[i_y]}
        
        w = self.w
        
        n = self.n_var*self.n_ics
        
        return NP.reshape([y, -w**2*x +field.force(at, state)], n)
    
    def __getitem__(self, pid):
        return getattr(self, 'log'+str(pid))
    
    def plot(self, Ylab, Xlab='t',**kwargs):
        from matplotlib import pyplot as PLT
        
        for pid in self.ics.keys():
            PLT.plot(self[str(pid)][Xlab], self[str(pid)][Ylab],label=pid,**kwargs)
            
        PLT.xlabel(Xlab)
        PLT.ylabel(Ylab)
        PLT.legend()

    def track(self, FldSeq , ntimes):
        brks = 101
        
        ics = list()
        
        names = ['START']+[e.name for e in FldSeq]
        n = str(len(names[NP.argmax(names)]))
        EType = 'U'+n
        vartype = [('Turn',int),('Element',EType),('t', float), ('x', float), ('y', float)]
        
        self.__fLastPnt = -1
        
        nrow = ntimes*len(FldSeq)*brks
        
        for pid, ic in self.ics.items():
            setattr(self, 'log'+str(pid), NP.recarray(nrow,dtype=vartype))
            ics.append(ic)
        
        ics = NP.array(ics)
        ind = 0
        
        n_ics = self.n_ics
        n_var = self.n_var
        state = ics.reshape(n_ics*n_var)
        
        t0 = 0
        for n in range(1,ntimes+1):
            for i in range(len(FldSeq)):
                percent = int(ind/nrow*100)
                if (percent%10==0):
                    print('Complete {} %'.format(percent))
                    
                element = FldSeq[i]
                t1 = t0 + element.length
                at = NP.linspace(t0, t1, brks)
                t0 = t1
                
                element.frontKick(state)
                vals = odeint(self.__RHS, state, at, args=(element,))
                for k in range(brks):
                    element.rearKick(vals[k])
                    valsk = vals[k].reshape(n_ics, n_var, order='F')
                    for pid in self.ics.keys():
                        log = self[pid]
                        log[ind] = n,element.name, at[k], *valsk[pid]
                    ind += 1
                
        print('Complete 100 %')
        
        
#%%
if __name__ is '__main__':
    from scipy.integrate import odeint
    
    n = 10
    elist = NP.empty(n, dtype=Element)
    for i in range(10):
        A = 0
        w = 3
        phi = NP.pi/9*i
        L = 5
        Name = 'EL_'+str(i)
        f = lambda t, state: A*NP.cos(w*t+phi) - .0*state['y']
        elist[i] = Element(Name, L, f)
        
        
#%%

    state0 = NP.array([0,1], float)
    state1 = NP.array([0,-1], float)
    state2 = NP.array([0,0], float)
    states = NP.array([state0, state1, state2])
    E = Ensemble(states)

    E.track(elist, 5)
    E.plot('y')
#    E.plot('y','x',linewidth=.05)
    