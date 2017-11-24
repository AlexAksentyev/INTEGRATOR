import numpy as NP



class Element:
    
    def __init__(self, Name, Length, force):
        self.name = Name
        self.length = Length
        self.force = force



class Particle:
    
    def __init__(self, mass, rigid):
        self.mass = mass
        self.rigid = rigid
        
    def __RHS(self, state, at, field):
        x,y = state
        
        k = self.rigid
        m = self.mass
        
        return [y, k/m*x + field.force(at)]
    
    def __getitem__(self, name):
        return self.log[name]
    
    def plot(self, Ylab, Xlab='t',**):
        from matplotlib import pyplot as PLT
        
        PLT.plot(self[Xlab], self[Ylab],**kwargs)
        PLT.xlabel(Xlab)
        PLT.ylabel(Ylab)

    def track(self, state0, FldSeq , ntimes):
        brks = 101
        
        names = ['START']+[e.name for e in FldSeq]
        n = str(len(names[NP.argmax(names)]))
        EType = 'U'+n
        vartype = [('Turn',int),('Element',EType),('t', float), ('x', float), ('y', float)]
        
        self.__fLastPnt = -1
        
        nrow = ntimes*len(FldSeq)*brks
        self.log = NP.recarray(nrow,dtype=vartype)
        ind = 0
        
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
                
                vals = odeint(self.__RHS, state0, at, args=(element,))
                for k in range(brks):
                    self.log[ind] = n,element.name, at[k], *vals[k]
                    ind += 1
                
        print('Complete 100 %')
        
        
#%%
if __name__ is '__main__':
    from scipy.integrate import odeint
    
    n = 10
    elist = NP.empty(n, dtype=Element)
    for i in range(10):
        A = 10
        w = 3
        phi = NP.pi/9*i
        L = 5
        Name = 'EL_'+str(i)
        f = lambda t: A*NP.cos(w*t+phi)
        elist[i] = Element(Name, L, f)
        
    
    p = Particle(10, 1.13)
        
#%%
    p.track([0,.1],elist, 10)
    p.plot('x')
    
    
    #%%
    x = NP.arange(0,1e6,.1)
    start = clock()
    y = NP.sin(x)
    ntime = clock() - start
    print('Time passed {} sec'.format(ntime))
    
    start = clock()
    y = NP.empty(len(x))
    for i in range(len(x)):
        y[i] = NP.sin(x[i])
        
    ntime_loop = clock() - start
    print('Loop time {} sec'.format(ntime_loop))
    
    print('Ratio: {}'.format(ntime_loop/ntime))