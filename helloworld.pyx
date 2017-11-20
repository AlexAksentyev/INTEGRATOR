# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 16:56:35 2017

@author: Аксентьев
"""
from scipy.integrate import odeint
from CElement import Element
#from CElement cimport Element

def primes(int kmax):
    cdef int n, k, i
    cdef int p[1000]
    result = []
    if kmax > 1000:
        kmax = 1000
    k = 0
    n = 2
    while k < kmax:
        i = 0
        while i < k and n % p[i] != 0:
            i = i + 1
        if i == k:
            p[k] = n
            k = k + 1
            result.append(n)
        n = n + 1
    return result

def f(double[:] x, double t, Element element):
#    Ex,Ey,Es = element.EField(x)
#    Bx,By,Bs = element.BField(x)
    
    km = -1.13
    
    return [x[1], km*x[0]]

#def track(Element[:] ElSeq, int ntimes, double[:] state):
#    
#    cdef float brks = 101
#    
#    cdef:
#        int n,i,k,ind = 0
#        double[:] at = NP.empty([brks])
#        double[:] vals = NP.empty([brks])
#        
#        vartype = [('Turn',int),('Element',object),('Point', object)]
#        vartype += list(zip(StateVars, NP.repeat(float, len(StateVars))))
#        
#        nrow = ntimes*len(ElementSeq)*self.fIntBrks
#        fStateLog = NP.recarray(nrow,dtype=vartype)
#        
#    for n in range(1,ntimes):
#        for i in range(len(ElSeq)):
#            element = ElSeq[i]
#            at = NP.linspace(0, element.fLength, brks)
#            vals = odeint(f, state, at, args=(element,))
#            state = vals[brks-1]
#            for k in range(brks-1):
#                fStateLog[ind] = n,element.fName,k, *vals[k]
#                ind += 1