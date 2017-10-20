# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 21:40:11 2017

@author: Аксентьев
"""

def fun(f, *arg):
    s = f + '('
    order = 0
    for e in arg:
        if type(e) is not Symbolic: 
                    e = str(e)
                    e = Symbolic(e,e,0)
        s += e.S(False) + ','
        order = max(e.fOrder, order)
    s = s[0:len(s)-1] + ')'
    return Symbolic(s,s,order)

def exp(arg): return fun('exp', arg)
def sqrt(arg): return fun('sqrt',arg)
def pow(arg, p): return fun('pow',arg, p)


class Symbolic:
    
    def __init__(self, Def, Symb, Order, Val = None):
        self.fDef = Def
        self.fSymb = Symb
        self.fOrder = Order
        self.fVal = Val
        
    @classmethod
    def from_symb(cls, Symbolic):
        Def = Symbolic.fDef
        Symb = Symbolic.fSymb
        Order = Symbolic.fOrder
        
        return cls(Def, Symb, Order)
    
    def S(self, append=True):
        res = self.fSymb
        if append: res = 'v'+ str(self.fOrder) + '_' + res
        return res
    
    def D(self): return self.fDef   

    def V(self): return self.fVal
    
    def O(self): return self.fOrder
    
    @staticmethod
    def phi(operation , *w):
        sdef = '('
        ssymb = '('
        order = 0
        for e in w: 
            if type(e) is not Symbolic: 
                try:
                    e1 = float(e)
                    e = str(e)
                    e = Symbolic(e,e,0,e1)
                    if e1 < 0: e.enclose()
                except ValueError:
                    pass
            
            sdef += e.fDef + operation
            ssymb += e.fSymb + operation
            order = max(order, e.fOrder)
            
        sdef = sdef[0:len(sdef)-1] + ')'
        ssymb = ssymb[0:len(ssymb)-1] + ')'
        return Symbolic(sdef,ssymb,order)
    
    def enclose(self):
        self.fDef = '('+self.fDef+')'
        self.fSymb = '(' + self.fSymb+')'
        
    def __add__(self, other):
        return Symbolic.phi('+', self, other)
    
    def __radd__(self, other):
        return Symbolic.phi('+', other, self)
        
    def __sub__(self, other):
        return Symbolic.phi('-', self, other)
    
    def __rsub__(self, other):
        return Symbolic.phi('-', other, self)
        
    def __mul__(self, other):
        return Symbolic.phi('*', self, other)
    
    def __rmul__(self, other):
        return Symbolic.phi('*', other, self)
        
    def __truediv__(self,other):
        return Symbolic.phi('/', self, other)
    
    def __rtruediv__(self,other):
        return Symbolic.phi('/', other, self)
    
    def redefine(self, Def):
        self.fDef = Def
    
    def __repr__(self):
        return str({'Def':self.fDef, 'Symb':self.fSymb, 'Order':self.fOrder, 'Val':self.fVal})
    
    def __call__(self, *arg):
        s = '('
        order = 0
        for e in arg:
            if type(e) is not Symbolic: 
                e = str(e)
                e = Symbolic(e,e,0)
            s += e.S(False) + ','
            order = max(order, e.fOrder)
            
        s = self.S(False) + s[0:len(s)-1] + ')'
        
        return Symbolic(s,s,order)

        
if __name__ is '__main__':
    
    Mass0 = Symbolic('Mass0','Mass0',0,1876)
    dK = Symbolic('dK','dK',0)
    KinEn = Symbolic('KinEn','KinEn',0)
    KinEn0 = KinEn(0); KinEn0.fSymb = 'KinEn0'
    KinEn.redefine((KinEn0*(1+dK)).S(False))
    Pc = sqrt(pow(Mass0+KinEn,2)-pow(Mass0,2))