#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 14:53:02 2017

@author: alexa
"""

class Base:
    
    def __init__(self, name):
        self.fName = name
        
class DerivedA(Base):
    fCount = 0
    def __init__(self, name):
        self.fName = (name,DerivedA.fCount)
        DerivedA.fCount +=1
        
class DerivedB(Base):
    fCount = 0
    def __init__(self, name):
        self.fName = (name, DerivedB.fCount)
        DerivedB.fCount +=1
        

if __name__ == '__main__':
    
    print("None: "+str(DerivedA.fCount))
    print("None: "+str(DerivedB.fCount))
    first = DerivedA("Johnny")
    print("+A: "+str(DerivedA.fCount))
    print("+A: "+str(DerivedB.fCount))
    second = DerivedB("Mary")
    print("+B: "+str(DerivedA.fCount))
    print("+B: "+str(DerivedB.fCount))
    third = DerivedA("Carl")
    print("+A: "+str(DerivedA.fCount))
    print("+A: "+str(DerivedB.fCount))
    
    