#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:53:58 2017

@author: alexa
"""

def phi(operation,*w):
    s = '('
    for e in w: 
        try: 
            e1 = float(e)
            e = str(e)
            if e1 < 0: e = '('+e+')'
        except ValueError:
            pass
        s += e+operation
    s = s[0:len(s)-1] + ')'
    return s