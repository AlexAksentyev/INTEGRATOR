#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:53:58 2017

@author: alexa
"""


class Bundle(dict):
    """Bundle serves as an interface for easy access
    to a bundle of ensemble particle data.
    http://code.activestate.com/recipes/52308-the-simple-but-handy-collector-of-a-bunch-of-named/
    """
    def __init__(self, **kw):
        dict.__init__(self, kw)
        self.__dict__ = self

    def __deepcopy__(self, memo):
        res = self.__class__(**self.__dict__)
        return res

    def __repr__(self):
        return str(list(self.keys()))
    
from collections import Sequence

class MutableNamedTuple(Sequence): 
    """Abstract Base Class for objects as efficient as mutable
    namedtuples. 
    Subclass and define your named fields with __slots__.
    https://codereview.stackexchange.com/questions/173045/mutable-named-tuple-or-slotted-data-structure
    """
    __slots__ = ()
    def __init__(self, *args):
        for slot, arg in zip(self.__slots__, args):
            setattr(self, slot, arg)
    def __repr__(self):
        return type(self).__name__ + repr(tuple(self))
    # more direct __iter__ than Sequence's
    def __iter__(self): 
        for name in self.__slots__:
            yield getattr(self, name)
    # Sequence requires __getitem__ & __len__:
    def __getitem__(self, index):
        return getattr(self, self.__slots__[index])
    def __len__(self):
        return len(self.__slots__)
    