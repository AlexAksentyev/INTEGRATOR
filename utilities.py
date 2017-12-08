#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:53:58 2017

@author: alexa
"""


class Bundle(dict):
    """ Bundle serves as an interface for easy access
        to a bundle of ensemble particle data
    """
    def __init__(self, **kw):
        dict.__init__(self, kw)
        self.__dict__ = self

    def __deepcopy__(self, memo):
        res = self.__class__(**self.__dict__)
        return res

    def __repr__(self):
        return str(list(self.keys()))
    