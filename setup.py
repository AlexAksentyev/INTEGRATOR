# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:39:14 2017

@author: Аксентьев
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "helloworld",
    ext_modules = cythonize('helloworld.pyx'),  # accepts a glob pattern
)