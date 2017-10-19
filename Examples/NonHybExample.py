#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 16:01:40 2017

@author: alexa
"""

import PyDSTool as dst
from PyDSTool import args
import numpy as np
from matplotlib import pyplot as plt
import cProfile
from time import clock

pars = {'eps': 1e-2, 'a': 0.5}
icdict = {'x': pars['a'],
          'y': pars['a'] - pars['a']**3/3}
xstr = '(y - (x*x*x/3 - x))/eps'
ystr = 'a - x'

event_x_a = dst.makeZeroCrossEvent('x-a', 0,
                            {'name': 'event_x_a',
                             'eventtol': 1e-6,
                             'term': False,
                             'active': True},
                    varnames=['x'], parnames=['a'],
                    targetlang='python')  # targetlang is redundant (defaults to python)

DSargs = args(name='vanderpol')  # struct-like data
DSargs.events = [event_x_a]
DSargs.pars = pars
DSargs.tdata = [0, 35]
DSargs.algparams = {'max_pts': 3000, 'init_step': 0.02, 'stiff': True}
DSargs.varspecs = {'x': xstr, 'y': ystr}
DSargs.xdomain = {'x': [-2.2, 2.5], 'y': [-2, 2]}
DSargs.fnspecs = {'Jacobian': (['t','x','y'],
                                """[[(1-x*x)/eps, 1/eps ],
                                    [    -1,        0   ]]""")}
DSargs.ics = icdict
vdp = dst.embed(dst.Vode_ODEsystem(DSargs),name="VDP")

#%%
start = clock()
vdp.compute(trajname='test',ics=icdict,tdata=[0,35])
print('Computed trajectory in %.3f seconds.\n' % (clock()-start))
