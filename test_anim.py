#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 13:57:09 2017

@author: alexa
"""
from matplotlib import pyplot as PLT
from matplotlib import animation
import numpy as NP

import matplotlib.cm as cm
from matplotlib.colors import Normalize
#%%
fig, ax = PLT.subplots(1,1)
colors = list(range(E.count()))
colormap = cm.inferno
norm = Normalize()
norm.autoscale(colors)
Q = ax.quiver(0,0,0,0, angles='xy', scale_units='xy', scale=1,width=5e-3, color=colormap(norm(colors)))
PLT.xlabel('$S_z$')
PLT.ylabel('$S_x$')

#ax.set_xlim(-1.2, 1.2)
#ax.set_ylim(-1.2, 1.2)
ax.set_xlim(-2e-4, 2e-4)
ax.set_ylim(-2e-4, 2e-4)

def update_quiver(i, Q, E, pref, num = 1, diff=True):
    """updates the horizontal and vertical vector components by a
    fixed increment on each frame
    """
    
    if diff:
        Ss = [p['Ss'][num*i] - pref['Ss'][num*i] for p in E]
        Sx = [p['Sx'][num*i] - pref['Sx'][num*i] for p in E]
    else:
        Ss = [p['Ss'][num*i] for p in E]
        Sx = [p['Sx'][num*i] for p in E]

    Q.set_UVC(Ss,Sx)

    return Q,

# you need to set blit=False, or the first set of arrows never gets
# cleared on subsequent frames
anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q, E, E.getReference(), 10, True),
                               interval=1, blit=False)

PLT.show()
PLT.grid()

anim.save('./img/decoh.gif')

#%%

p=E[3]

fig, ax = PLT.subplots(1,1)
Q = ax.quiver(0,0,p['Ss'][0],p['Sx'][0],angles='xy', scale_units='xy', scale=1)

ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)

def update_quiver(i, Q, p):
    """updates the horizontal and vertical vector components by a
    fixed increment on each frame
    """

    Q.set_UVC(p['Ss'][i],p['Sx'][i])

    return Q,

# you need to set blit=False, or the first set of arrows never gets
# cleared on subsequent frames
anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q, p),
                               interval=1, blit=False)

PLT.show()
PLT.grid()
