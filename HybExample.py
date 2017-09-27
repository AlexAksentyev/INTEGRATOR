#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 10:11:07 2017

@author: alexa
"""

import PyDSTool as DST
from matplotlib import pyplot as PLT

event_args={'name':'threshold',
            'eventtol':1e-3,
            'eventdelay':1e-4,
            'term':True}

thresh_ev = DST.makeZeroCrossEvent('V-threshval',
        0, event_args, varnames=['V'],parnames=['threshval'])

leak_args = DST.args(name='leak')
leak_args.pars = {'I': 1.3, 'gl': 0.1, 'vl': -67, 'threshval': -65}
leak_args.xdomain = {'excited': 0} #'V': [-100,50],
leak_args.varspecs = {'V': 'I - gl*(V-vl)',
                      'excited': '0'}
leak_args.algparams = {'init_step': 0.02}
leak_args.events = thresh_ev
leak_args.abseps = 1e-7

DS_leak = DST.embed(DST.Generator.Vode_ODEsystem(leak_args),name='leak',tdata=[0,200])

spike_args = DST.args(name='spike')
spike_args.tdomain = [0.0, 200]
spike_args.varspecs = {'V': 'a-t', 'excited': '1'}

spike_args.pars = {'splen': 0.75,'threshval':-90,'a':-55}
spike_args.xdomain = {'excited': 1} # 'V': [-97, 101], 
spike_args.events = thresh_ev

DS_spike = DST.embed(DST.Generator.ExplicitFnGen(spike_args),name='spike')

allnames=['leak','spike']
leak_MI=DST.intModelInterface(DS_leak)
leak_info=DST.makeModelInfoEntry(leak_MI,allnames,[('threshold','spike')])
spike_MI=DST.intModelInterface(DS_spike)
spike_info=DST.makeModelInfoEntry(spike_MI,allnames,[('threshold','leak')])

modelInfoDict= DST.makeModelInfo([leak_info, spike_info])
mod_args={'name':'IF_model', 'modelInfo':modelInfoDict}

IFmodel=DST.Model.HybridModel(mod_args)

IFmodel.compute(trajname='test',tdata=[0,600],ics={'V':-80,'excited':0})
pts = IFmodel.sample('test',dt=.05)

PLT.plot(pts['t'],pts['V'])
