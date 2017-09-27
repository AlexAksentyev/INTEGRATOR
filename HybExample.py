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

leak_thresh_ev = DST.makeZeroCrossEvent('V-threshval',
        1, event_args, varnames=['V'],parnames=['threshval'])

leak_args = DST.args(name='leak')
leak_args.pars = {'I': 1.3, 'gl': 0.1, 'vl': -67, 'threshval': -65}
leak_args.xdomain = {'V': [-100,50], 'excited': 0}
leak_args.varspecs = {'V': 'I - gl*(V-vl)',
                      'excited': '0'}
leak_args.algparams = {'init_step': 0.02}
leak_args.events = leak_thresh_ev
leak_args.abseps = 1e-7

DS_leak = DST.embed(DST.Generator.Vode_ODEsystem(leak_args),name='leak',tdata=[0,200])

spike_args = DST.args(name='spike')
spike_args.tdomain = [0.0, 1.5]
spike_args.varspecs = {'V': 'if(t<splen,50,-95)', 'excited': '1'}
spike_args.pars = {'splen': 0.75}
spike_args.xdomain = {'V': [-97, 51], 'excited': 1}

DS_spike = DST.embed(DST.Generator.ExplicitFnGen(spike_args),name='spike')

allnames=['leak','spike']
leak_MI=DST.intModelInterface(DS_leak)
leak_info=DST.makeModelInfoEntry(leak_MI,allnames,[('threshold','spike')])
spike_MI=DST.intModelInterface(DS_spike)
spike_info=DST.makeModelInfoEntry(spike_MI,allnames,[('time','leak')])

modelInfoDict= DST.makeModelInfo([leak_info, spike_info])
mod_args={'name':'IF_model', 'modelInfo':modelInfoDict}

IFmodel=DST.Model.HybridModel(mod_args)

IFmodel.compute(trajname='test',tdata=[0,60],ics={'V':-80,'excited':0})
pts = IFmodel.sample('test',dt=.05)

PLT.plot(pts['t'],pts['excited'])
