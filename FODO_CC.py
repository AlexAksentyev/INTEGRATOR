#from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid, geom_point, geom_vline
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload
import numpy as NP
from matplotlib import pyplot as PLT

reload(ENT)
reload(PCL)

xs = NP.linspace(-5e-3,5e-3,3)
ys = NP.linspace(-5e-3,5e-3,3)


StateList = list()
for x in xs:
    for y in ys:
        StateList.append([x,y]+[0]*6+[1, 0, 0])

#%%

tLat = [ENT.MQuad(5e-2,-.82,"QD"), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,.736,"QF"),
        ENT.MQuad(5e-2,.736,"QF"), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,-.82,"QD"),
        ENT.MQuad(5e-2,-.82,"QD"), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,.736,"QF")]
#%%

E = PCL.Ensemble.from_state(StateList)
E.track(tLat,5)
    

def pos(data):
    if data['Element'] == "QF": return 'R'
    elif data['Element'] == "QD": return 'B'
    else: return 'Black'

df = E.getDataFrame()
#df['Quad']=df.apply(pos, axis=1)
#df = PDS.melt(df, id_vars=['PID','t','s', 'Turn','Element','Quad'])
#%%
dfs = df.loc[df['PID']==8]
PLT.plot(dfs['s'],dfs['x'],label='x')
PLT.plot(dfs['s'],dfs['y'],label='y')
#%%
dat = df.loc[df['variable'].isin(['x','y'])&df['PID'].isin([8])]
ggplot(dat,aes(x='s',y='value',color='variable')) + \
     geom_line() + geom_point(color=dat['Quad']) + theme_bw()
     
    
