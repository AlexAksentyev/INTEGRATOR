from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload

theme_bw()

state = [1e-3, -1e-3, 0, 0, 0, 1e-4, 0, 0, 1, 0]
ds_25 = ENT.Drift(.25)
ds_15 = ENT.Drift(.15)
ds2_2 = ENT.Drift(2.2)
fquad = ENT.MQuad(5,.831)
dquad = ENT.MQuad(5,-.86)
fsext = ENT.MSext(5,1.023)
dsext = ENT.MSext(5,-1.34)

p = PCL.Particle(state)

#%%

FODO = [ds_25, fquad, ds_25, dquad, ds_25]

B0 = .46
R = ENT.MDipole.computeRadius(p,B0)
mdip = ENT.MDipole(1.8, R, (0,B0,0))
sol = ENT.Solenoid(1.8, .46)

V = ENT.Wien.computeVoltage(p,R,.05)
B1 = ENT.Wien.computeBStrength(p,R,.05)
wa = ENT.Wien(1.808,R,.05,V,B0)

#%%

p.track(FODO,10,FWD=True)
x = [p.fStateLog[i][0] for i in p.fStateLog]
y = [p.fStateLog[i][1] for i in p.fStateLog]
t = [p.fStateLog[i][2] for i in p.fStateLog]
dW = [p.fStateLog[i][5] for i in p.fStateLog]
Sx = [p.fStateLog[i][6] for i in p.fStateLog]
Sy = [p.fStateLog[i][7] for i in p.fStateLog]
Ss = [p.fStateLog[i][8] for i in p.fStateLog]
H = [p.fStateLog[i][9] for i in p.fStateLog]

df = PDS.DataFrame({'x':x,'y':y,'Sx':Sx,'Sy':Sy,'t':t,'H':H,'dW':dW})
df = PDS.melt(df, id_vars=['t','H'])
ggplot(df.loc[df['variable'].isin(['x','y','Sx','Sy'])],aes(x='t',y='value'))+\
    geom_line() + facet_wrap('variable',scales='free')
