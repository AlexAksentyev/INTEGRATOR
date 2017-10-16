from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid, geom_point, geom_vline
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload
import numpy as NP
from matplotlib import pyplot as PLT

reload(ENT)
reload(PCL)

def form_state_list(Nx,Ny):
    xs = NP.linspace(-5e-3,5e-3,Nx)
    ys = NP.linspace(-5e-3,5e-3,Ny)
    
    
    StateList = list()
    for x in xs:
        for y in ys:
            StateList.append([x,y]+[0]*6+[1, 0, 0])
    
    return StateList

#%%

tLat = [ENT.MQuad(5e-2,-8.2,"QD"), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,7.36,"QF"),
        ENT.MQuad(5e-2,7.36,"QF"), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,-8.2,"QD"),
        ENT.MQuad(5e-2,-8.2,"QD"), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,7.36,"QF")]
    
#tLat = [ENT.MDipole(2,8,.46)]
#%%
StateList = form_state_list(3,3)
E = PCL.Ensemble.from_state(StateList)
E.track(tLat,10)
    
#%%

def pos(data):
    if data['Element'] == "QF": return 'Red'
    elif data['Element'] == "QD": return 'Blue'
    else: return 'Black'

df = E.getDataFrame()
df['Quad']=df.apply(pos, axis=1)
df = PDS.melt(df, id_vars=['PID','t','s', 'Turn','Element','Quad'])
dat = df.loc[df['variable'].isin(['Sx','Sy'])&df['PID'].isin([8])]
print(ggplot(dat,aes(x='s',y='value',color='variable')) + 
     geom_line() + geom_point(color=dat['Quad']) + theme_bw())
    

#%%
# fields 
QD = ENT.MQuad(5e-2,-8.2,"QD")
QF = ENT.MQuad(5e-2,7.36,"QF")

QList = [QD, QF]

StateList = form_state_list(5,5)

fldDF = PDS.DataFrame()
i=0
v = E[0].GammaBeta(E[0].fKinEn0)[1]
clight = E[0].CLIGHT()
for q in QList:
    for state in StateList:
        x,y = state[0:2]
        Px,Py,dW = state[3:6]
        Pc = E[0].Pc(dW+E[0].fKinEn0)
        Ps = NP.sqrt(Pc**2 - Px**2 - Py**2)
        vx = v*Px/Pc
        vy = v*Py/Pc
        vs = v*Ps/Pc
        Bx,By,Bs = q.BField(state)
        Fx = (vy*Bs-By*vs)*clight
        Fy = (vs*Bx-vx*Bs)*clight
        Fs = (vx*By-Bx*vy)*clight
        d = {'Quad':q.fName,'x':x,'y':y,'Bx':Bx,'By':By,'Bs':Bs,'Vx':vx,'Vy':vy,'Vs':vs,'Fx':Fx, 'Fy': Fy, 'Fs':Fs,'Px':Px,'Py':Py,'Pc':Pc}
        fldDF = fldDF.append(PDS.DataFrame(data=d,index=[i]))
        i += 1
        
#%%
names = NP.unique(fldDF['Quad'])
fldDF1 = fldDF[fldDF['Quad']==names[0]]
PLT.subplot(211)
PLT.quiver(fldDF1['x'], fldDF1['y'],fldDF1['Fx'],fldDF1['Fy'])
PLT.title(names[0])
fldDF1 = fldDF[fldDF['Quad']==names[1]]
PLT.subplot(212)
PLT.quiver(fldDF1['x'], fldDF1['y'],fldDF1['Fx'],fldDF1['Fy'])
PLT.title(names[1])
