import re
import pandas as PDS
import CElement as ENT
import PyDSTool as DST

class Particle:
    
    def __init__(self,Mass0 = 1876.5592, KinEn0 = 270.11275, G = -.142987):
        clight = 2.99792458e8
        q = 1.602176462e-19

        self.pardict = {'Mass0': Mass0, 'KinEn0': KinEn0, 'G':G,
                        'q':q, 'clight':clight,
                        'm0': q*1e6*Mass0/clight**2}
        
        self.fndict = {'KinEn':(['dK'],'KinEn0*(1+dK)'), 
                       'Lgamma':(['dK'],'KinEn(dK)/Mass0 + 1'),
                       'Lbeta':(['dK'],'sqrt(pow(Lgamma(dK),2)-1)/Lgamma(dK)'),
                       'Pc':(['dK'],'sqrt(pow(Mass0 + KinEn(dK),2) - pow(Mass0,2))')
                       }
        
        self.reuse = {'Pc(0)':'v0_P0c','Pc(dK)':'v0_Pc', 
                      'Lbeta(dK)':'v0_Lbeta', 'Lgamma(dK)':'v0_Lgamma', 'KinEn(dK)':'v0_KinEn',
                      'v0_Lbeta*clight':'v1_V',
                      'v0_P0c*px':'v1_Px', 'v0_P0c*py':'v1_Py', 
                      'q*1e6/clight*v1_Px' : 'v2_Px', 'q*1e6/clight*v1_Py': 'v2_Py',
                      'sqrt(pow(v0_Pc,2)- pow(v1_Px,2) - pow(v1_Py,2))':'v2_Ps',
                      'v1_V*v1_Px/v0_Pc':'v2_Vx', 'v1_V*v1_Py/v0_Pc':'v2_Vy',
                      'v1_V*v2_Ps/v0_Pc':'v3_Vs', 'q*1e6/clight*v2_Ps' : 'v3_Ps'
                    }
        
        self.defs = dict()
        ell = list(self.reuse.items())
        for key,value in ell:
            self.defs.update({re.sub('.*_','',value):key})
            
            
class Ensemble:
    
    def __init__(self, StateDict):
        self.fIniStateDict = {key:value for key,value in StateDict.items()}
        self.fCount = len(self.fIniStateDict)
        
    @classmethod
    def from_state(cls, StateList):
        names = ['X'+str(i) for i in range(len(StateList))]
        StateList = [dict(zip(ENT.Element.fArgList,e)) for e in StateList] # name state vars
        d = dict(zip(names,StateList)) # name initial conditions
        return cls(d)
        
    def __getitem__(self, index):
        return list(self.fIniStateDict.items())[index]
    
    def __repr__(self):
        return str(self.fIniStateDict)
        
    def getDataFrame(self):
            rval = PDS.DataFrame()
            for name, traj in self.fTrajectories.items():
                pts = traj.sample()
                pd = PDS.DataFrame(dict(zip(pts.coordnames, pts.coordarray)))
                pd['s'] = pts.indepvararray
                pd['PID'] = name
                rval=rval.append(pd)
                
            return rval
        
    def set(self, name, **pdict):
        self.fIniStateDict[name].update(pdict)
        
    def rename(self, oldname, newname):
        self.fIniStateDict[newname] = self.fIniStateDict.pop(oldname)
        
    def listNames(self):
        return list(self.fIniStateDict.keys())
        
        
        