import re
import pandas as PDS
import CElement as ENT
import copy

class DB:
    def __init__(self):
          self._id2name_map = {}
          self._name2id_map = {}
    
    def add(self, ind, name):
      self._name2id_map[name] = ind
      self._id2name_map[ind] = name
      
    def map(self, namelist, indlist):
        self._name2id_map.update(dict(zip(namelist,indlist)))
        self._id2name_map.update(dict(zip(indlist,namelist)))

    def by_name(self, name):
        ind = self._name2id_map[name]
        assert self._id2name_map[ind] == name
        return ind
    
    def by_index(self, ind):
        name = self._id2name_map[ind]
        assert self._name2id_map[name] == ind
        return name

class Particle:
    
    def __init__(self,Mass0 = 1876.5592, KinEn0 = 270.11275, G = -.142987):
        clight = 2.99792458e8
        q = 1.602176462e-19

        self.fPardict = {'Mass0': Mass0, 'KinEn0': KinEn0, 'G':G,
                        'q':q, 'clight':clight,
                        'm0': q*1e6*Mass0/clight**2}
        
        self.fFndict = {'KinEn':(['dK'],'KinEn0*(1+dK)'), 
                       'Lgamma':(['dK'],'KinEn(dK)/Mass0 + 1'),
                       'Lbeta':(['dK'],'sqrt(pow(Lgamma(dK),2)-1)/Lgamma(dK)'),
                       'Pc':(['dK'],'sqrt(pow(Mass0 + KinEn(dK),2) - pow(Mass0,2))')
                       }
        
        self.fReuse = {'Pc(0)':'v0_P0c','Pc(dK)':'v0_Pc', 
                      'Lbeta(dK)':'v0_Lbeta', 'Lgamma(dK)':'v0_Lgamma', 'KinEn(dK)':'v0_KinEn',
                      'v0_Lbeta*clight':'v1_V',
                      'v0_P0c*px':'v1_Px', 'v0_P0c*py':'v1_Py', 
                      'q*1e6/clight*v1_Px' : 'v2_Px', 'q*1e6/clight*v1_Py': 'v2_Py',
                      'sqrt(pow(v0_Pc,2)- pow(v1_Px,2) - pow(v1_Py,2))':'v2_Ps',
                      'v1_V*v1_Px/v0_Pc':'v2_Vx', 'v1_V*v1_Py/v0_Pc':'v2_Vy',
                      'v1_V*v2_Ps/v0_Pc':'v3_Vs', 'q*1e6/clight*v2_Ps' : 'v3_Ps'
                    }
        
        self.fDefs = dict()
        ell = list(self.fReuse.items())
        for key,value in ell:
            self.fDefs.update({re.sub('.*_','',value):key})
            
            
class Ensemble:
    
    def __init__(self, StateDict, Name = "X"):
        self.fName = Name
        self.fIniStateDict = {Name+key:value for key,value in StateDict.items()}
        self.fCount = len(self.fIniStateDict)
        self.__fDB = DB()
        self.__fDB.map(self.fIniStateDict.keys(),list(range(self.fCount)))
        
    @classmethod
    def from_state(cls, StateList, Name = "X"):
        names = [str(i) for i in range(len(StateList))]
        StateList = [dict(zip(ENT.Element.fArgList,e)) for e in StateList] # name state vars
        d = dict(zip(names,StateList)) # name initial conditions
        return cls(d, Name)
    
    def __add__(self, other):
        res = copy.deepcopy(self)
        for name,state in other.fIniStateDict.items():
            if name in res.fIniStateDict: name += "_1"
            res.fIniStateDict.update({name:state})
            res.__fDB.add(res.fCount-1,name)
            res.fCount += 1
        return res
        
    def __getitem__(self, pid):
        if type(pid) is int: name,index = self.__fDB.by_index(pid), pid
        else: name,index = pid,self.__fDB.by_name(pid)
            
        res = {key:value for key,value in self.fIniStateDict[name].items()}
        cols = list(res.keys())
        res.update({'Name':name})
        cols = ['Name']+cols
        return PDS.DataFrame(res,index=[index])[cols]
    
    def __repr__(self):
        pd = PDS.DataFrame(self.fIniStateDict).T
        return str(pd)
        
    def getDataFrame(self):
        rval = PDS.DataFrame()
        traj0 = self.fTrajectories[self.__fDB.by_index(0)]
        np = len(traj0.timePartitions)
        evt = traj0(list(range(np+1)),asmap=True)['t']
        for name, traj in self.fTrajectories.items():
            pts = traj.sample()
            pd = PDS.DataFrame(dict(zip(pts.coordnames, pts.coordarray)))
            pd['s'] = pts.indepvararray
            pd['PID'] = name
            rval=rval.append(pd)
           
        rval.fTransitions = evt
        return rval
        
    def set(self, name, **pdict):
        self.fIniStateDict[name].update(pdict)
        
    def rename(self, oldname, newname):
        self.fIniStateDict[newname] = self.fIniStateDict.pop(oldname)
        self.__fDB.map([newname],[self.__fDB.by_name(oldname)])
        
    def listNames(self):
        return list(self.fIniStateDict.keys())
        
        
        