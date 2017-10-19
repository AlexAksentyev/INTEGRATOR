import re

class Particle:
    
    def __init__(self,Mass0 = 1876.5592, KinEn0 = 270.11275, G = -.142987):
        q = 1.602176462e-19
        clight = 2.99792458e8
        self.pardict = {'Mass0':Mass0, 'KinEn0':KinEn0, 'G':G,
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
                      'sqrt(pow(v0_Pc,2)- pow(v1_Px,2) - pow(v1_Py,2))':'v2_Ps',
                      'v1_V*v1_Px/v0_Pc':'v2_Vx', 'v1_V*v1_Py/v0_Pc':'v2_Vy',
                      'v1_V*v2_Ps/v0_Pc':'v3_Vs'
                    }
        
        self.defs = dict()
        ell = list(self.reuse.items())
        for key,value in ell:
            self.defs.update({re.sub('.*_','',value):key})
            
            
#    def track(Lattice, ntimes, FWD=True):
        