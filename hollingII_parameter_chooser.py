
import numpy as np
import random as rnd
from hollingII_simulator import hollingII_simulator_parameter_tester
import sys, os

class parameter_chooser():
    
    def __init__(self):
        
        self.a_max = 100.0
        self.b_max = 100.0
        self.c_max = 100.0
        self.d_max = 100.0

        self.a_min = 0.100
        self.b_min = 0.100
        self.c_min = 0.100
        self.d_min = 0.100
        
        self.a = rnd.uniform(self.a_min, self.a_max)
        self.b = rnd.uniform(self.b_min, self.b_max)
        self.c = rnd.uniform(self.c_min, self.c_max)
        self.d = rnd.uniform(self.d_min, self.d_max)
        
        self.x0 = self.d/(self.c - 1)
        self.x1 = (self.a * self.d * self.c * (self.c - 1.0 - self.d)) / (self.b* (self.c - 1.0)**2 )
        
        self.tr = None
        self.det = None
        self.jacobian = np.zeros((2,2))
        self.eigen = None
        self.evec = None
        
        self.T2P = None
        
        
    def prnt(self):
        print("a = %f" %self.a)
        print("b = %f" %self.b)
        print("c = %f" %self.c)
        print("d = %f" %self.d)
        
    def check_tr(self):
        
        self.tr = self.a*(1.0-2.0*self.x0) - self.b*self.d*self.x1/((self.x0+self.d)**2) - 1.0 + self.c*self.x0/(self.x0+self.d)
        if self.tr <0:
            return True
        else:
            return False
        
    def check_det(self):
        
        self.jacobian[0,0] = self.a*(1.0-2.0*self.x0) - self.b*self.d*self.x1/((self.x0+self.d)**2)
        self.jacobian[0,1] = - self.b*self.x0/(self.x0+self.d)
        self.jacobian[1,0] = self.c*self.d*self.x1/((self.x0+self.d)**2)
        self.jacobian[1,1] = - 1.0 + self.c*self.x0/(self.x0+self.d)
        
        self.det = np.linalg.det(self.jacobian)
        if self.det > 0:
            return True
        else:
            return False
        
    def check_stable_spiral(self):
        
        if 4.0*self.det > self.tr**2:
            return True
        else:
            return False
        
    def check_eigenvalues(self):
        
        self.jacobian[0,0] = self.a*(1.0-2.0*self.x0) - self.b*self.d*self.x1/((self.x0+self.d)**2)
        self.jacobian[0,1] = - self.b*self.x0/(self.x0+self.d)
        self.jacobian[1,0] = self.c*self.d*self.x1/((self.x0+self.d)**2)
        self.jacobian[1,1] = - 1.0 + self.c*self.x0/(self.x0+self.d)
        
        self.eigen , self.evec = np.linalg.eig(self.jacobian)
        return self.eigen        
            
    def check_dynamics(self, outfile):
        ls = hollingII_simulator_parameter_tester(self.a, self.b, self.c, self.d, self.x0/2.0, self.x1/2.0, self.x0, self.x1)
        ls.run()
        
        if ls.params_ok:
            D = ls.get_dynamics()
            np.savetxt(outfile, D, delimiter=",")
            self.T2P = D[0,-1]
        
        return ls.params_ok
        
    def add_to_array(self,array, id):
        array[id,0] = id
        array[id,1] = self.a
        array[id,2] = self.b
        array[id,3] = self.c
        array[id,4] = self.d
        array[id,5] = self.x0
        array[id,6] = self.x1
        array[id,7] = self.tr
        array[id,8] = self.det
        
        self.check_eigenvalues()
        array[id,9] = np.real(self.eigen[0])
        array[id,10] = np.real(self.eigen[1])
        array[id,11] = np.imag(self.eigen[0])
        array[id,12] = np.imag(self.eigen[1])
        array[id,13] = self.T2P
        array[id,14] = self.x0/2.0
        array[id,15] = self.x1/2.0
        
        
        

if __name__=='__main__':
    
    n_params = 100
    param_id = 0
        
    outdir = "./"
    logfile = outdir + "parameter_log.txt"
    log_array = np.zeros((n_params, 16))
    

    while param_id < n_params:
        
        pc = parameter_chooser()
    
        if (pc.check_tr() and pc.check_det() and pc.check_stable_spiral()):
            
            dynamics_outfile = outdir + "dynamics/param_%d.dynamics" %param_id
            if pc.check_dynamics(dynamics_outfile):
                # save the good parameters in log file
                #pc.prnt()
                pc.add_to_array(log_array, param_id)
                print("pID = %d" %param_id)
                param_id += 1
    
    np.savetxt(logfile, log_array, delimiter=',')            
        
