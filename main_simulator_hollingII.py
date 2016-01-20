# THIS SCRIPT RUNS LINEAR FR SIMULATIONS USING PRE-SELECTED PARAMETERS
# Repeated simulation for each parameter set.
# STORES RESULTS IN DIRECTORIES : noise_ni/paramI/
# results: >dynamics
#          >logfile

from configure import NOISE_VALUE, NUMBER_OF_PARAMETERS, PARAMETER_FILE, DT, REPEATS_PER_PARAMETER_SET
from hollingII_simulator import hollingII_simulator
import numpy as np
import os, sys


#if not os.path.exists("./noise_%f" %NOISE_VALUE) :
#    os.makedirs("./noise_%f" %NOISE_VALUE)

#os.chdir("./noise_%f" %NOISE_VALUE)    
parameters = np.genfromtxt(PARAMETER_FILE, delimiter=',')
    
for p in range(NUMBER_OF_PARAMETERS):
    
    if not os.path.exists("./param_%d" %p) :
        os.makedirs("./param_%d" %p)
    os.chdir("./param_%d" %p)

    summary_array = np.zeros((REPEATS_PER_PARAMETER_SET, 10))  # to store and save summary results for each simulation
    
    for r in range(REPEATS_PER_PARAMETER_SET):
        
        a = parameters[p,1]
        b = parameters[p,2]
        c = parameters[p,3]
        d = parameters[p,4]
        
        x00 = parameters[p,14]
        x10 = parameters[p,15]
        T2P = parameters[p,13]
        
        ls = hollingII_simulator(a, b, c, d, x00, x10, DT, T2P, NOISE_VALUE, True)
        ls.run()
        D = ls.get_dynamics()
        E_prey = np.asarray(ls.ext_prey)
        E_pred = np.asarray(ls.ext_pred)
        np.savetxt("rep_%d.dynamics" %r, D, delimiter=',')
        np.savetxt("rep_%d_prey.extinctions" %r, E_prey, delimiter=',')
        np.savetxt("rep_%d_pred.extinctions" %r, E_pred, delimiter=',')
        
        summary_array[r,0] = np.min(D[1,:])
        summary_array[r,1] = np.max(D[1,:])
        summary_array[r,2] = np.mean(D[1,:])
        summary_array[r,3] = np.var(D[1,:])
        
        summary_array[r,4] = np.min(D[2,:])
        summary_array[r,5] = np.max(D[2,:])
        summary_array[r,6] = np.mean(D[2,:])
        summary_array[r,7] = np.var(D[2,:])
        
        summary_array[r,8] = len(ls.ext_prey)
        summary_array[r,9] = len(ls.ext_pred)
    
    np.savetxt("param_%d.summary" %p, summary_array, delimiter=',')
        
    os.chdir("../")
    
