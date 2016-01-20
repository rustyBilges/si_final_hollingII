# THIS SCRIPT:
# Tests repeated simulation for a single parameter set. Does Timme inference, saves result. Times whole process.
## Varies: noise intensity

from datetime import datetime
from configure import NOISE_VALUE, NUMBER_OF_PARAMETERS, PARAMETER_FILE, DT, REPEATS_PER_PARAMETER_SET, NUMBER_OF_SAMPLES, NUMBER_OF_BINS
from hollingII_simulator import hollingII_simulator
from sampler import sampler
from timme_calculator import timme_calculator
import numpy as np
import os, sys

##################################
##configure:
PARAMETER_FILE = './parameter_log.txt' #'/home/rusty/Documents/phd files/write_along/second_year/functional_response_paper/final_code_version/linearFR/test/src/parameter_log.txt'

REPEATS_PER_PARAMETER_SET = 2 #1000
DT = 0.0001
NUMBER_OF_SAMPLES = 1000 #10000

p = 0 #87  ## pID
printout = True #False
plot_dynamics = False
##################################

start_sim = datetime.now()
parameters = np.genfromtxt(PARAMETER_FILE, delimiter=',')

noise = range(10)  # CHANGE! 

RESULTS = np.zeros((len(noise), 1 + (6 + 2 + 4)*2))  ## as below. *2 for mean and variance over repeats.
RESULTS[:,0] = noise
nid = 0

for ni in noise:
    print "noise = ", ni
    summary_array = np.zeros((REPEATS_PER_PARAMETER_SET, 6 + 2 + 4))  # to store results for each simulation, used to calculate mean and variance of estimates
								  ## a0 | J00 | J01 | a1 | J10 | J11 | err1 | err2 | meA00 | meA01 | meA10 | meA11 (mean IM elements over simulation)
    
    for r in range(REPEATS_PER_PARAMETER_SET):
        
        a = parameters[p,1]
        b = parameters[p,2]
        c = parameters[p,3]
        d = parameters[p,3]
        
        x00 = parameters[p,14]
        x10 = parameters[p,15]
        T2P = parameters[p,13]

        ls = hollingII_simulator(a, b, c, d, x00, x10, DT, T2P, ni, plot_dynamics)
        ls.run()
        D = ls.get_dynamics()
        
	E_prey = np.asarray(ls.ext_prey)
        E_pred = np.asarray(ls.ext_pred)

	## now do inference:
        S = sampler(NUMBER_OF_SAMPLES, D)
        S.sample()

        calc = timme_calculator(S, NUMBER_OF_BINS)
	results, err = calc.calculate()
	results = results[0:6]
	results = np.append(results, err[0])
	results = np.append(results, err[1])
	## evaluate mean IM elements over sim:
	x0 = D[1,:]  # prey time series
	x1 = D[2,:]  # prey time series
	L  = np.ones(len(x0))
	results = np.append(results, np.mean(-a*L + b*x1/((x0+d)**2)) )  ## equala to \alpha_{00}
	results = np.append(results, np.mean(-b/(x0+d) ))  ## equala to \alpha_{00}
	results = np.append(results, np.mean(c*d/((x0+d)**2)) )  ## equala to \alpha_{00}
	results = np.append(results, 0 ) 				 ## equala to \alpha_{11} !!
        summary_array[r,:] = results


	if printout==True:
		print("prey extinctions = %d" %(len(E_prey)))
		print("pred extinctions = %d" %(len(E_pred)))

		print(np.shape(D))

		print("\nRESULTS:\n")
		print(results)
   
    RESULTS[nid,1:] = np.append(np.mean(summary_array,0), np.var(summary_array,0))
    nid += 1

np.savetxt('single_params_vs_noise_pID_%d_nsamples_%d_dt_%f_reps_%d.HII.results' %(p,NUMBER_OF_SAMPLES, DT, REPEATS_PER_PARAMETER_SET), RESULTS, delimiter=',')
    
stop_sim = datetime.now()
elapsed_sim = stop_sim-start_sim
print 'time for simulation' , elapsed_sim

