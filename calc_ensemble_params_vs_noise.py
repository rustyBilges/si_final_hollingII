# THIS SCRIPT:
# Tests repeated simulation for the ensemble of parameter sets. Does Timme inference, saves result. Times whole process.
## Varies: noise intensity
## saves mean and variance in relative error (excepting J11 - raw error since this parameter is always zero)

from datetime import datetime
from configure import NOISE_VALUE, NUMBER_OF_PARAMETERS, PARAMETER_FILE, DT, REPEATS_PER_PARAMETER_SET, NUMBER_OF_SAMPLES, NUMBER_OF_BINS
from hollingII_simulator import hollingII_simulator
from sampler import sampler
from timme_calculator import timme_calculator
import numpy as np
import os, sys

##################################
##configure:
PARAMETER_FILE = './parameter_log.txt' 

REPEATS_PER_PARAMETER_SET = 1
DT = 0.0001
NUMBER_OF_SAMPLES = 10000

params = range(10) ## change!!
printout = False
plot_dynamics = False
##################################

start_sim = datetime.now()
parameters = np.genfromtxt(PARAMETER_FILE, delimiter=',')

noise = range(10)  ## change!! 

RESULTS = np.zeros((len(noise), 1 + (6 + 2)*2))  ## as below. *2 for mean and variance over repeats.
RESULTS[:,0] = noise
nid = 0

Nparams = len(params)

for ni in noise:
    print "noise = ", ni
    summary_array = np.zeros((Nparams * REPEATS_PER_PARAMETER_SET, 6 + 2))  # to store results for each simulation, used to calculate mean and variance of estimates
					               			    ## a0 | J00 | J01 | a1 | J10 | J11 | err1 | err2
									    ## relative error in all of the above, excepting J11, err1, err2		 
									    ## therefore extra +4 not required as in single param calculators
    sa_id = 0
    for pi in params:    
    	
	for r in range(REPEATS_PER_PARAMETER_SET):
        
		a = parameters[pi,1]
		b = parameters[pi,2]
		c = parameters[pi,3]
		d = parameters[pi,4]
		
		x00 = parameters[pi,14]
		x10 = parameters[pi,15]
		T2P = parameters[pi,13]

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

        	x0 = D[1,:]  # prey time series
	        x1 = D[2,:]  # prey time series
        	L  = np.ones(len(x0))

	        a00 = np.mean(-a*L + b*x1/((x0+d)**2))   ## equala to \alpha_{00}
        	a01 = np.mean(-b/(x0+d) )                ## equala to \alpha_{00}
	        a10 = np.mean(c*d/((x0+d)**2))           ## equala to \alpha_{00}
        	a11 = 0                                  ## equala to \alpha_{11} !!

		## calculate relative error:
		re = []
		re.append(np.abs((results[0]-a)/a))
		re.append(np.abs((results[1]-a00)/a00))
		re.append(np.abs((results[2]-a01)/a01))
		re.append(np.abs((results[3]+1)/1.0))
		re.append(np.abs((results[4]-a10)/a10))
		re.append(np.abs(results[5]))
		re.append(np.abs(results[6]))
		re.append(np.abs(results[7]))
		
		summary_array[sa_id,:] = re  ## saving relative errors
		sa_id += 1

		if printout==True:
			print("prey extinctions = %d" %(len(E_prey)))
			print("pred extinctions = %d" %(len(E_pred)))

			print(np.shape(D))

			print("\nRESULTS:\n")
			print(results)
	   
    RESULTS[nid,1:] = np.append(np.mean(summary_array,0), np.var(summary_array,0))
    nid += 1

np.savetxt('ensemble_params_vs_noise_nsamples_%d_dt_%f_reps_%d.HII.results' %(NUMBER_OF_SAMPLES, DT, REPEATS_PER_PARAMETER_SET), RESULTS, delimiter=',')
    
stop_sim = datetime.now()
elapsed_sim = stop_sim-start_sim
print 'time for simulation' , elapsed_sim

