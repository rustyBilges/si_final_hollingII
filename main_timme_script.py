## this script travels through all the simulation results for a certain NOISE_VALUE (defined in configure.py)
## performs timme gradient fitting calculations for each simualtion rep

## results saved in :  "paramI_sampleS_binB.timme_fit"
## where 
## I = param_id
## S = NUMBER_OF_SAMPLES taken from dynamics
## B = NUMBER_OF_BINS used for range sampling

## results format:
## each row corresponds to a single run (i.e. row_num = rep_id)
##
## row elements: 6xJhat_entries (a0, J00, J01, a1, J10, J11), bin_min_prey, bin_min_pred, bin_max_prey, bin_max_pred, bin_centre_prey, bin_centre_pred, bin_mean_prey, bin_mean_pred, bin_var_prey, bin_var_pred, num_points_in_bin
## concatenated B times (one for each bin)

## !! THIS SCRIPT SHOULD ALSO GENERATE SOME SUMMARY RESULTS !!

import glob, os, sys
import numpy as np
from sampler import sampler
from timme_calculator import timme_calculator
from configure import NOISE_VALUE, NUMBER_OF_PARAMETERS, REPEATS_PER_PARAMETER_SET, NUMBER_OF_SAMPLES, NUMBER_OF_BINS, RESULTS_FOLDER


#parameters = np.genfromtxt(PARAMETER_FILE, delimiter=',')
    
for p in range(NUMBER_OF_PARAMETERS):
    
    if not os.path.exists(RESULTS_FOLDER+"noise_%f/results_p%d" %(NOISE_VALUE,p)) :
        os.makedirs(RESULTS_FOLDER+"noise_%f/results_p%d" %(NOISE_VALUE,p))
    #os.chdir("./results_p%d" %p)
    save_to_dir = RESULTS_FOLDER+"noise_%f/results_p%d/" %(NOISE_VALUE,p)

    #summary_array = np.zeros((REPEATS_PER_PARAMETER_SET, 10))  
    
    ## do no use this because some repeats may be missing!
    #for r in range(REPEATS_PER_PARAMETER_SET):
    
    file_list = glob.glob(RESULTS_FOLDER+'noise_%f/param_%d/*.dynamics' %(NOISE_VALUE,p))
    
    results_array = np.zeros((len(file_list), 17 * NUMBER_OF_BINS))  
    
    f_id = 0
    for file in file_list:
        
        S = sampler(NUMBER_OF_SAMPLES, file)
        S.sample()

        calc = timme_calculator(S, NUMBER_OF_BINS)
        results_array[f_id,:] = calc.calculate()
        f_id += 1

    np.savetxt(save_to_dir + "param%d_sample%d_bin%d.timme_fit" %(p,NUMBER_OF_SAMPLES,NUMBER_OF_BINS), results_array, delimiter=',')
        
    #os.chdir("../")
    


