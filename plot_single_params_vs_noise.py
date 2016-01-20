## produces plot of estimated parameter values; variance; and error functions versus noise intensity
## FOR A SINGLE PARAMETER SET:
## 	p0_versus_noise.txt

import numpy as np
import matplotlib.pyplot as plt

pID = 0 
nsamples = 1000 
parameters = np.genfromtxt("parameter_log.txt", delimiter=',')
#results = np.genfromtxt("single_params_vs_noise_pID_%d_nsamples_%d_dt_0.000100_reps_1000.results" %(pID,nsamples), delimiter=',')
results = np.genfromtxt("single_params_vs_noise_pID_%d_nsamples_%d_dt_0.000100_reps_2.HII.results" %(pID,nsamples), delimiter=',')

fsa = 15

fig, axes = plt.subplots(1,3, figsize=(15,6))
AX = axes.flatten()


AX[0].plot(results[:,0], results[:,1], 'ro', label='$\hat{r}_0$')
AX[0].axhline(parameters[0,1], color='r', linestyle='--')

AX[0].plot(results[:,0], results[:,2], 'bo', label='$\hat{J}_{00}$')
#AX[0].axhline(-parameters[pID,1], color='b', linestyle='--')
AX[0].plot(results[:,0], results[:,2+7], color='b', linestyle='--')

AX[0].plot(results[:,0], results[:,3], 'go', label='$\hat{J}_{01}$')
#AX[0].axhline(-parameters[pID,2], color='g', linestyle='--')
AX[0].plot(results[:,0], results[:,3+7], color='g', linestyle='--')

AX[0].plot(results[:,0], results[:,4], 'ko', label='$\hat{r}_{1}$')
AX[0].axhline(-1.0, color='k', linestyle='--')

AX[0].plot(results[:,0], results[:,5], 'mo', label='$\hat{J}_{10}$')
#AX[0].axhline(parameters[pID,3], color='m', linestyle='--')
AX[0].plot(results[:,0], results[:,5+6], color='m', linestyle='--')

AX[0].plot(results[:,0], results[:,6], 'co', label='$\hat{J}_{11}$')
#AX[0].axhline(0, color='c', linestyle='--')
AX[0].plot(results[:,0], results[:,6+6], color='c', linestyle='--')

AX[0].grid()
AX[0].set_xlabel('noise intensity ($\sigma_{noise}$)', fontsize = fsa)
AX[0].set_ylabel('mean estimated parameter value', fontsize = fsa)
AX[0].annotate("A", (0,0), (0.9,0.95), color='black', fontsize= fsa, fontweight='bold', xycoords='data', textcoords='axes fraction')

AX[1].plot(results[:,0], results[:,1+8], 'r-', label='$\hat{r}_0$')
AX[1].plot(results[:,0], results[:,2+8], 'b-', label='$\hat{J}_{00}$')
AX[1].plot(results[:,0], results[:,3+8], 'g-', label='$\hat{J}_{01}$')
AX[1].plot(results[:,0], results[:,4+8], 'k-', label='$\hat{r}_{1}$')
AX[1].plot(results[:,0], results[:,5+8], 'm-', label='$\hat{J}_{10}$')
AX[1].plot(results[:,0], results[:,6+8], 'c-', label='$\hat{J}_{11}$')
AX[1].legend(loc=2, fontsize=fsa)
AX[1].grid()
AX[1].set_xlabel('noise intensity ($\sigma_{noise}$)', fontsize = fsa)
AX[1].set_ylabel('variance in estimates', fontsize = fsa)
AX[1].annotate("B", (0,0), (0.9,0.95), color='black', fontsize= fsa, fontweight='bold', xycoords='data', textcoords='axes fraction')

AX[2].errorbar(results[:,0], results[:,7], yerr=np.sqrt(results[:,7+8]), fmt='ro', label='error 1')
AX[2].errorbar(results[:,0], results[:,8], yerr=np.sqrt(results[:,8+8]) ,fmt='bo', label='error 2')
AX[2].legend(loc=2, fontsize=fsa)
AX[2].grid()
AX[2].set_xlabel('noise intensity ($\sigma_{noise}$)', fontsize = fsa)
AX[2].set_ylabel('error function value', fontsize = fsa)
AX[2].annotate("C", (0,0), (0.9,0.95), color='black', fontsize= fsa, fontweight='bold', xycoords='data', textcoords='axes fraction')

plt.tight_layout()
#plt.savefig("plots/single_params_v_noise_pID_%d_nsamples_%d.png" %(pID, nsamples))
plt.show()
