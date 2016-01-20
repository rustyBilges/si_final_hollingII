## this script find the minimum number of timesteps for simulations using a given parameter set
## 12 column = T2P 
import numpy as np

parameters = np.genfromtxt('parameter_log.txt', delimiter=',')

dt = 0.0001  # "hi res"

T2P = parameters[:,12]

print 'T2P = ', np.min(T2P)

print 'for dt = 0.0001 -> ', np.min(T2P)/dt
