## class to take in dynamics array
## and sample to timeseries, for use by inference calculators

## NEED TO ADD EXTINCTION HANDLING OT THIS CLASS, OR TO TIMME CALCULATE CLASS
## probably best post-interpolation
## i.e. extinction at t=2 => discard interpolate data point at i=1

import numpy as np
#import matplotlib.pyplot as plt

class sampler():
    
    def __init__(self, n_samples, dynamics, prey_extinction_file=None, pred_extinction_file=None, plot=False):
        
        self.n_samples = n_samples
        
        self.dynamics = dynamics
    
        #self.ext_prey = np.genfromtxt(prey_extinction_file, delimiter=',')
        #self.ext_pred = np.genfromtxt(pred_extinction_file, delimiter=',')
        
        self.sampled_dynamics = np.zeros((3, self.n_samples))
        
        self.plot = plot
        
    def sample(self):
        
        length = np.shape(self.dynamics)[1]
        skip = np.floor( length / self.n_samples)
        #print(skip)
        if skip==0:
            print("warning not enough data for this resolution of sampling")
            return False
       
	self.sampled_dynamics = self.dynamics[:,0::skip]

	#self.sampled_dynamics = self.dynamics[:,0:self.n_samples]

        
        #if self.plot==True:
        #    
        #if False:  #True:
	#    plt.subplot(1,2,1)
        #    plt.plot(self.dynamics[0,:], self.dynamics[1,:], 'g')
        #    plt.plot(self.dynamics[0,:], self.dynamics[2,:], 'r')
        #    plt.subplot(1,2,2)
        #    plt.plot(self.sampled_dynamics[0,:], self.sampled_dynamics[1,:], 'go-')
        #    plt.plot(self.sampled_dynamics[0,:], self.sampled_dynamics[2,:], 'ro-')
        #    plt.show()
            
if __name__=='__main__':
    
    S = sampler(100, 'example_run.dynamics', 'example_prey.extinctions', 'example_pred.extinctions', True)
    S.sample()
        
        
