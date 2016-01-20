## this class does gradient fitting a la Shandylia and Timme
## with the addition of range sampling (default number of bins=1)

## NEED TO IMPROVE BINNING ALGORITHM?
## currently points are split evenly between bins, with remainder added to last bin

## NEED TO ADD PARAMETER COMPARISON TO THIS? - OR USE SEPARATE ANAYLSIS SCRIPTS

## method calculate() returns row array with elements :
## 6xJhat_entries (a0, J00, J01, a1, J10, J11), bin_min, bin_max, bin_centre, bin_mean, num_points_in_bin

## also interested in estimated biomass flows??

from sampler import sampler
import numpy as np

class timme_calculator():
    
    def __init__(self, sampler, n_bins=1):
    
        self.n_bins = n_bins
        self.sampler = sampler
        
        self.n_samples = self.sampler.n_samples
        
        self.interpolated = np.zeros((5,self.n_samples-1))
        
        self.boundaries = None
        self.binned_data = []
        
    def interpolate(self):
        
        for s in range(self.n_samples-1):
            
            self.interpolated[0,s] = ( self.sampler.sampled_dynamics[0,s+1] + self.sampler.sampled_dynamics[0,s] ) /2.0
            self.interpolated[1,s] = ( self.sampler.sampled_dynamics[1,s+1] + self.sampler.sampled_dynamics[1,s] ) /2.0
            self.interpolated[2,s] = ( self.sampler.sampled_dynamics[2,s+1] + self.sampler.sampled_dynamics[2,s] ) /2.0
            
            self.interpolated[3,s] = ( self.sampler.sampled_dynamics[1,s+1] - self.sampler.sampled_dynamics[1,s] ) / ( self.sampler.sampled_dynamics[0,s+1] - self.sampler.sampled_dynamics[0,s] )
            self.interpolated[4,s] = ( self.sampler.sampled_dynamics[2,s+1] - self.sampler.sampled_dynamics[2,s] ) / ( self.sampler.sampled_dynamics[0,s+1] - self.sampler.sampled_dynamics[0,s] )        
            
    def binning(self, test=None):
        #according to prey density
        if test==None:
            sorted = np.sort(self.interpolated[1,:])
            data = self.interpolated
        else:
            sorted = np.sort(test[1,:])
            data = test
            
        n_per_bin = int(np.floor((len(sorted))/self.n_bins))
        
        self.boundaries = np.zeros((1,self.n_bins+1))
        self.boundaries[0,0] = sorted[0]
        
        bin_id = 0
        sort_id = 0
        while bin_id<self.n_bins:
           # find the next boundary
            sort_id += n_per_bin
            bin_id += 1
            self.boundaries[0,bin_id] = sorted[sort_id - 1]
        
        self.boundaries[0,-1] = sorted[-1]
        #print(self.boundaries)
        
        # select data for each bin:
        binned = data[:,(data[1,:]>=self.boundaries[0,0]) * (data[1,:]<=self.boundaries[0,1])]
        self.binned_data.append(binned)
        
        bin_id = 1
        while bin_id<self.n_bins:
            binned = data[:,(data[1,:]>self.boundaries[0,bin_id]) * (data[1,:]<=self.boundaries[0,bin_id+1])]
            self.binned_data.append(binned)
            bin_id += 1
        
        #print(self.binned_data)
            
    
    
    
 
    def bin_stats(self):
        bin_stats = []
        for i in range(self.n_bins):
            stats = np.zeros((1,11))
            stats[0,0] = np.min(self.binned_data[i][1,:])
            stats[0,1] = np.min(self.binned_data[i][2,:])
            stats[0,2] = np.max(self.binned_data[i][1,:])
            stats[0,3] = np.max(self.binned_data[i][2,:])
            stats[0,4] = (stats[0,1]+stats[0,2])/2.0
            stats[0,5] = (stats[0,3]+stats[0,4])/2.0
            stats[0,6] = np.mean(self.binned_data[i][1,:])
            stats[0,7] = np.mean(self.binned_data[i][2,:])
            stats[0,8] = np.var(self.binned_data[i][1,:])
            stats[0,9] = np.var(self.binned_data[i][2,:])
            stats[0,10] = len(self.binned_data[i][2,:])
                                    
            bin_stats.append(stats)
        
        #2xbin_min, 2xbin_max, 2xbin_centre, 2xbin_mean, 2xbin_var, 2xest_biomass_flow, num_points_in_bin
        #print(bin_stats)
        return bin_stats
        
    def calculate(self):
        self.interpolate()
	#print(self.interpolated)
        self.binning()
        stats = self.bin_stats()
        #jHat = self.jHat(self.interpolated, self.n_samples - 1)
        #print(jHat)
        
        results = np.zeros(0)
        
        for b in range(self.n_bins):
            
            jHat = self.jHat(self.binned_data[b], stats[b][0,-1])
            results = np.append(np.append(results, jHat), stats[b])
        
        #print(results)
        return results
    
    def jHat(self, binned_data, M):
    # calculates the estimated interaction matrix
    # M is number of points in this bin
        X0 = binned_data[3,:]
        G0 = np.zeros((3, M))
        G0[0,:] = binned_data[1,:]
        G0[1,:] = binned_data[1,:]*binned_data[1,:]
        G0[2,:] = binned_data[1,:]*binned_data[2,:]
        
        X1 = binned_data[4,:]
        G1 = np.zeros((3, M))
        G1[0,:] = binned_data[2,:]
        G1[1,:] = binned_data[2,:]*binned_data[1,:]
        G1[2,:] = binned_data[2,:]*binned_data[2,:]
        
        J0 =np.dot( np.dot(X0, np.transpose(G0)), np.linalg.inv (np.dot(G0, np.transpose(G0))))
        J1 =np.dot( np.dot(X1, np.transpose(G1)), np.linalg.inv (np.dot(G1, np.transpose(G1))))
        
        J = np.asarray([J0,J1])
        #print("shape of J = ")
        #print(np.shape(J))
        return J
    
        #print(J0)
        #print(J1)
        
        
        
        
        
        
        
if __name__=='__main__':
    
    S = sampler(100, 'example_run.dynamics', 'example_prey.extinctions', 'example_pred.extinctions')
    S.sample()

    calc = timme_calculator(S, 2)
    calc.calculate()
    #calc.binning(np.asarray([[1,2,3,4,5,6,7,8, 9],[1, 7, 5, 3, 4, 6, 2, 8, 9], [1, 7, 5, 3, 4, 6, 2, 8, 9]]))  # testing
    #calc.bin_stats()
