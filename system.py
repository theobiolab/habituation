#Copyright © 2020 Ziyuan Zhao 
#Copyright © 2023 Lina Eckert, Rosa Martinez-Corral, Maria Sol Vidal <rosamcorral@hotmail.com>
#License GPLv3

import os
import numpy as np

from adaint import integrate, integrate2, default_h_max, default_int_threshold, default_steps_per_time
from itertools import product
from calc_util import sliding_maxnorm_ht



class System:
    """Class used to store the ODE system and parameters, run the calculations for habituation and store the associated data."""
       
    def __init__(self, f, parameter_set, X0,output_var_idx=-1,
                 steps_per_time=default_steps_per_time,hmax=0):
        self.f = f
        self.parameter_set = list(parameter_set) #should be only the rates and other ODE parameters
        self.X0 = X0 #initial condition for habituation protocol
        self.output_variable=output_var_idx

        self.steps_per_time=steps_per_time
        self.step_size=1/self.steps_per_time
        self.hmax=hmax

        # This is a dictionary for holding outputs from each call to computational methods,
        # before it was filled based on full_output=True, now always
        self.computational_data = {}
        

    def compute(self, T=None, Ton=None, Amin=0, Amax=None, verbose=0,
                ht_threshold=default_int_threshold,
                recovery_threshold=0.95,
                **kwargs):
        """Compute habituation and recovery times and trajectories. 
        - T: period.
        - Ton: duration of the stimulus at level Amax.
        - Amax: level of stimulus during the peak (variable S in f)
        - Amin: level of basal stimulus. Usually 0.
        - verbose: set to 1 for printing output along the way.
        - See documentation of integrate for what the rest of arguments mean."""

        
        if T is None:
            raise ValueError("Please specify period. Exiting...")
        else:
            self.T=T
        if Ton is None:
            ValueError("Please specify Ton. Exiting...")
        else:
            self.Ton=Ton
        if Amax is None:
            ValueError("Please specify Amax. Exiting...")
        else:
            self.Amax=Amax
        self.Amin=Amin


        temp_parameter_set=[T,Ton,Amin,Amax]+list(self.parameter_set)
        
        # find ss before integration
        #X0_ss = find_ss(self.f, self.X0, temp_parameter_set)
        #if sum(X0_ss) < 0: # find_ss returns vector with -1 for error and positive entries else
            #return zero_filter_array, 0, 0
        #else:
            #self.X0_ss = X0_ss
            
        ### compute habituation trajectory using function integrate() #################
        
        state, trajectory, peaks_time, peaks_level, troughs_time, troughs_level,tvec = integrate(self.f, self.X0, temp_parameter_set, \
            int_threshold=ht_threshold, hmax=self.hmax,steps_per_time=self.steps_per_time,
            output_var_idx=self.output_variable, **kwargs)
        

        if verbose > 0:
            print("peaks level are: ", peaks_level)
        
        
        ### store data in dictionary #####################################################
        ht = sliding_maxnorm_ht(peaks_level,ht_threshold=ht_threshold) #calculate habituation time

        self.computational_data = {'trajectory': trajectory, 
                                       'tvec': tvec,
                                       'parameter_set': temp_parameter_set, 
                                       'peaks_time': peaks_time,
                                       'peaks_level': peaks_level, 
                                       'troughs_time': troughs_time,
                                       'troughs_level': troughs_level,
                                       'habituation_time':ht*T,
                                       'habituation_time_step':int(ht*T*self.steps_per_time) - 1 #index of the trajectory vectors corresponding to the timepoint of habituation. Useful to integrate without stimulation for recovery
                                       }
            
        
            
        ### compute recovery trajectory and time #########################################
        if ht > 0:
            rt = self.find_recovery_time(recovery_threshold=recovery_threshold)
        else:
            rt = 0
        
        
        return ht, rt 

    def integrate_posthabituation_atAmin(self, tend=10):
        """integrate for tend time units after habituation time. Initial condition is the trajectory state at habituation time.
        - param_set: [ k1,...,kn] 
            Rates and other ODE parameters. In the same order they are defined in f.
        """
        allpars=[self.T,self.Ton,self.Amin,self.Amax]+self.parameter_set
        
        ### initialize X0
        try:
            trajectory=self.computational_data["trajectory"]
            habituate_time_step=self.computational_data["habituation_time_step"]
            habituation_time=self.computational_data["habituation_time"]
            X0 = trajectory[habituate_time_step, :]
            X0[X0 < 0] = 0  # prevent negative blowup   # update per 0813
        except:
            print(self.computational_data)
            raise ValueError("Wrong data in computational_data")
        
        
        
        # integrate2 with negative t will integrate with stimulus set to Amin
        tvec=np.arange(-tend,0,self.step_size)
        trajectory_ = integrate2(self.f, X0, tvec, tuple(allpars), hmax=self.hmax)

        tvec=tvec+tend+habituation_time #offset the tvec so that the first value is the habituation time
        return [tvec,trajectory_]


    def apply_singlestimulus(self,X0,T=None,Ton=None, Amin=None, Amax=None,tvec=None):

        if T is None:
            T=self.T
        if Ton is None:
            Ton=self.Ton
        if Amin is None:
            Amin=self.Amin
        if Amax is None:
            Amax=self.Amax
        if tvec is None:
            tvec=np.arange(0, T, self.step_size)
        allpars=[T,Ton,Amin,Amax]+self.parameter_set
        trajectory=integrate2(self.f, X0, tvec,tuple(allpars), hmax=self.hmax)
        return [tvec, trajectory]

            
    def find_recovery_time(self, recovery_threshold=0.95):
        
        # parameters for the binary search
        max_depth = 12
        max_duration = self.T * 2 ** max_depth #maximum time allowed to test for recovery

        # set reference peak
        peaks_level=self.computational_data["peaks_level"]
        first_peak = peaks_level[0] # set to first, but not necessarily highest peak

        #first let system relax for long time
        r_tvec, recovery_trajectory=self.integrate_posthabituation_atAmin(tend=max_duration)

        # Binary algorithm for determining recovery time  0 < r.t. < max_period
        # under the assumption that the inhibitory species is monotonic decreasing
        perturbational_tvec = np.arange(0, self.T, self.step_size)  # time for the stimulation
        #if full_output:
        self.computational_data['perturbational_tvec'] = perturbational_tvec

        dt = int(max_duration * self.steps_per_time / 2) #time at which perturbation is applied. Start at half the time for which the system has been integrated
        t = 0
        while dt > 0:
            X0 = recovery_trajectory[t + dt - 1, :]
            ptvec,perturbational_trajectory = self.apply_singlestimulus(X0,tvec=perturbational_tvec)
            # normalize peak
            post_recovery_peak = np.max(perturbational_trajectory[:, self.output_variable]) / first_peak 
            
            # update if not yet recovered
            if post_recovery_peak < recovery_threshold: #check recovery peak height
                t = t + dt # update t
                #if full_output:
                if not 'perturbational_trajectory' in self.computational_data:
                    self.computational_data['perturbational_trajectory'] = {}
                #self.computational_data['perturbational_trajectory'][habituation_time + t / steps_per_time] = perturbational_trajectory
                self.computational_data['perturbational_trajectory'][t / self.steps_per_time] = perturbational_trajectory #at time t, it gives the corresponding trajectory of the system corresponding to the one period with stimulation and off
            dt = int(dt / 2) # update dt. If already recovered, this will reduce the time at which perturbation is applied by half
            
        #if full_output:
        self.computational_data['recovery_trajectory'] = recovery_trajectory
        self.computational_data['recovery_end_time_step'] = self.computational_data["habituation_time_step"] + (t + 1)
            
        return ((t + 1)/self.steps_per_time)
