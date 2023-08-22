#Copyright © 2020 Ziyuan Zhao 
#Copyright © 2023 Lina Eckert, Rosa Martinez-Corral, Maria Sol Vidal <rosamcorral@hotmail.com>
#License GPLv3

from scipy.integrate import odeint
import numpy as np


from calc_util import pulse, is_close, is_far

### set default integration steps and parameters. 
default_int_threshold = 0.01
default_h_max = 0
default_steps_per_time=100


def integrate(f, X0, args,
              output_var_idx=-1,
              int_threshold=default_int_threshold,
              min_output_level=1e-4,
              hmax=default_h_max, 
              steps_per_time=default_steps_per_time,
              num_periods_per_expansion_attempts=10,
              max_expansion_attempts = 3, 
              steady_counter_threshold = 4,
              increasing_counter_threshold = 10):
    """
    Function to integrate an habituating system until it has habituated. The system is defined as a set of ODEs and 
    stimulation is defined as a square-wave input signal with period T, 
    Ton time units at amplitude Amax and T-Ton time units at amplitude Amin. To ensure accuracy,
    the integration proceeds in a stepwise manner for each period, and for each of the two input levels.
    The system is initially integrated for a given amount of time (T*num_periods_per_expansion_attempts) and tested for habituation. If the system has not 
    habituated, it continues to be integrated up to max_expansion_attempts times. 
    
    - f: f(X,t, S, k1, k2, ..., kn) 
        The function that defines the ODE system to integrate, to be used with odeint. 
        The function takes as first argument a vector of state variables (X), the time (t), the input level (S), and the rest of rates and parameters of the ODEs k1,...,kn. 

    - X0: initial condition
    - args: [T, Ton, Amin, Amax, k1,...,kn] 
            The relevant parameters. The first four parameters must be period, Ton, Amin, Amax. Then the rest of the rates, in the same order they are defined in f.
    - output_var_idx: index that corresponds to the variable considered to be the output to test for habituation. If -1 it means the last one. Otherwise, specify index according to the definition of f.
    - int_threshold: when the relative difference between two consecutive peaks is less than this threshold for 4 times in a row, the integration stops. 
      The integration also stops when the absolute levels of the output variable is less than min_output_level for 4 times in a row, or a combination of the two conditions for 4 times in a row.
    - hmax: maximum integration time step passed to odeint. If 0: solver-determined.
    - steps_per_time: number of data points of the trajectory to be retrieved for each time unit of integration.
    - num_periods_per_expansion_attempts: number of periods over which the system is integrated each time (to avoid integrating excessively long times in non-habituating systems)
    - max_expansion_attempts: integration proceeds periodwise, for up to num_periods_per_expansion * max_expansion_attempts. 
    - steady_counter_threshold: number of peaks where the output is less than min_output_level (very close to 0) or not changing by more than int_threshold before integration stops
    - increasing_counter_threshold: number of peaks that are allowed to be increasing instead of decreasing before aborting integration. This allows to simulate systems that exhibit sensitization at the beginning.
   
    """

    
    num_periods=num_periods_per_expansion_attempts

    step_size = 1 / steps_per_time
    half_step_size = step_size / 2
    

   
    ### constants
    # cut period, Ton, Amin, Amax out of arguments 
    period, Ton, Amin, Amax = args[0:4]  # argument should be a list instead of tuple now
    arguments = np.insert(args[4:],0,0)  # parameters taken by f. The first is the level of the input signal. Set here to 0 and below is set to the appropriate level, either Amin or Amax
    
    integration_duration = int(period * steps_per_time)  # the duration has unit of a time step
    Ton_duration = int(Ton * steps_per_time) #number of time steps in ON
    Toff_duration= int((period-Ton)*steps_per_time) #number of time steps in OFF
    max_integration_time = period * num_periods  # The initial maximal duration of integration time
    
    # set index of output variable
    variable_count = len(X0)
    if int(output_var_idx)==-1:
        output_variable = variable_count - 1  # last variable of the model must be output!
    else:
        output_variable = output_var_idx
    
    
    ### intermediate variables
    integration_start_time = 0
    integration_time = int(max_integration_time * steps_per_time)  # also in the unit of steps
    expansion_counter = max_expansion_attempts
    steady_counter = 0
    increasing_counter = 0

    ### output initialization
    peaks_time = []
    peaks_level = []
    troughs_time = [0]  # count the starting value as a trough so peaks and troughs can be paired
    troughs_level = [0]
    trajectory = np.empty([integration_time, variable_count]) # set dimension of output trajecotry
    
    # set starting concentrations X0
    xout = np.expand_dims(X0, axis=0)

    ### expand integration time and integrate (loop)
    
    integration_attempts = 0
    while steady_counter < steady_counter_threshold: # while not habituated yet or output not too low
        if integration_start_time >= integration_time: #at the end of the integration, the end-time replaces the start-time, so at some point if not habituation this will be true
            if expansion_counter > 0: # expansions still allowed 
                trajectory = np.concatenate([trajectory, 
                                             np.empty([integration_time, variable_count])])
                expansion_counter = expansion_counter - 1
                integration_time = integration_time * 2
            else:
                # end integration after max number of expansions
                break

        # Taking min of two values prevent unnecessary integration
        # at the end of the time range that could not be fitted into 
        # the trajectory array; 
        integration_end_time = min(integration_start_time + integration_duration, integration_time) #integrate over one period or less if it does not fit 

        
        # note, np.arange have weird end behavior when using nonint
        # step size due to floating point error, see
        # https://stackoverflow.com/questions/17531961/numpy-arange-strage-behaviour
        # so nudge the end of the array to make it between steps
        
        tvec = np.arange(integration_start_time / steps_per_time,
                        integration_end_time / steps_per_time - half_step_size, step_size) #end time will be either end of period or a bit less if reaching the end of the time allowed to integrate
        
        # integrate ODE periodwise
        arguments[0] = Amax
        xout_1 = odeint(f, xout[-1], tvec[:Ton_duration+1], args=tuple(arguments), hmax=hmax)
        arguments[0] = Amin
        xout_2 = odeint(f, xout_1[-1], tvec[Ton_duration:], args=tuple(arguments), hmax=hmax)

        ### check whether integration failed and retry with lower hmax #######################
        # if the current integration results in NaN values,
        # then we will reattempt the integration at smaller hmax values
        # until we still fail at hmax=1e-6
        for var_hmax in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
            if np.isnan(np.sum(xout[-1])) or np.any(xout_1<0) or np.any(xout_2<0):
                arguments[0] = Amax
                xout_1 = odeint(f, xout[-1], tvec[:Ton_duration+1], args=tuple(arguments),
                                hmax=var_hmax)
                arguments[0] = Amin
                xout_2 = odeint(f, xout_1[-1], tvec[Ton_duration:], args=tuple(arguments),
                                hmax=var_hmax)
            else: 
                break
                
        ### output
        # concatenate output for Ton and Toff
        xout = np.concatenate((xout_1, xout_2[1:,:]), axis=0) # crop doubly calculated intersect
        
        # add output of current period to trajectory
        trajectory[integration_start_time:integration_end_time, :] = xout
        
        # add peaks_time and peaks_level
        if len(xout[:, output_variable]) > 1: # This prevents producing a fake peak when we haven't really integrated
            # index of max output value
            peaks_time.append(np.argmax(xout[:, output_variable]) + integration_start_time)
            # extract peak level from trajectory
            peaks_level.append(trajectory[peaks_time[-1], output_variable])

        # break if too many increasing peaks/irregularities or habituated
        if len(peaks_time) > 2:
            if peaks_level[-1] >= peaks_level[-2]:
                increasing_counter += 1
                if increasing_counter >= increasing_counter_threshold:
                    break
                    
            # stops integration if peaks do not decrease substantially for n times in a row
            if peaks_level[-2] < min_output_level: # prevent division by 0
                steady_counter +=1
            elif (1 - peaks_level[-1]/peaks_level[-2]) < int_threshold:   
                steady_counter +=1 #update counter
            else:
                steady_counter = 0 #reset counter to 0 if too big decrease

        
        # Update the start time for next round of integration, 
        # expand the trajectory array if necessary
        # should be close to performance 
        integration_start_time = integration_end_time
       
    # Integration ended normally, now remove the empty parts of the trajectory array
    trajectory = trajectory[0:integration_end_time, :]

    #fill trough time and level list
    for i in range(len(peaks_time) - 1):
        troughs_time.append(peaks_time[i] + np.argmin(trajectory[peaks_time[i]:peaks_time[i + 1],
                                                                 output_variable]))
    troughs_level = trajectory[troughs_time, output_variable]
    
    tvecfull=np.arange(trajectory.shape[0])*step_size  #time corresponding to each trajectory output point

    return 0, trajectory, np.array(peaks_time), np.array(peaks_level),\
           np.array(troughs_time), np.array(troughs_level),\
           tvecfull 


### integrate2 used for recovery time ##########################################################
def integrate2(f, X0, tvec, args, hmax=0):
    """Function used to calculate the recovery time and when we want to let the system relax and then apply one single pulse. 
    Use with negative t for no input, and with positive t for the pulse."""

    def f2(X, t, T, T_on, Amin, Amax, *varargs):
        A = pulse(t, Amin, Amax, T, T_on)
        return f(X, t, A, *varargs)
    return odeint(f2, X0, tvec, tuple(args), hmax=hmax)


###############################################################################################




def find_ss(f, X0_random, args, pre_ss_duration=2000, hmax = default_h_max):
    """Find the steady-state of the system so that the habituation protocol starts from there.
    This function has not been heavily tested so use with caution."""

    # initialize
    no_var = len(X0_random)
    pre_ss_tvec = np.arange(-pre_ss_duration, 0, step_size) # tvec for integration
    
    period, Ton, Amin, Amax = args[0:4]  # argument should be a list instead of tuple now
    arguments = np.insert(args[4:],Amin,Amin)  # parameters taken by f. The first is the level of the input signal. Set here to Amin because we want to integrate when the stimulus is at its basal level.

    
    # integrate
    xout = odeint(f, X0_random, pre_ss_tvec, args=tuple(args), hmax=hmax)
    
    # check whether integration failed and reintegrate with lower hmax if necessary
    for hmax in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
        #print(hmax)
        if np.any(xout < 0) or np.any(np.isnan(xout)): 
            xout = odeint(f, X0_random, pre_ss_tvec, args=tuple(args), hmax=hmax)
        else:
            break
            
    # check whether ss was reached in final integration
    if np.any(xout < 0) or np.any(np.isnan(xout)):
        X0_ss = np.repeat(-1, no_var)
    else:
        X0_ss = xout[-1,:]
    
    ### check whether values compatible with total concentrations!!!!!!
    return X0_ss
