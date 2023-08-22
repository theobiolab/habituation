#Copyright © 2020 Ziyuan Zhao 
#Copyright © 2023 Lina Eckert, Rosa Martinez-Corral, Maria Sol Vidal <rosamcorral@hotmail.com>
#License GPLv3

import numpy as np

default_sliding_ht_threshold = 0.01 # should match int_threshold!

def pulse(t, Amin, Amax, T, ton):
    toff = T - ton
    if t >= 0 and t % T < ton:
        return Amax
    else:
        
        return Amin

def is_close(A, B, threshold=1e-8):
    if np.max(np.abs(A - B)) <= threshold:
        return True
    else:
        return False

def is_far(A, B, threshold=1e5):
    if np.max(np.abs(A - B)) >= threshold:
        return True
    else:
        return False


def sliding_maxnorm_ht(peaks_level, ht_threshold=default_sliding_ht_threshold):
    """Function to find habituation time on the basis of an array with the max level of the response at each output peak.""" 
    
    ### set general parameters #########################################
    # ht_threshold must not be lower than threshold of adaptive integration, otherwise the integration may have stopped prematurely
    
    ### find max peak index & normalize ###############################
    max_peak_index = np.argmax(peaks_level) # starts from 0
    
    # minmax normalization for cropped peak sequence
    if max_peak_index > 0: #select peaks including and after max peak
        peaks = peaks_level[max_peak_index::]
    else: #select all peaks
        peaks = peaks_level
    
    ### find habituated peak #########################################################
    # check peaks difference starting from tail for max efficiency
    # select peak which is less than ht_treshold lower than peak before
    i = len(peaks)-1

    # error: not enough peaks or peak difference too big
    if i <=0 or ht_threshold < (1 - peaks[i]/peaks[i-1]): 
        #print('last peaks have too high difference')
        return 0
    # error: difference of first peaks is already below treshold
    if (1 - peaks[1]/peaks[0] < ht_threshold): ###
        #print('threshold not met, 1st peaks too close')
        return 0
    
    ### calculate habituated cropped peak index
    # !!! monotonicity is a requirement for this algorithm to work !!!
    habituated_cropped_peak_index = -1 # initialize
    while i > 0:  
        habituated_cropped_peak_index = i # update selected peak
        if (1 - peaks[i]/peaks[i-1] > ht_threshold) or (peaks[i-1]>1e-5 and peaks[i]<1e-5 and peaks[0]>1e-3): #Lina please explain. Should these values perhaps be passed as keyword arguments with these as default?
            break
        i -=1
    
    # error: difference of first peaks is already below treshold
    if i == 0 :
        return 0
    
    # error: while loop was omitted
    if habituated_cropped_peak_index == -1:
        return 0
    else: 
    ### calculate habituation time #################################################
        # habituation time in number of peaks (not peak index!)
        ht = max_peak_index + habituated_cropped_peak_index+1
    
    return ht