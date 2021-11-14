import numpy as np
import sys
from scipy.io import loadmat

def load_data(xRay_csvFile, flareIR_mFile, delimiter=','):

    Xray_data_all = []
    
    try:
        file = open(xRay_csvFile, 'r')
    except FileNotFoundError:
        print(f"FileNotFoundError: File '{xRay_csv_File}' doesn't exist.",
              file=sys.stderr)
        sys.exit(1)
        
    # Append lines to the list containing all lines.
    for line in file:
        linestring = line
        Xray_data_all.append(linestring.strip().split(delimiter))
    file.close()
    
    try:
        dict304 = loadmat(flareIR_mFile, variable_names=['vec304','vectime','windowstart2','filt304','vecerr'])
    except FileNotFoundError:
        print(f"FileNotFoundError: File '{flareIR_mFile}' doesn't exist.",
              file=sys.stderr)
        sys.exit(1)
    
    #error checking - if these don't exist?
    vec304 = dict304['vec304']
    vectime = dict304['vectime']
    wind_st = dict304['windowstart2']
    sq304 = dict304['filt304']
    vecerr = dict304['vecerr']
    
    return Xray_data_all, vec304, vectime, vecerr, wind_st, sq304
    
def determine_flares(Xray_data_all, vec304, vectime, wind_st, sq304, vecerr, ind):

    #takes input of window start time and window start time of next flare
    #run this once for each flare
    
    wst = wind_st[0][ind]
    wst_n = wind_st[0][ind + 1]
    
    #expand window - 2 hours before, start, 3.5 hours after start
    wst_b = wst - (2/24)
    
    wend = wst + (3.5/24)
    
    #find indices of times corresponding to window
    eventi = np.where(np.logical_and(vectime >= wst_b, vectime < wend))
    
    #identify larger window
    eventi_larger = np.where(np.logical_and(vectime >= (wst_b-(12/24)), vectime < (wend + 12/24)))
    
    #error check - flares overlap? If so, go to next
    if wend > wst_n:
        eventi = np.where(np.logical_and(vectime >= wst_b, vectime < wstn))
    
        #identify larger window
        eventi_larger = np.where(np.logical_and(vectime >= (wst_b-(12/24), vectime < (wstn + 12/24))))
                                   
    # define events
    irrev = vec304[eventi]
    # need to find where nan?
    timeev = vectime[eventi]
    #define solar quiet
    sqev = sq304[eventi]
    #define error series
    errev = vecerr[eventi]
    
    #larger window definition, just to get mean values
    passage = vec304[eventi_larger]
    passaget = vectime[eventi_larger]
    
    #some statistics
    sqstd = np.nanstd(sqev)
    irrstd = np.nanstd(irrev)
    
    irrmean = np.nanmean(passage)
    
    #difference array between light curve and solar quiet
    diff = irrev - sqev
    
    #identify/remove large spikes
    for i in range(1,len(irrev)-1):
        if irrev[i] > (irrev[i-1]+2*irrstd) and irrev[i] > (irrev[i+1]+2*irrstd):
            irrev[i] = np.nan
            
    for i in range(1,len(irrev)-1):
        if irrev[i] < (irrev[i - 1] - irrstd):
            irrev[i] = np.nan
    
    return irrev, irrstd, sqev, errev, sqstd, timeev, diff, irrmean
            
    
                                       
                                       
    
    