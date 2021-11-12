import numpy
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
        dict304 = loadmat(flareIR_mfile, variablenames=['vec304','vectime','windowstart2','filt304','errev'])
    except FileNotFoundError:
        print(f"FileNotFoundError: File '{flareIR_mfile}' doesn't exist.",
              file=sys.stderr)
        sys.exit(1)
    
    #error checking - if these don't exist?
    vec304 = dict304['vec304']
    vectime = dict304['vectime']
    wind_st = dict304['windowstart2']
    sq304 = dict304['filt304']
    err304 = dict304['errev']
    
    return Xray_data_all, vec304, vectime, wind_st, sq304, err304
    
def determine_flares(Xray_data_all, vec304, vectime, wind_st, st_in, end_in, ind)

    #takes input of window start time and window start time of next flare
    #run this once for each flare
    
    wst = wind_st[ind]
    wst_n = wind_st[ind + 1]
    
    #expand window - 2 hours before, start, 3.5 hours after start
    wst_b = wst - (2/24)
    
    wend = wst + (3.5/24)
    
    #find indices of times corresponding to window
    eventi = np.where(np.logical_and(vectime >= wst_b, vectime < wend))
    
    #identify larger window
    eventi_larger = np.where(np.logical_and(vectime >= (wst_b-(12/24), vectime < (wend + 12/24))))
    
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
    sqev = filt304[eventi]
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

    # LAST STEP IN THIS FUNCTION - CHECK NAN FUNCTIONALITY
    # I.E. DOES PYTHON USE NAN THE SAME WAY? DO WE HAVE TO ACCOUNT FOR THEM?
    # ALSO, CAN WE USE NUMPY? I THINK I CAN WRITE ALL FUNCTIONS IN BUILT-IN PYTHON.
    
    return irrev, irrstd, sqev, sqstd, timeev, diff, irrmean
            
    
def find_flare_start_time(diff, first, irrev, irrstd, j, lastj, starti, smSubTime, subTime, timeev, tst, num, second):
    """
    Parameters
    ----------
    first : boolean
        This bool is used to run only certain lines of code if it's the absolute first time this function has been called.
    irrev : list
    irrstd : int  
    j : int
    lastj : int
    starti : int
    smSubTime : int
    subTime : int
    timeev : list
    tst : int
    num : int
    second : boolean
        This bool is used to run certain lines of code for the second set of times this function is called.
    
    Returns
    -------
    
    """
    while j < len(irrev): 
        if diff[j] > (2 * irrstd):
            num += 1
            if first && num == 1:
                lastj = j
            if num == points and j > (num + points): # 40 -> points 
                if first:
                    lastj = j
                tst = timeev[j - subTime]
                starti = j - subTime 
                break
            elif num == points and j < (num + points):
                if first:
                    lastj = j
                tst = timeev[j - smSubTime]
                starti = j - smSubTime
                break
        elif second:
            n = 0
    return irrev, irrstd, j, starti, timeev, tst

def find_other_parameters(j, timeev, tst, irrev, sqev, irrstd, num):
    """
    Parameters
    ----------
    j : int
    timeev : list
    tst : int
    irrev : list
    sqev : list
    irrstd : int  
    num : int
    
    Returns
    -------
    
    """
    if starti > len(timeev) or endj > len(timeev):
        print('No classification of flare!')
    elif tst < timeev[1]:
        print('Start too early')
    else:
        tst = timeev[starti]
        m = 0
        j = starti + 40
        for j in range(len(diff)):
            if diff(j) < num:
                m += 1
                if m == 50:
                    endj = j
                    break
        if np.nanstd(endj):
            print('endj is NaN')
        else:
            tend = timeev[endj]
            window = irrev[starti:endj]
            windowt = timeev[starti:endj]
            sqwindow = sqev[starti:endj]
            sq = sqev[1:starti]
            avgmeans = NaN(len(window)-7,1)
            k = 8
            maxind1 = float("-inf")
            for k in range (len(window)-7):
                avgmeans[k - 7] = np.nanmean(window[k-7:k+7])
                maxind1 = max(maxind1, avgmeans[k - 7])
            # maxind1 = find(avgmeans == max(avgmeans))
            maxind = starti+maxind1
            maxt = timeev[maxind]

    return tst, endj, starti, timeev, tend, maxt, tend

def store_times(endj, starti, timeev, tst, tend, maxt):
    """
    Parameters
    ----------
    endj : int
    starti : int
    timeev : list
    tst : int
    tend : int
    maxt : list
    
    Returns
    -------
    
    """
    # Critical times for that flare have been found (well, assuming it was nicely behaved enough for the code)! 
    # Store the times in arrays initialized above
    if ~np.nanstd(endj) and starti<len(timeev) and endj < len(timeev) && starti<endj:
        starttimes.append(tst)
        endtimes.append(tend)
        maxtimes.append(maxt)
    return startTimes, endTimes, maxTimes
    