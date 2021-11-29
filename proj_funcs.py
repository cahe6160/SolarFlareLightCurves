import numpy as np
import sys
from scipy.io import loadmat


def load_data(xRay_csvFile, flareIR_mFile, delimiter=','):
    """

    This function loads all xRay and EUV data for processing later.

    Parameters
    ----------
    xRay_csvFile : string
        This contains the pathname and filename for the xRay parameters.

    flareIR_mFile : string
        This contains the pathname and filename for the EUV data.

    delimiter : string
        The delimiter for the reading of the .csv file.

    Returns
    -------

    Xray_data_all : list
        Contains an array of all parameters in xRay.

    vec304 : list
        Contains all flux values in the 2010-2014 time range, in mW/m^2/s.

    vectime : list
        Times for SDO/EVE data, in MATLAB datetime.

    vecerr : list
        Error values for SDO/EVE data.

    wind_st : list
        Start times for the windows corresponding to different flares.

    sq304 : list
        Values for the solar quiet, in mW/m^2/s.

    """

    # Initialize Xray data.
    Xray_data_all = []

    # Error checking for the existence of the Xray file.
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

    # Exception handing for the existence of the EUV data file.
    try:
        dict304 = loadmat(flareIR_mFile, variable_names=['vec304', 'vectime',
                                                         'windowstart2',
                                                         'filt304', 'vecerr'])
    except FileNotFoundError:
        print(f"FileNotFoundError: File '{flareIR_mFile}' doesn't exist.",
              file=sys.stderr)
        sys.exit(1)

    # Storage of the variables contained in the EUV data file.
    vec304 = dict304['vec304']
    vectime = dict304['vectime']
    wind_st = dict304['windowstart2']
    sq304 = dict304['filt304']
    vecerr = dict304['vecerr']

    return Xray_data_all, vec304, vectime, vecerr, wind_st, sq304


def determine_flares(Xray_data_all, vec304, vectime, wind_st, sq304, vecerr,
                     ind):
    """

    This function determines the parameters and flux values corresponding to
    the specific flare identified by the index "num" and returns a smaller
    window for processing and parameter estimation.

    Parameters
    ----------
    Xray_data_all : list
        Contains an array of all parameters in xRay.

    vec304 : list
        Contains all flux values in the 2010-2014 time range, in mW/m^2/s.

    vectime : list
        Times for SDO/EVE data, in MATLAB datetime.

    wind_st : list
        Start times for the windows corresponding to different flares.

    sq304 : list
        Values for the solar quiet, in mW/m^2/s.

    vecerr : list
        Error values for SDO/EVE data.

    ind : int
        The index corresponding to the flare to process - iterated over
        externally.


    Returns
    -------
    irrev : list
        Contains all of the flux values corresponding to this specific flare.

    irrstd : float
        Contains the standard deviation (spread) of all irradiance values.

    sqev : list
        Contains all of the flux values correponding to solar quiet during
        event.

    errev : list
        Contains all of the error values corresponding to the event.

    sqstd : list
        Contains the standard deviation of the solar quiet values.

    timeev : list
        Contains the times corresponding to the particular event (ind).

    diff : list
        Contains the difference between the flux values and sq for this event.

    irrmean : float
        Contains the mean for a window surrouding the specific event, in a
        larger window than the event itself.
    """

    # Identify the start of the window within which the flare exists.
    wst = wind_st[0][ind]
    wst_n = wind_st[0][ind + 1]

    # Expand window - 2 hours before, start, 3.5 hours after start.
    wst_b = wst - (2/24)
    wend = wst + (3.5/24)

    # Find indices of times corresponding to window
    eventi = np.where(np.logical_and(vectime >= wst_b, vectime < wend))

    # Identify larger window
    eventi_larger = np.where(np.logical_and(vectime >= (wst_b-(12/24)),
                                            vectime < (wend + 12/24)))

    # Error check - do the nearby flares overlap? If so, go to next flare.
    if wend > wst_n:
        eventi = np.where(np.logical_and(vectime >= wst_b, vectime < wstn))

        # Identify larger window, as above.
        eventi_larger = np.where(np.logical_and(vectime >= (wst_b-(12/24),
                                                            vectime <
                                                            (wstn + 12/24))))

    # Define irradiance, time, solar quiet, error values for specific flare.
    irrev = vec304[eventi]
    timeev = vectime[eventi]
    sqev = sq304[eventi]
    errev = vecerr[eventi]

    # Larger window definition, just to get mean values.
    passage = vec304[eventi_larger]
    passaget = vectime[eventi_larger]

    # Process some statistics.
    sqstd = np.nanstd(sqev)
    irrstd = np.nanstd(irrev)
    irrmean = np.nanmean(passage)

    # Define difference array between light curve and solar quiet
    diff = irrev - sqev

    # Identify and remove large spikes in data (erroneous)
    for i in range(1, len(irrev)-1):
        if irrev[i] > (irrev[i-1]+2*irrstd) and irrev[i] > \
                (irrev[i+1]+2*irrstd):
            irrev[i] = np.nan

    for i in range(1, len(irrev)-1):
        if irrev[i] < (irrev[i - 1] - irrstd):
            irrev[i] = np.nan

    return irrev, irrstd, sqev, errev, sqstd, timeev, diff, irrmean


def find_flare_start_time(diff, irrev, irrstd, j, smSubTime, subTime, timeev,
                          tst, num, points):
    """

    This function iteratively and autonomously finds the start time of the
    flare in question (defined by "ind" above).

    Parameters
    ----------
    diff: list
        Contains the difference between the flux values and sq for this event.

    irrev : list
        Contains all of the flux values corresponding to this specific flare.

    irrstd : float
        The standard deviation of the flare flux values.

    j : int
        An initialized index, to iterate over to find the start time.

    subtime : int
        An initialized index - how many times do we count back from the found
        start time to where the flare actually starts?

    smSubTime : int
        The same as subTime, but a smaller version to avoid indexing error if
        necessary.

    timeev : list
        Contains the times corresponding to the particular event (num).

    tst : float
        Initialized value for the start time.

    num : int
        Initialized value for the number of points exceeding standard
        deviation.

    points : int
        The number of points which need to be found above the flare standard
        deviation consecutively before the start time is identified.

    Returns
    -------

    j : int
        The number of times iterated over before finding the start time.

    tst : float
        The start time of the flare, in MATLAB datetime.

    starti : int
        The index for the flare start time.
    """

    # Begin the iterative process to find flare start time.
    while j < len(irrev):

        # Add to "num" if the difference between irradiance and solar quiet
        # exceeds twice the standard deviation of irradiance for the flare.
        if diff[j] > (2 * irrstd):
            num += 1

            # If "num" has hit the criterion, identify the start time.
            if num == points and j > (num + points):
                tst = timeev[j - subTime]
                starti = j - subTime
                break
            elif num == points and j < (num + points):
                tst = timeev[j - smSubTime]
                starti = j - smSubTime
                break
            j += 1
        else:
            j += 1

    return j, starti, tst


def find_other_parameters(timeev, tst, irrev, sqev, irrstd, starti, endj,
                          diff):
    """

    This function iteratively and autonomously finds the peak and end times
    for the flare of interest.

    Parameters
    ----------
    timeev : list
        Contains the times corresponding to the particular event (num).

    tst : float
        The start time of the flare, in MATLAB datetime.

    irrev : list
        Contains all of the flux values corresponding to this specific flare.

    sqev : list
        Contains all of the flux values correponding to solar quiet during
        event.

    irrstd : float
        The standard deviation of the flare flux values.

    starti : int
        The index for the flare start time.

    endj : int
        Initialized value for the end index of the flare.

    diff: list
        Contains the difference between the flux values and sq for this event.

    Returns
    -------

    endj : int
        The found index of the flare end time.

    tend : float
        The found flare end time, MATLAB datetime.

    maxind : int
        The found index of the flare peak time.

    maxt : float
        The found flare peak time, MATLAB datetime.

    """

    # Initialize remaining parameters to be found.
    tend, maxt = 0, 0

    # If there is overlap, or an indexing issue, exit loop and print error.
    if starti > len(timeev) or endj > len(timeev):
        print('No classification of flare!')

    # Otherwise, if the start time is before the time array, exit loop and
    # print error.
    elif tst < timeev[1]:
        print('Start too early')

    # Otherwise, iterative over values to find the end time.
    else:

        # Initialize the end time as the start time, to begin.
        tend = timeev[starti]

        # Initialize a count.
        m = 0

        # Find the first time a large number of values have returned to solar
        # quiet
        for j in range(starti + 40, len(diff)):
            if diff[j] < 0.5*irrstd:
                m += 1
                if m == 50:
                    endj = j
                    print('endj', endj)
                    break

        # Account for if this point is NaN.
        if np.isnan(endj):
            print('endj is NaN')

        # If a time has been found, define the window of the flare between
        # start and end times.
        else:
            # The flare end time and flare window.
            tend = timeev[endj]
            window = irrev[starti:endj]

            # Initialize array of filtered values.
            avgmeans = np.zeros(len(window)-7)

            # Find the peak time by searching within a window of filtered
            # flux values.
            for k in range(8, len(window)):
                avgmeans[k - 7] = np.nanmean(window[k-7:k+7])

            # Find the peak index and the max time.
            maxind1 = np.where(avgmeans == max(avgmeans))
            maxind = starti + maxind1[0]
            maxt = timeev[maxind]

    return endj, tend, maxind, maxt


def store_times(endj, starti, timeev, tst, tend, maxt, startTimes, endTimes,
                maxTimes):
    """

    This function stores the indices and time parameters for the flare in
    question.

    Parameters
    ----------
    endj : int
        The index of the flare end time.

    starti : int
        The index of the flare start time.

    timeev : list
        The times corresponding to the flare window.

    tst : float
        The start time of the flare.

    tend : float
        The end time of the flare.

    maxt : float
        The peak time of the flare.

    startTimes: list
        List containing the known start times of flares.

    endTimes: list
        List containing the known end times of flares.

    maxTimes: list
        List containing the known peak times of flares.

    Returns
    -------

    startTimes: list
        List containing the known start times of flares, with the most recent
        appended.

    endTimes: list
        List containing the known end times of flares, with the most recent
        appended.

    maxTimes: list
        List containing the known peak times of flares, with the most recent
        appended.
    """

    # Append parameters for the most recent flare.
    startTimes.append(tst)
    endTimes.append(tend)
    maxTimes.append(maxt)

    return startTimes, endTimes, maxTimes
