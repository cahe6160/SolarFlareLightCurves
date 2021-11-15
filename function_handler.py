import proj_funcs as p

def main():
    """
    This function calls proj_funcs and run.sh.
    Grabs input from the shell file.
    Then, calls necessary functions to generate output.
    """
    Xray_data_all, vec304, vectime, wind_st, sq304, err304 = p.load_data(xRay_csvFile, flareIR_mFile, delimiter=',')
    irrev, irrstd, sqev, sqstd, timeev, diff, irrmean = p.determine_flares(Xray_data_all, vec304, vectime, wind_st, st_in, end_in, ind)
    
    ################ lines 174 -> 200 ################
    num = 0
    lastj = 0
    tst = 0
    
    points = 40 # points -> 40, 30, 20, 10
    subTime = 80 # subTime -> 80, 60, 40, 20
    smSubTime = 39 # smSubTime -> 39, 29, 19, 9
    
    irrev, irrstd, j, starti, timeev, tst = p.find_flare_start_time(diff, True, irrev, irrstd, j, lastj, starti, smSubTime, subTime, timeev, tst, num, False)
    ################ lines 205 -> ################
    while points > 10:
        j = 0
        if tst == 0 or starti > 1500:
            points -= 10
            subTime -= 20
            smSubTime -= 10
            if points == 20:
                num = 0# %3
            irrev, irrstd, j, starti, timeev, tst = p.find_flare_start_time(diff, False, irrev, irrstd, j, starti, smSubTime, subTime, timeev, tst, num, False)
            
    startj = j
    ################ lines 275 ################
    
    ################ lines 283 -> ################
    num = 0
    tst, endj, starti, timeev, tend, maxt, tend = p.find_other_parameters(j, timeev, tst, irrev, sqev, irrstd, num)
    ################ lines 352 ################
    tst = 0
    
    points = 40
    subTime = 80
    smSubTime = 39
    flag, irrev, irrstd, j, starti, timeev, tst = p.find_flare_start_time(diff, False, flag, irrev, irrstd, j, lastj, starti, smSubTime, subTime, timeev, tst, num, True)
    
    while points > 10: # 6
        if tst == 0:
            points -= 10
            subTime -= 20
            smSubTime -= 10
            j = endj
            flag, irrev, irrstd, j, starti, timeev, tst = p.find_flare_start_time(diff, False, flag, irrev, irrstd, j, lastj, starti, smSubTime, subTime, timeev, tst, num, True)
    
    if tst == 0:
        starti = 100000
    startj = j 
    ################ lines 450 ################
    
    ################ lines 466 -> ################
    num = 0.5 * irrstd
    tst, endj, starti, timeev, tend, maxt, tend = p.find_other_parameters(j, timeev, tst, irrev, sqev, irrstd, num)
    if len(maxt) > 1:
        maxt = maxt[1]
    ################ lines 537 ################
    
    ################ lines 541 -> ################
    startTimes, endTimes, maxTimes = p.store_times(endj, starti, timeev, tst, tend, maxt)
    ################ lines 550 ################
if __name__ == "__main__":
    main()
