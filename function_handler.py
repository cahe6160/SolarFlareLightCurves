import proj_funcs

xr = '/home/jovyan/SolarFlareLightCurves/ribbondb_v1.0.csv'

evef = '/home/jovyan/SolarFlareLightCurves/csci_proj_arr.mat'

Xray_data_all, vec304, vectime, vecerr, wind_st, sq304 = proj_funcs.load_data(xr,evef)
ind = 0
irrev, irrstd, sqev, errev, sqstd, timeev, diff, irrmean = proj_funcs.determine_flares(Xray_data_all, vec304, vectime, wind_st, sq304, vecerr, ind)

################ lines 174 -> 200 ################
starti = 100000
num = 0
tst = 0
j = 0

points = 40 # points -> 40, 30, 20, 10
subTime = 80 # subTime -> 80, 60, 40, 20
smSubTime = 39 # smSubTime -> 39, 29, 19, 9

irrev, irrstd, j, starti, timeev, tst = proj_funcs.find_flare_start_time(diff, irrev, irrstd, j, starti, smSubTime, subTime, timeev, tst, num, points)

print('points', points)
print('tst', tst)
print('starti', starti)

############### line 205 -> ################
if tst == 0 or starti > 1500:
    while points > 10:
        if tst == 0 or starti > 1500:
            j = 0
            points -= 10
            subTime -= 20
            smSubTime -= 10
            if points == 20:
                num = 0# %3
            irrev, irrstd, j, starti, timeev, tst = proj_funcs.find_flare_start_time(diff, irrev, irrstd, j, starti, smSubTime, subTime, timeev, tst, num, points)

startj = j
################ lines 275 ################
################ lines 283 -> ################
num = 0
endj = 0
tst, endj, starti, timeev, tend, maxt, tend = proj_funcs.find_other_parameters(j, timeev, tst, irrev, sqev, irrstd, num, starti, endj, diff)
################ lines 352 ################
    