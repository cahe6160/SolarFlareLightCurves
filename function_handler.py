import proj_funcs
import matplotlib.pyplot as plt

# Establish X-ray and SDO/EVE 304 Angstrom files

xr = '/home/jovyan/final_project/SolarFlareLightCurves/ribbondb_v1.0.csv'

evef = '/home/jovyan/final_project/SolarFlareLightCurves/csci_proj_arr.mat'

# Run function which loads data
Xray_data_all, vec304, vectime, vecerr, wind_st, sq304 = \
    proj_funcs.load_data(xr, evef)

# Establish the flare index (or indices) to iterate over
ind = 0

# For the flare index above, determine the corresponding light curve
irrev, irrstd, sqev, errev, sqstd, timeev, diff, irrmean = \
    proj_funcs.determine_flares(Xray_data_all, vec304,
                                vectime, wind_st, sq304, vecerr, ind)

# Initialize number of points exceeding large value - adding to this increases
# probability of finding the start of the flare
num = 0

# Initialize start time variable
tst = 0

# Initialize iterative variable, which will be changed with every loop
j = 0

# Initialize criteria for flare detection - how many successive points
# exceed the solar quiet before we find a flare? Reduce this value for less
# stringent criteria, if necessary, while looping (to 30, 20, 10).
points = 40

# Number of points to "backpedal" from the found number of exceeded points,
# to identify actual start time. Reduce as necessary (to 30, 20, 10) for less
# stringent detection criteria.
subTime = 80

# A smaller backpedal, if we're very close to the start of the window, to
# avoid indexing issues. Reduce as necessary (to 29, 19, 9) for less stringent
# criteria ('points').
smSubTime = 39

# Run the function to find the flare start time.
j, starti, tst = proj_funcs.find_flare_start_time(diff, irrev, irrstd,
                                                  j, smSubTime,
                                                  subTime, timeev,
                                                  tst, num, points)

# As a sanity check, print the number of exceeding points found, the start
# time, and the start index.
print('points', points)
print('tst', tst)
print('starti', starti)

# Lower the criteria to find a flare, if the start time has not been found
# or the start time index is very large (unlikely to have found a flare).
if tst == 0 or starti > 1500:
    while points > 10:
        if tst == 0 or starti > 1500:
            j = 0
            points -= 10
            subTime -= 20
            smSubTime -= 10
            if points == 20:
                num = 0

            # Run find_flare_start_time for less stringent criteria
            j, starti, tst = proj_funcs.find_flare_start_time(diff, irrev,
                                                              irrstd, j,
                                                              smSubTime,
                                                              subTime, timeev,
                                                              tst, num, points)

# Now we have found the start time, the index of which is provided by 'j'
startj = j

# Begin the process of finding the rest of the flare parameters - initilize
# end index first.
endj = 0

# Run function to find end and peak times.
endj, tend, maxind, maxt = \
    proj_funcs.find_other_parameters(timeev, tst, irrev, sqev, irrstd,
                                     starti, endj, diff)

# Plot the light curve with the parameters which have been found.
fig, ax = plt.subplots()
ax.scatter(timeev, irrev*1000, c='blue', s=2, label='Flux')
ax.axvline(tst, c='g', label='Start')
ax.axvline(tend, c='r', label='End')
ax.axvline(maxt, c='black', label='Peak')
ax.grid()
ax.legend()
ax.set_xlabel('Time [datetime]', fontsize=12)
ax.set_ylabel('Flux '+r'$[\mu W/m^2/s]$', fontsize=12)
ax.set_title('Detected Parameters, Flare '+str(ind), fontsize=20)
plt.savefig('lctest'+str(ind)+'.png')
