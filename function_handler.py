import proj_funcs
import matplotlib.pyplot as plt
import argparse
import configparser
import sys


def input_parser():
    """
    This function parses the filename for the UFO file, specifies the column
    number for the state abbreviation, and identifies the start and end year
    to filter by year and select the years desired to count events while also
    identifying events corresponding to hours of the day.

    Returns
    -------
    inputs : class 'argparse.Namespace'
        Parsed inputs from config file containing filenames.
    """

    # Initialize the parser.
    my_parser = argparse.ArgumentParser()

    # Identify the config file for the parser.
    my_parser.add_argument(
        '-c',
        '--config',
        type=str,
        action='store',
        help="This option allows you to specify a config file.",
        required=False
    )

    # Do not allow for cli and config together.
    if len(sys.argv) > 3:
        raise ValueError(f"User tried to use both config and input args. "
                         f"Expected <script> --config <config_filename>. "
                         f"Got {sys.argv=}")
        sys.exit()

    # Collect the parser inputs.
    inputs = my_parser.parse_args()

    return inputs


def main():
    """
    This function calls proj_funcs and run.sh.
    Grabs input from the shell file.
    Then, calls necessary functions to generate output.
    """
    args = input_parser()

    # Collect config file inputs
    config_filename = args.config
    config = configparser.ConfigParser()
    config.read(config_filename)
    xr = config['FILES']['fileReadCSV']
    evef = config['FILES']['mReadFile']

    # Establish X-ray and SDO/EVE 304 Angstrom files

    # Run function which loads data
    Xray_data_all, vec304, vectime, vecerr, wind_st, sq304 = \
        proj_funcs.load_data(xr, evef)

    # Establish the flare index (or indices) to iterate over
    ind = 0
    
    for ind in range(24):
        print('################', ind, '################')
        # For the flare index above, determine the corresponding light curve
        irrev, irrstd, sqev, errev, sqstd, timeev, diff, irrmean = \
            proj_funcs.determine_flares(Xray_data_all, vec304,
                                        vectime, wind_st, sq304, vecerr, ind)

        # Initialize number of points exceeding large value - adding to
        # this increases probability of finding the start of the flare
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
        # to identify actual start time. Reduce as necessary (to 30, 20, 10)
        # for less stringent detection criteria.
        subTime = 80

        # A smaller backpedal, if we're very close to the start of the window, to
        # avoid indexing issues. Reduce as necessary (to 29, 19, 9) for less
        # stringent criteria ('points').
        smSubTime = 39

        # Run the function to find the flare start time.
        j, starti, tst = proj_funcs.find_flare_start_time(diff, irrev, irrstd,
                                                          j, smSubTime,
                                                          subTime, timeev,
                                                          tst, num, points, 
                                                         ind)

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
                                                                      subTime,
                                                                      timeev,
                                                                      tst, num,
                                                                      points,
                                                                     ind)
                else:
                    break

        # Now we have found the start time, the index of which is provided by 'j'
        startj = j

        # Begin the process of finding the rest of the flare parameters - initilize
        # end index first.
        endj = 0

        # Run function to find end and peak times.
        endj, tend, maxind, maxt = \
            proj_funcs.find_other_parameters(timeev, tst, irrev, sqev, irrstd,
                                             starti, endj, diff, ind)
        
        if tst > 0 and tend > 0 and maxt > 0:
            # Print the start time to the flare.
            print('Start time: ', tst)

            # Print max time
            print('Max time: ',maxt)
            
            # Print end time
            print('End time: ',tend)

        if endj == 0 and tend == 0 and maxind == 0 and maxt == 0:
            continue

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
        plt.savefig('/home/jovyan/final_project/SolarFlareLightCurves/flare_plots/lctest'+str(ind)+'.png')
        print(" ")

        
if __name__ == "__main__":
    main()
