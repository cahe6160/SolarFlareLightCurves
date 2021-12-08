* **CSCI 6118/4118 Final Project** *
* *Parametrizing Solar Flares* *
* *Science Lead: Cole Tamburri* *
* *Computational Lead: Caroline Hernandez* *

With this software, we develop an autonomous method for detection of solar flare light curves in the Solar Dynamics Observatory Extreme Ultraviolet Experiment (SDO/EVE) 304 Angstrom line.  Flares correspond to those identified by the RibbonDB database (Kazachenko et al. 2017) within the 2010-2014 time period.  Roughly 2048 flares are included in the SDO/EVE log, though the software here has been applied to a smaller subset for the purposes of development.  Using the methods described below, the start, peak, and end times for flare light curves from these data are computed and printed to a file. 

csci_6118_proj.m : The original framework upon which this process is based; this code, originally written in MATLAB, is converted into python and improved for performance and usability.  *This is a very large file containing four years of data, and it is necessary to install GitHub's large file service (LFS) if contributing to code and attempting to push to the GitHub repository.*

csci_proj_arr.mat : The SDO/EVE 304 Angstrom data.

ribbondb_v1.0.csv : The parameters corresponding to the Xray line in the RibbonDB database, used to guide the identification of parameters in the 304 Angstrom (EUV) line.

In proj_funcs.py :

- load_data() takes as input (1) a .csv file containing Xray flare parameters from RibbonDB and (2) a .mat file containing a time series of SDO/EVE 304 data as well as error values obtained from the instrument. The function returns lists for useful time series from which flares will be extracted.

- determine_flares() uses the parameters from the RibbonDB database and identifies a rough window within which a flare is most likely to exist.  This window, a subset of the overall flare time series, is isolated and output to a list.

- find_flare_start_time() takes the output of the isolated window and uses an iterative process to identify the beginning of a flare, based on a certain number of values exceeding the solar quiet (undisturbed sun) flux.

- find_other_parameters() identifies the maximum time (peak flux) and end time of the flare (when it returns to the solar quiet flux.

- store_times() stores the start, peak, and end times for flares identified by the above function into arrays.

In function_handler.py :

- input_parser() uses config.ini to parse inputs from run.sh which determine which flares will be studied within a particular set.

- main() collects the inputs from the config.ini file, loads the data, and for each flare: isolates the window where it is likely to occur, runs through the iterative process to find temporal parameters, and saves to an array.  Additionally, the identified light curves with temporal parameters are printed to a file in the flare_plots folder, and parameters for each flare are output to the parameters.csv file, which shows the parameters in datenum units, a general timing structure often used in solar physics and by MATLAB which can be converted easily to Universal Time or any other format.  Finally, the results are printed in an easy-to-read flaresummaries.txt file, which identifies flares for which the parameters could be identified, as well as those for which the algorithm failed (due to low SNR or other reasons).

run.sh identifies the output and config files and runs "function_handler.py," outputting error logs to "errors.err."

Along the way, non-fatal error checks are included for (1) gaps in data (returning a RuntimeWarning), (2) no identification of flare, (3) a disagreement between the isolated window and the algorithmically determined start time of the flare.  All of these are printed into the flaresummaries.txt file, except for case (1), when it is still usually possible to identify flare parameters.  The RuntimeWarning is printed to the errors.err file.

Fatal error checks such as FileNotFoundError are included as well.

For version v1.0, we include the display output (flaresummaries.txt) file, flare plots (flare_plots) and the science output (paramaters.csv) file for the first 24 flares in the sample, and see that the method works well for flare identification.  Future work will include application of this code to a larger sample of flares, for identification of the flare impulsiveness index and other interesting solar flare characteristics.  However, these 24 flares show the successful application and limits of this process, both identifying flare parameters for most and showing examples of when the algorithm doesn't work.  In these cases there is no useable flare signature in the 304 Angstrom EUV chromospheric line, though one may have existed in XRay).