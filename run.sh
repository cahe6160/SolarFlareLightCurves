#!/bin/bash
chmod +x run.sh

flaresumfile='/home/jovyan/final_project/SolarFlareLightCurves/flaresummaries.txt'

config_file='/home/jovyan/final_project/SolarFlareLightCurves/config.ini'

python function_handler.py --config $config_file > $flaresumfile 2> errors.err
    