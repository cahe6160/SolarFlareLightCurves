#!/bin/bash
chmod +x run.sh

config_file='/home/jovyan/SolarFlareLightCurves/config.ini'

python function_handler.py --config $config_file
    