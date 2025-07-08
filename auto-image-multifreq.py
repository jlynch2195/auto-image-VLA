from casatools import regionmanager
import numpy as np
import os
import pandas as pd
import shutil
import sys
import time
import yaml
from datetime import datetime

# import homemade functions 
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
from helper_functions import *

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

# read params from config.yaml
measurement_set = config["measurement_set"]
source_name = config["source_name"]
image_size = config["image_size"]
try_point_source = config["try_point_source"]
split = config["split"]
use_single_band = config["use_single_band"]
single_band = config["single_band"]

# optional imfit params
source_free_region = None #config["source_free_region"] 
print_results = config["print_results"]          
write_results = config["write_results"]          
write_regions = config["write_regions"]          
override_sfr_request = config["override_sfr_request"]
print(write_regions)

# where am I and what am I looking at
root_dir = os.path.dirname(measurement_set)                   # should be program/source/obs_date/
ms_name = os.path.splitext(measurement_set.split("/")[-1])[0] # should be program.source.obs_date.ms/
ms_prefix = os.path.splitext(measurement_set)[0]              # should be program.source.obs_date

# set up listfile and scrape it for information
listfile = f"{root_dir}/listfile.txt"
listobs(vis=measurement_set, listfile=listfile, overwrite=True)
print(f"Created listfile \n")
df_store, field = scrape_listfile(listfile, source_name, split, use_single_band, single_band)

# set up filesystem for image and flux report generation
set_up_filesystem(df_store, root_dir, try_point_source)

# create the images 
print(df_store)
image_paths = tclean_helper(df_store, measurement_set, ms_name, image_size, field, root_dir, try_point_source)

# do flux measurement
if try_point_source: 
    print("\n")
    logfile_path = f"{root_dir}/imfit_results/{ms_name}.{image_size}px.flux_measurements.csv"
    run_fit_point_source(image_paths, source_free_region, print_results, write_results, write_regions, override_sfr_request, 
                         logfile_path, df_store)

    
    
    
    


