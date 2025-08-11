import pandas as pd
import numpy as np
from datetime import datetime
import sys
import shutil
import argparse
import os
import warnings
import time
from casatasks import tclean, imfit, imstat, imhead

# import homemade functions 
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
from helper_functions import fit_point_source

prefix = "/Users/jimmylynch/Desktop/radio/observations/23A-241/23A-241.ASASSN-14ae.2023-09-17/EVLA_L"
image_paths = [f"{prefix}/1.25GHz/23A-241.ASASSN-14ae.2023-09-17.1.25GHz.2700px.image.tt0",
               f"{prefix}/1.5GHz/23A-241.ASASSN-14ae.2023-09-17.1.5GHz.2700px.image.tt0",

# optional toggles
combine_results = True        # If multiple image paths, will append fit results for each and write to one big csv file for easy copy/pasting
source_free_region = None     # Specify casa region file or string to use in calculating the near-source RMS. If NONE, default annulus used.
print_results = True          # if True, prints results as a dictionary in terminal
write_results = True          # if True, writes results as an excel sheet in imfit_files subdirectory
write_regions = True          # if True, writes casa region files for beam, source region of fit, fit result, and source-free region
override_sfr_request = True   # if True, bypasses the error that will be thrown if the source-free region has a suspected source

# run the script
results_dicts = []
for image_path in image_paths:
    results_dict = fit_point_source(image_path, 
                                    source_free_region=source_free_region, 
                                    print_results=print_results, 
                                    write_results=write_results,
                                    write_regions=write_regions,
                                    override_sfr_request=override_sfr_request)
    results_dicts.append(results_dict)

if combine_results:
    df_save = pd.DataFrame(results_dicts)
    df_save.to_csv(f"{prefix}/all_fit_results.csv")
    