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
from helper_functions import fit_point_source, fit_point_source_basic

'''
prefix = "/Users/jimmylynch/Desktop/radio/observations/20A-492/ASASSN-14ae/final_files"
image_paths = [f"{prefix}/5.0GHz/Target_ASASSN_14ae_EVLA_C_final.image.tt0",
               f"{prefix}/6.0GHz/Target_ASASSN_14ae_EVLA_C_final.image.tt0",
               f"{prefix}/7.0GHz/Target_ASASSN_14ae_EVLA_C_final.image.tt0"]
'''
prefix = "/Users/jimmylynch/Desktop/radio/observations/23A-241/iPTF16fnl/EVLA_C/lower-5.0GHz"
image_paths = [f"{prefix}/iPTF16fnl_23A-241_5.0GHz.image.tt0"]

# fitting procedure
fitting_procedure = "basic"
print_results = True          # if True, prints results as a dictionary in terminal
write_results = False          # if True, writes results as an excel sheet in imfit_files subdirectory
combine_results = False        # If multiple image paths, will append fit results for each and write to one big csv file for easy copy/pasting

# ======================================================================
# doesn't matter if using basic procedure

# override the command to fit at the middle of the image
ra_pix = None
dec_pix = None

# optional toggles
source_free_region = None     # Specify casa region file or string to use in calculating the near-source RMS. If NONE, default annulus used.
write_regions = True          # if True, writes casa region files for beam, source region of fit, fit result, and source-free region
override_sfr_request = True   # if True, bypasses the error that will be thrown if the source-free region has a suspected source

# ======================================================================
# run the script
results_dicts = []
for image_path in image_paths:
    if fitting_procedure == "basic":

        results_dict = fit_point_source_basic(image_path, 
                                              print_results=print_results, 
                                              write_results=write_results)
        
    else:
        
        results_dict = fit_point_source(image_path, 
                                        source_free_region=source_free_region, 
                                        print_results=print_results, 
                                        write_results=write_results,
                                        write_regions=write_regions,
                                        override_sfr_request=override_sfr_request,
                                        ra_pix=ra_pix,
                                        dec_pix=dec_pix)
    results_dicts.append(results_dict)

if combine_results:
    df_save = pd.DataFrame(results_dicts)
    df_save.to_csv(f"{prefix}/all_fit_results.csv")
    