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
               f"{prefix}/1.75GHz/23A-241.ASASSN-14ae.2023-09-17.1.75GHz.2700px.image.tt0"]

'''
image_paths = [f"{prefix}/EVLA_L/1.25GHz/23A-241.ASASSN-14ae.2023-07-23.1.25GHz.2700px.image.tt0",
               f"{prefix}/EVLA_L/1.5GHz/23A-241.ASASSN-14ae.2023-07-23.1.5GHz.2700px.image.tt0",
               f"{prefix}/EVLA_L/1.75GHz/23A-241.ASASSN-14ae.2023-07-23.1.75GHz.2700px.image.tt0",
               f"{prefix}/EVLA_S/2.5GHz/23A-241.ASASSN-14ae.2023-07-23.2.5GHz.2700px.image.tt0",
               f"{prefix}/EVLA_S/3.0GHz/23A-241.ASASSN-14ae.2023-07-23.3.0GHz.2700px.image.tt0",
               f"{prefix}/EVLA_S/3.5GHz/23A-241.ASASSN-14ae.2023-07-23.3.5GHz.2700px.image.tt0"]


image_paths = [f"{prefix}/5.0GHz/Target_ASASSN_14ae_EVLA_C_final.image.tt0",
               f"{prefix}/6.0GHz/Target_ASASSN_14ae_EVLA_C_final.image.tt0",
               f"{prefix}/7.0GHz/Target_ASASSN_14ae_EVLA_C_final.image.tt0",
               f"{prefix}/9.0GHz/Target_ASASSN_14ae_EVLA_X_final.image.tt0",
               f"{prefix}/10.0GHz/Target_ASASSN_14ae_EVLA_X_final.image.tt0",
               f"{prefix}/11.0GHz/Target_ASASSN_14ae_EVLA_X_final.image.tt0"]
               #f"{prefix}/13.55GHz/Target_ASASSN_14ae_EVLA_KU_final.image.tt0",
               #f"{prefix}/15.08GHz/Target_ASASSN_14ae_EVLA_KU_final.image.tt0",
               #f"{prefix}/16.62GHz/Target_ASASSN_14ae_EVLA_KU_final.image.tt0",
               #f"{prefix}/20.06GHz/Target_ASASSN_14ae_EVLA_K_final.image.tt0",
               #f"{prefix}/22.0GHz/Target_ASASSN_14ae_EVLA_K_final.image.tt0",
               #f"{prefix}/23.94GHz/Target_ASASSN_14ae_EVLA_K_final.image.tt0"]
'''

# optional toggles
combine_results = True
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
    