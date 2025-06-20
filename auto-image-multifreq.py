from casatools import regionmanager
import numpy as np
import os
import pandas as pd
import shutil
import sys
import time
import yaml

from datetime import datetime

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

# read params from config.yaml
measurement_set = config["measurement_set"]
source_name = config["source_name"]
image_size = config["image_size"]
try_point_source = config["try_point_source"]
best_guess_flux = config["best_guess_flux [Jy/beam]"]
split = config["split"]
use_single_band = config["use_single_band"]
single_band = config["single_band"]

# Function to scrape a listfile for information needed for tclean ================================================================
# Inputs:
#     listfile (str): path/to/listfile.txt, which was automatically created in location of measurement set in main
#     source_name (str): user-specified name of source, like "ASASSN-14ae"
#     band (str): user-specified VLA band to image, like "C"
#     use_manual_spws (boolean): trigger in main to use own spectral windows and overwrite those in the band
#     manual_spws (str): manual spectral windows, either in a range using ~ or all comma separated: combinations not yet supported
# Returns:
#     field (str): index of source in listfile
#     cell_size (float): resolution [arcsec/pixel] of image, calculated by dividing the synthesized beamwidth by a factor of 4, as is recommended by NRAO
#     spw_range (str): spectral window range to image, given in format 'start_spw~stop_spw'
#     central_freq (float): central frequency of VLA band given
#     ra (str) and dec (str): coordinates of source
def scrape_listfile(listfile, source_name):

    # open listfile
    with open(listfile) as f:
        lines = f.read().splitlines()
    
    # open VLA configuration schedule and resolution tables
    df_resolution = pd.read_csv("vla-resolution.csv")
    df_resolution = df_resolution.loc[:, ~df_resolution.columns.str.contains('^Unnamed')]
    df_schedule = pd.read_csv("vla-configuration-schedule.csv")
    
    # loop through lines to find important lines
    for i, line in enumerate(lines):
        if "Fields" in line:
            field_line = line
            field_indx = i
    
        if "Spectral Windows" in line:
            spw_line = line
            spw_indx = i
    
        # determine on which line the observation datetimes are listed
        if "Observed from" in line:
            time_line = line
            time_indx = i
    
    # FIELD ===============================================
    nfields = int(field_line.split(" ")[-1][0])
    ls = [lines[field_indx+1+i] for i in range(nfields+1)]
    for l in ls:
        if source_name in l:
            field = l.split()[0]
            ra = l.split()[3]
            dec = l.split()[4]
    
    # CONFIGURATION =======================================
    # find start and finish time for observations
    t0 = time_line.split()[2]
    t1 = time_line.split()[4]
    
    # convert to datetime objects
    lf_date_format = "%d-%b-%Y/%H:%M:%S.%f"
    date0 = datetime.strptime(t0, lf_date_format)
    date1 = datetime.strptime(t1, lf_date_format)
    
    df_date_format = "%Y %b %d"
    # find the row in the schedule dataframe that encapsulates the observation
    for i, row in df_schedule.iterrows():
        start_epoch = datetime.strptime(row["observing_start"], df_date_format)
        end_epoch = datetime.strptime(row["observing_end"], df_date_format)
    
        if (start_epoch <= date0) and (date0 < end_epoch):
            configuration = row["configuration"]
    
    # CELL SIZE ==============================================
    #central_freq = (df_resolution[df_resolution["band"] == band]["central_freq"].values[0]).item()
    #synthesized_beamwidth = df_resolution[df_resolution["band"] == band][configuration].values[0].item()
    #cell_size = synthesized_beamwidth/4
    
    # SPECTRAL WINDOWS =================================================
    # determine how many spectral windows there are
    nspws = int(spw_line.split(' ')[3].split('(')[-1])
    ls = [lines[spw_indx+1+i] for i in range(nspws+1)]
    
    # get formatting right
    result = []
    for line in ls:
        row = list(filter(None, line.split(' ')))
        result.append(row)
    
    # save as dataframe
    cols = result[0]
    cols = cols[0:8]+["BBC-Num", "Corr1", "Corr2", "Corr3", "Corr4"]
    data = result[1:]
    df = pd.DataFrame(data, columns=cols)

    # get list of bands in listfile
    bands = list(set([df["Name"].iloc[i].split("#")[0] for i in range(df.shape[0])]))
    
    rows_list = []
    for i, band in enumerate(bands):
    
        # cell size
        #central_freq = (df_resolution[df_resolution["band"] == band]["central_freq"].values[0]).item()
        central_freq = (df_resolution[df_resolution["band"] == band]["central_freq"]).item()
        synthesized_beamwidth = (df_resolution[df_resolution["band"] == band][configuration]).item()
        cell_size = synthesized_beamwidth/4
    
        # get section of df for the band
        in_band = [band+"#" in b for b in list(df["Name"].values)]
        indxs = np.where(in_band)[0]
        df_band = df.iloc[indxs]
    
        # remove two cal spws from X-band
        if band == "EVLA_X":
            df_band = df_band.iloc[2:]
    
        # split into frequency bands
        nspws = df_band.shape[0]
        df_lower = df_band.iloc[0:int(nspws/2)]
        df_upper = df_band.iloc[int(nspws/2):df_band.shape[0]]
        df_all = df_band
    
        # lower
        freq_ghz = round(df_lower["CtrFreq(MHz)"].values.astype(float).mean()/1000, 2)
        spws = df_lower["SpwID"].values.astype(int)
        spw_range = f"{min(spws)}~{max(spws)}"
        rows_list.append({"band":band, "split":"lower", "freq [GHz]":freq_ghz, "spws":spw_range, "cell size [arcsec/pixel]":cell_size})
    
        # upper
        freq_ghz = round(df_upper["CtrFreq(MHz)"].values.astype(float).mean()/1000, 2)
        spws = df_upper["SpwID"].values.astype(int)
        spw_range = f"{min(spws)}~{max(spws)}"
        rows_list.append({"band":band, "split":"upper", "freq [GHz]":freq_ghz, "spws":spw_range, "cell size [arcsec/pixel]":cell_size})

        # all
        freq_ghz = round(df_all["CtrFreq(MHz)"].values.astype(float).mean()/1000, 2)
        spws = df_all["SpwID"].values.astype(int)
        spw_range = f"{min(spws)}~{max(spws)}"
        rows_list.append({"band":band, "split":"all", "freq [GHz]":freq_ghz, "spws":spw_range, "cell size [arcsec/pixel]":cell_size})

    df_store = pd.DataFrame(rows_list)#columns=["band", "split", "freq [GHz]", "spws", "cell size [arcsec/pixel]"])
    df_store = df_store.sort_values(by=["freq [GHz]"])
    df_store = df_store.reset_index(drop=True)

    return df_store, field, ra, dec

def get_theoretical_rms(listfile):
    return None

def fit_point_source(image_name, estimates_file, region_flux, region_rms, region_non_detection, imfit_logfile):

    near_source_imstat_results = imstat(imagename=image_name, region=region_rms)
    near_source_rms = near_source_imstat_results["rms"][0]
    
    # if imfit works, it found something to fit to, but not necessarily your source
    # it's good to check manually by opening the image in CARTA
    try:
        rms = near_source_rms
        print(region_flux)
        imfit_results = imfit(imagename=image_name, region=region_flux, rms=rms, logfile=imfit_logfile)
    
        flux = imfit_results["results"]["component0"]["flux"]["value"][0]
        flux_err = imfit_results["results"]["component0"]["flux"]["error"][0]
        strength = flux/rms
    
        # just in case imfit does work but there isn't a 3 sigma detection
        if strength >= 3:
            detection = True
            result = "imfit success"
        else:
            flux = 3*rms
            flux_err = 0
            detection = False
            result = "imfit success, <3 RMS"

    # if there is nothing there, imfit will fail. Reporting flux as 3x the rms from imstat
    except:
        print(f"Imfit failed: check results carefully!")
        imstat_results = imstat(imagename=image_name, region=region_non_detection)
        rms = imstat_results["rms"][0]
        flux = 3*rms
        flux_err = 0
        strength = 0
        detection = False
        result = "imfit fail: flux is 3*RMS"

    # returning values in mJy
    return round(flux*1000, 3), round(flux_err*1000, 3), round(rms*1000, 3), round(strength, 1), detection, result

directory = os.path.dirname(measurement_set)
ms_prefix = os.path.splitext(measurement_set)[0]
ms_name = os.path.splitext(measurement_set.split("/")[-1])[0]
image_directory = f"{directory}/images"

listfile = directory+"/listfile.txt"

# create listfile and scrape for tclean parameters
listobs(vis=measurement_set, listfile=listfile, overwrite=True)
print(f"Created listfile {listfile} \n")
df_store, field, ra, dec = scrape_listfile(listfile, source_name)

# trim down df_store to just what user wants
if split == "whole":
    df_store = df_store[df_store["split"] == "all"].reset_index(drop=True)
elif split == "halves":
    df_store = df_store[df_store["split"].isin(["upper", "lower"])].reset_index(drop=True)

if use_single_band:
    df_store = df_store[df_store["band"] == single_band].reset_index(drop=True)

# store flux results
if try_point_source:
    fluxes, flux_errs, rmses, strengths, detections, results, dynamic_ranges = [], [], [], [], [], [], []
    logfile_path = f"{ms_prefix}.{image_size}px.flux_measurements.csv"

print(df_store)

for i, row in df_store.iterrows():
    
    spws = row["spws"]
    freq = f"{row['freq [GHz]']}GHz"
    band = row["band"]
    cell_size = row["cell size [arcsec/pixel]"]

    # set up a folder to dump all the tclean outputs for this band
    band_directory = f"{directory}/{band}"
    if not os.path.exists(band_directory):
        os.makedirs(band_directory)
    
    image_name = f"{band_directory}/{ms_name}.{freq}.{image_size}px"

    # check if the image already exists to avoid overwriting
    files_in_directory = os.listdir(band_directory)
    full_paths = [band_directory+"/"+f for f in files_in_directory]
    
    # if the image exists, can still choose to open it up and fit a source there
    if image_name+".image.tt0" in full_paths:
        
        if try_point_source:
            print(f"Image exists and try_point_source=True: opening image to fit point source")
        else:
            print(f"Image exists and try_point_source=False: stopping to avoid overwrite")
            sys.exit(0)

    # if image does not exist, make it
    else:

        t0 = time.time()
        print(f"Begin imaging {image_name} \n")
        
        # run tclean (user can change or specify additional tclean parameters here)
        tclean(vis=measurement_set, 
               imagename=image_name, 
               field=field, 
               spw=spws, 
               specmode='mfs', 
               nterms=2, 
               deconvolver='mtmfs', 
               imsize=[image_size, image_size], 
               cell=[cell_size], 
               weighting='briggs', 
               threshold='0.0mJy', 
               niter=2000, 
               pblimit=-0.01, 
               savemodel='modelcolumn')
        
        print(f"Finished imaging {image_name} in {round((time.time()-t0)/60)}mins.")
    
    # fit to point source and, if fails, return upper limit as 3*RMS in region where source should be
    if try_point_source:

        # logfile
        imfit_logfile = f"{band_directory}/imfit_results.txt"
        estimates_file = f"{band_directory}/estimates.txt"

        # write estimates file
        # format is peak intensity, peak x-pixel value, peak y-pixel value, major axis, minor axis, position angle, fixed.
        midx, midy = image_size/2, image_size/2
        beam_size = 4*cell_size
        estimates = [f"{best_guess_flux}, {midx}, {midy}, {beam_size}pix, {beam_size}pix, {0}deg"]

        with open(estimates_file, "w") as f:
            for line in estimates:
                f.write(line + "\n")

        # global results
        imstat_results = imstat(imagename=image_name+".image.tt0")
        rms_image = imstat_results["rms"][0]
        max_image = imstat_results["max"][0]
        dynamic_range = round(max_image/rms_image, 1)
        dynamic_ranges.append(dynamic_range)

        # define region for detection fit: using 2.5 times the synthesized beamwidth centered on source
        radius_of_fit = 10*cell_size
        region_flux = f"circle[[{ra}, {dec}], {radius_of_fit}arcsec]"

        # define region for near-source RMS measurement: using annulus with 100 synthesized beams squared area
        inner_rad_annulus = 10*cell_size
        outer_rad_annulus = 22.5*cell_size
        region_rms = f"annulus[[{ra}, {dec}], [{inner_rad_annulus}arcsec, {outer_rad_annulus}arcsec]]"

        # define region for on-source RMS measurement if non-detection: using circle with outer annulus size
        region_non_detection = f"circle[[{ra}, {dec}], {outer_rad_annulus}arcsec]"
    
        # try the fit and print values
        flux, flux_err, rms, strength, detection, result = fit_point_source(image_name+".image.tt0", 
                                                                            estimates_file, 
                                                                            region_flux, 
                                                                            region_rms, 
                                                                            region_non_detection, 
                                                                            imfit_logfile)
        if detection:
            print(f"Detection at {freq}: {flux} ± {flux_err} mJy.")
            print(f"RMS: {rms} mJy/beam")
        else:
            print(f"Non-detection at {freq}: <{flux} mJy")

        # append values to 
        fluxes.append(flux)
        flux_errs.append(flux_err)
        rmses.append(rms)
        strengths.append(strength)
        detections.append(detection)
        results.append(result)

    # move image to images folder
    source = f"{image_name}.image.tt0"
    destination = image_directory 
    shutil.copytree(source, os.path.join(destination, os.path.basename(source)), dirs_exist_ok=True)

# write logfile
if try_point_source:
    df_store["Flux [mJy]"] = fluxes
    df_store["Flux_err [mJy]"] = flux_errs
    df_store["RMS [mJy]"] = rmses
    df_store["Flux/RMS [-]"] = strengths
    df_store["Detection"] = detections
    df_store["Fit result"] = results
    df_store["Dynamic range"] = dynamic_ranges
    
    df_store.to_csv(logfile_path)
        
    
    
    


