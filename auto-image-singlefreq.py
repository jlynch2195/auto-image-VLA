'''
Jimmy Lynch
University of Oregon, Department of Physics

Script to create an image from a single-frequency VLA observation using CASA's tclean method

Running the script: 
    (base) cd place/where/this/code/lives
    (base) casa
    (CASA) <1>: execfile("auto-image-singlefreq.py")

User inputs are defined below in main section. Change these accordingly. 

The script proceeds in the following way:
    1. Creates and scrapes a listfile for information to tclean
        a. Finds the observation date and uses the VLA schedule to determine the array configuration
        b. Uses the user-specified band or manual_spws to calculate the resolution and central frequency from the VLA sensitivity table
        c. Finds the source's field index
        d. Finds the source's coordinates
    2. Creates an image using tclean, using the parameters above
    3. Optionally fits a point source in region near source coordinates using imfit and imstat
        a. Returns the flux for a detection or an upper limit for a non-detection


A couple of notes:
    1. There's undoubtedly a better way to scrape the listfile, as locating specific lines via keywords is subject to failure if VLA changes listfile structure
    2. For manual_spws, the central frequency is calculated as the mean of the center frequencies for the manual_spws. This may or may not be correct
    3. For manual_spws, the cell_size (resolution) is calculated using theta ~ lambda/D, with D being the maximum baseline of the configuration

'''
import numpy as np
import os
import pandas as pd
import shutil
import sys
import time
import yaml

from datetime import datetime

'''
Function to scrape a listfile for information needed for tclean ================================================================
Inputs:
    listfile (str): path/to/listfile.txt, which was automatically created in location of measurement set in main
    source_name (str): user-specified name of source, like "ASASSN-14ae"
    band (str): user-specified VLA band to image, like "C"
    use_manual_spws (boolean): trigger in main to use own spectral windows and overwrite those in the band
    manual_spws (str): manual spectral windows, either in a range using ~ or all comma separated: combinations not yet supported
Returns:
    field (str): index of source in listfile
    cell_size (float): resolution [arcsec/pixel] of image, calculated by dividing the synthesized beamwidth by a factor of 4, as is recommended by NRAO
    spw_range (str): spectral window range to image, given in format 'start_spw~stop_spw'
    central_freq (float): central frequency of VLA band given
    ra (str) and dec (str): coordinates of source
'''
def scrape_listfile(listfile, source_name, band, use_manual_spws, manual_spws):

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
    central_freq = (df_resolution[df_resolution["band"] == band]["central_freq"].values[0]).item()
    synthesized_beamwidth = df_resolution[df_resolution["band"] == band][configuration].values[0].item()
    cell_size = synthesized_beamwidth/4
    
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
    df_spws = pd.DataFrame(data, columns=cols)
    
    # get section of df_spws corresponding to the band
    in_band = ["EVLA_"+band+"#" in b for b in list(df_spws["Name"].values)]
    indxs = np.where(in_band)[0]
    df_band = df_spws.iloc[indxs]
    
    # get spectral window range for tclean
    spws = df_band["SpwID"].values.astype(int)
    spw_range = f"{min(spws)}~{max(spws)}"

    # (IF) MANUAL SPECTRAL WINDOWS ==================================================
    if use_manual_spws:

        # check formatting to be able to select appropriate rows from df_band
        if "~" in manual_spws:
            start_spw = int(manual_spws.split("~")[0])
            if "," not in manual_spws:
                end_spw = int(manual_spws.split("~")[1])
                spw_indices = np.arange(start_spw, end_spw+1)
                spw_indices = [str(i) for i in spw_indices]
    
            # currently not making the time to allow this string configuration
            else:
                print(f"Spectral windows with both ~ ranges and comma-separated values not yet supported: please list all with commas")
                sys.exit(0)
    
        # but will enable all comma-separated values, so can just list them all if needed
        elif "," in manual_spws:
            spw_indices = [spw.strip() for spw in manual_spws.split(',')]
    
        else:
            print(f"Invalid manual_spws format: {manual_spws}, check https://casaguides.nrao.edu/index.php/Selecting_Spectral_Windows_and_Channels")
            sys.exit(0)

        # select only the df_band rows in manual_spws to calculate central frequency
        spw_range = manual_spws
        df_manual_spws = df_band[df_band['SpwID'].isin(spw_indices)]
        central_freq = round(df_manual_spws["CtrFreq(MHz)"].values.astype(float).mean()/1000, 2)

        # need to update cell_size based on central_freq using the maximum baseline of the VLA
        B_max_arr = {"A": 36.4, "B": 11.1, "C": 3.4, "D": 1.03}
        B_max = B_max_arr[configuration]*1000

        # find cell_size based on theta ~ lambda/B_max
        cell_size = (3*10**8)/(central_freq*10**9*B_max)*206265

    return field, cell_size, spw_range, central_freq, ra, dec


'''
Function to fit an image created with tclean in a region specific ==============================================================
Inputs:
    image_name (str): image_name, including the .image.tt0 suffix
    region (str): CASA region to fit, default defined in main as a 2.5*synthesized_beamwidth arcsec circle centered at the source coordinates
Returns:
    flux (float): flux measured by imfit in the case of a detection, or 3*rms measured by imstat in the case of a non-detection 
    rms (float): root-mean square error measured by imstat (our current method of flux measurement)
    detection (boolean): True if there was a 3 sigma detection, False if not (turns False either by imfit fail or flux<3*rms)
'''
def fit_point_source(image_name, region):

    # if imfit works, it found something to fit to, but not necessarily your source
    # it's good to check manually by opening the image in CARTA
    try:
        imfit_results = imfit(imagename=image_name, region=region)
        imstat_results = imstat(imagename=image_name, region=region)
    
        flux = imfit_results["results"]["component0"]["flux"]["value"][0]
        rms = imstat_results["rms"][0]
        strength_of_detection_in_sigma = flux/rms
    
        # just in case imfit does work but there isn't a 3 sigma detection
        if strength_of_detection_in_sigma >= 3:
            detection = True
        else:
            flux = 3*rms
            detection = False

    # if there is nothing there, imfit will fail. Reporting flux as 3x the rms from imstat
    except:
        print(f"Imfit failed: check results carefully!")
        imstat_results = imstat(imagename=image_name, region=region)
        rms = imstat_results["rms"][0]
        flux = 3*rms
        flux_err = 0
        detection = False

    # returning values in mJy
    return round(flux*1000, 3), round(rms*1000, 3), detection


'''
Import user inputs from config file here: =======================================================================================================
    measurement_set (str): path/to/measurement_set.ms
    source_name (str): name of source as defined by VLA observation
    band (str): VLA band, like "C"
    image_size (int): default 256

Note: if you want to change tclean parameters, scroll down to that command

'''
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

measurement_set = config["measurement_set"]
source_name = config["source"]
band = config["band"]
image_size = config["image_size"]
use_manual_spws = config["use_manual_spws"]
manual_spws = config["manual_spws"]
try_point_source = config["try_point_source"]


'''
Code executes below this comment block; no need to change ======================================================================
'''
# create listfile in same directory as measurement set (compared to sticking it where you run this code)
directory = os.path.dirname(measurement_set)
listfile = directory+"/listfile.txt"

# create listfile and scrape for tclean parameters
listobs(vis=measurement_set, listfile=listfile, overwrite=True)
print(f"Created listfile {listfile} \n")
field, cell_size, spw_range, central_freq, ra, dec = scrape_listfile(listfile, source_name, band, use_manual_spws, manual_spws)

# define image name (up to user, just a useful default)
image_name = directory+f"/{source_name}_{central_freq}GHz_{image_size}px"

# check if the image already exists to avoid overwriting
files_in_directory = os.listdir(directory)
full_paths = [directory+"/"+f for f in files_in_directory]

# if the image exists, can still choose to open it up and fit a source there
if image_name+".image.tt0" in full_paths:
    
    if try_point_source:
        print(f"Image exists and try_point_source=True: opening image to fit point source")
    else:
        print(f"Image exists and try_point_source=False: stopping to avoid overwrite")
        sys.exit(0)

# if not, image using tclean
else:

    t0 = time.time()
    print(f"Begin imaging {image_name} \n")
    
    # run tclean (user can change or specify additional tclean parameters here)
    tclean(vis=measurement_set, 
           imagename=image_name, 
           field=field, 
           spw=spw_range, 
           specmode='mfs', 
           nterms=2, 
           reffreq=central_freq,
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

    # up to debate how big the circle to fit within should be: using 2.5 times the synthesized beamwidth centered on source
    radius_of_fit = 10*cell_size
    
    # define region for fit
    region = f"circle[[{ra}, {dec}], {radius_of_fit}arcsec]"
    print(f"Fitting point source in region {radius_of_fit}arcsec centered at {ra}, {dec} \n")

    # try the fit and print values
    flux, rms, detection = fit_point_source(image_name+".image.tt0", region)
    if detection:
        print(f"Detection at {central_freq}GHz: {flux} +/- {rms} mJy")
    else:
        print(f"Non-detection at {central_freq}GHz: <{flux} mJy")
    
    
    










