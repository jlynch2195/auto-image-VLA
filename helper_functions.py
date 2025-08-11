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

# make proper directories
def set_up_filesystem(df_store, root_dir, try_point_source):

    # make root/image directory
    os.makedirs(f"{root_dir}/images", exist_ok=True)

    # make root/band/freq/ directories
    for i, row in df_store.iterrows():
        band = row["band"]
        freq = row["freq [GHz]"]
        subdir_name = f"{band}/{freq}GHz"
        
        full_path = os.path.join(root_dir, subdir_name)
        os.makedirs(full_path, exist_ok=True)

    # make root/imfit_results/ directory
    if try_point_source:
        os.makedirs(f"{root_dir}/imfit_results", exist_ok=True)

# make images
def tclean_helper(df_store, measurement_set, ms_name, image_size, field, root_dir, try_point_source):

    # save image paths
    image_paths = []
    
    # loop through each row and make an image for it
    for i, row in df_store.iterrows():

        # parse parameters for imaging from row
        spws = row["spws"]
        freq = f"{row['freq [GHz]']}GHz"
        band = row["band"]
        cell_size = row["cell size [arcsec/pixel]"]
    
        # set up a folder to dump all the tclean outputs for this band/freq combo
        tclean_out_directory = f"{root_dir}/{band}/{freq}"
        if not os.path.exists(tclean_out_directory):
            os.makedirs(tclean_out_directory)

        # make image name and set full path
        image_name = f"{ms_name}.{freq}.{image_size}px"
        image_path = f"{tclean_out_directory}/{image_name}"
    
        # check if the image already exists to avoid overwriting
        files_in_directory = os.listdir(tclean_out_directory)
        full_paths = [tclean_out_directory+"/"+f for f in files_in_directory]
        
        # if the image exists, can still choose to open it up and fit a source there
        if image_path+".image.tt0" in full_paths:
            
            if try_point_source:
                print(f"Image exists and try_point_source=True: will open image to fit point source after all imaging complete.")
            else:
                print(f"Image exists and try_point_source=False: stopping to avoid overwrite")
                sys.exit(0)
    
        # if image does not exist, make it
        else:
    
            t0 = time.time()
            print(f"Began imaging {image_name} \n")
            
            # run tclean (user can change or specify additional tclean parameters here)
            tclean(vis=measurement_set,  
                   imagename=image_path, 
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

        # copy image to images folder
        source = f"{image_path}.image.tt0"
        destination = f"{root_dir}/images" 
        shutil.copytree(source, os.path.join(destination, os.path.basename(source)), dirs_exist_ok=True)

        # append image path for output
        image_paths.append(image_path)

    return image_paths

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
def scrape_listfile(listfile, source_name, split, use_single_band, single_band):

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
    if use_single_band:
        bands = [single_band]
    
    rows_list = []
    for i, band in enumerate(bands):
    
        # cell size
        print(df_resolution[df_resolution["band"] == band])
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
        if band == "EVLA_L":
            freq_ghz = 1.25
        spws = df_lower["SpwID"].values.astype(int)
        spw_range = f"{min(spws)}~{max(spws)}"
        rows_list.append({"band":band, "split":"lower", "freq [GHz]":freq_ghz, "spws":spw_range, "cell size [arcsec/pixel]":cell_size})
    
        # upper
        freq_ghz = round(df_upper["CtrFreq(MHz)"].values.astype(float).mean()/1000, 2)
        if band == "EVLA_L":
            freq_ghz = 1.75
        spws = df_upper["SpwID"].values.astype(int)
        spw_range = f"{min(spws)}~{max(spws)}"
        rows_list.append({"band":band, "split":"upper", "freq [GHz]":freq_ghz, "spws":spw_range, "cell size [arcsec/pixel]":cell_size})

        # all
        freq_ghz = round(df_all["CtrFreq(MHz)"].values.astype(float).mean()/1000, 2)
        if band == "EVLA_L":
            freq_ghz = 1.5
        spws = df_all["SpwID"].values.astype(int)
        spw_range = f"{min(spws)}~{max(spws)}"
        rows_list.append({"band":band, "split":"all", "freq [GHz]":freq_ghz, "spws":spw_range, "cell size [arcsec/pixel]":cell_size})

    df_store = pd.DataFrame(rows_list)#columns=["band", "split", "freq [GHz]", "spws", "cell size [arcsec/pixel]"])
    df_store = df_store.sort_values(by=["freq [GHz]"])
    df_store = df_store.reset_index(drop=True)

    # trim down df_store to just what user wants
    if split == "whole":
        df_store = df_store[df_store["split"] == "all"].reset_index(drop=True)
    elif split == "halves":
        df_store = df_store[df_store["split"].isin(["upper", "lower"])].reset_index(drop=True)

    return df_store, field

# required inputs
'''
image_paths = ['/Users/jimmylynch/Desktop/radio/observations/25A-060/AT2018bsi/2025-07-01/images/25A-060.AT2018bsi.2025-07-01.5.0GHz.2700px.image.tt0', '/Users/jimmylynch/Desktop/radio/observations/25A-060/AT2018bsi/2025-07-01/images/25A-060.AT2018bsi.2025-07-01.7.0GHz.2700px.image.tt0']

# optional toggles
source_free_region = None     # Specify casa region file or string to use in calculating the near-source RMS. If NONE, default annulus used.
print_results = True          # if True, prints results as a dictionary in terminal
write_results = True          # if True, writes results as an excel sheet in imfit_files subdirectory
write_regions = True          # if True, writes casa region files for beam, source region of fit, fit result, and source-free region
override_sfr_request = False  # if True, bypasses the error that will be thrown if the source-free region has a suspected source
'''
def fit_point_source(image_path, source_free_region=None, print_results=True, write_results=True, write_regions=True, override_sfr_request=False):

    def write_region(region, region_name):
        with open(f"{region_name}.crtf", "w") as f:
            f.write("#CRTFv0\n")        
            f.write(region + "\n")

    def compare_fit_to_beam(bmaj, bmin, bpa, fmaj, fmin, fpa):
    
        barea = np.pi*bmaj*bmin
        farea = np.pi*fmaj*fmin
        dArea = farea/barea 
        dPhi = bpa-fpa
    
        return dArea, dPhi
    
    def get_region_area(image_path, region, cell_size):
    
        results = imstat(image_path, region=region)
        n_pixels = results["npts"][0]
    
        return n_pixels*cell_size**2

    # where am I
    directory = os.path.dirname(image_path)
    image_name = os.path.basename(image_path)

    # make directory 
    subdir_name = 'imfit_files'
    subdir_path = os.path.join(directory, subdir_name)
    os.makedirs(subdir_path, exist_ok=True)

    # scrape information from header
    header = imhead(image_path, mode="list")
    x0, y0 = header["crpix1"], header["crpix2"] # pixel number of center of field
    cell_size = header["cdelt2"] # arcseconds per pixel
    source = header["object"]
    freq = round(header["restfreq"][0]/10**9, 2)
    dt = datetime.strptime(header["date-obs"], "%Y/%m/%d/%H:%M:%S.%f")
    obs_date = dt.strftime("%Y-%m-%d")
    telescope = header["telescope"]
    source = header["object"]
    
    # beam info
    bmaj = header["beammajor"]["value"] # arcsec
    bmin = header["beamminor"]["value"] # arcsec
    bpa = header["beampa"]["value"] # degrees
    beam_region = f"ellipse[[{x0}pix, {y0}pix], [{bmaj/2}arcsec, {bmin/2}arcsec], {bpa}deg]"
    if write_regions: write_region(region=beam_region, region_name=f"{subdir_path}/{image_name}.beam_region")
    beam_radius = 0.5*(bmaj + bmin)

    # get global image stats
    imstat_results_image = imstat(image_path)
    rms_image = imstat_results_image["rms"][0]
    max_image = imstat_results_image["max"][0]
    dynamic_range = round(max_image/rms_image, 1)

    # define region for estimating near-source RMS
    if source_free_region is None:

        # try an annulus at the point source with area 50 beam areas
        inner_rad_annulus = 2.5*beam_radius
        outer_rad_annulus = 7.5*beam_radius
        source_free_region = f"annulus[[{x0}pix, {y0}pix], [{inner_rad_annulus}arcsec, {outer_rad_annulus}arcsec]]"

        # check to make sure the region is source free
        # if the max value in the region is greater than 5x the image SNR, there's possibly a source there
        region_SNR = imstat(image_path, region=source_free_region)["max"][0]/rms_image
        if region_SNR >= 5:
            is_source_free = False
        else:
            is_source_free = True

        if write_regions: write_region(region=source_free_region, region_name=f"{subdir_path}/{image_name}.source_free_region.annulus")

        # if there are no sources in the region, calculate the RMS
        if is_source_free:
            near_source_imstat_results = imstat(image_path, region=source_free_region)
            near_source_rms = near_source_imstat_results["rms"][0]
        
        # if there are, write a circle region with 50 beam areas and ask user to move it to a source free region
        else:

            # if overriding that, throw a warning but keep the annulus and throw a warning that this needs to be redone
            if override_sfr_request:

                near_source_imstat_results = imstat(image_path, region=source_free_region)
                near_source_rms = near_source_imstat_results["rms"][0]
                
                source_free_region_to_write = f"circle[[{x0}pix, {y0}pix], {(50**0.5)*beam_radius}arcsec]"
                write_region(region=source_free_region_to_write, region_name=f"{subdir_path}/{image_name}.source_free_region.circle.move_and_resave")
                warnings.warn(f"source-free region has potential source: {region_SNR} SNR found in source-free region. "
                              f"Overridden to keep code running, but this inaccurate RMS will be passed to imfit. "
                              f"The flux value for {image_name} will need to be retaken. "
                              f"Wrote circular region for you to move, resave, and pass into fit_point_source.")

            # if not overridden, raise error and prompt user to move the region that was written and pass it back in
            else:
                source_free_region = f"circle[[{x0}pix, {y0}pix], {(50**0.5)*beam_radius}arcsec]"
                write_region(region=source_free_region, region_name=f"{subdir_path}/{image_name}.source_free_region.circle.move_and_resave")
                raise ValueError(f"source-free region has potential source: {region_SNR} SNR found in source-free region. "
                                 f"Wrote region for you to move, resave, and pass into fit_point_source.")

    # if user provided source_free_region, find near-source RMS
    else:
        near_source_imstat_results = imstat(image_path, region=source_free_region)
        if override_sfr_request:
            near_source_rms = near_source_imstat_results["rms"][0]
            is_source_free = None
        else:
            region_SNR = imstat(image_path, region=source_free_region)["max"][0]/rms_image
            if region_SNR >= 5:
                is_source_free = False
                raise ValueError(f"The region you provided {source_free_region} has potential source: {region_SNR} "
                                 f"SNR found in source-free region")
            else:
                is_source_free = True
                near_source_rms = near_source_imstat_results["rms"][0]

    # get region area
    source_free_region_area = get_region_area(image_path, source_free_region, cell_size)
    source_free_region_area_beams = source_free_region_area/(np.pi*bmaj*bmin)

    # define region for fitting point source: an ellipse with axes {scale} times the beam size
    scale = 2.5
    region_flux = f"ellipse[[{x0}pix, {y0}pix], [{scale*bmaj/2}arcsec, {scale*bmin/2}arcsec], {bpa}deg]"
    if write_regions: write_region(region=region_flux, region_name=f"{subdir_path}/{image_name}.region_to_fit_point_source")

    # sometimes imfit will fail so going for try/except
    try:
        imfit_results = imfit(imagename=image_path, region=region_flux, rms=near_source_rms)
        result = "imfit success"
        
        # get imfit flux results
        flux = imfit_results["results"]["component0"]["flux"]["value"][0]
        flux_err = imfit_results["results"]["component0"]["flux"]["error"][0]
        flux_peak = imfit_results["results"]["component0"]["peak"]["value"]
        SNR = flux_peak/near_source_rms
    
        # get imfit shape results
        fx0, fy0 = imfit_results["results"]["component0"]["pixelcoords"]
        fmaj = imfit_results["results"]["component0"]["shape"]["majoraxis"]["value"]
        fmin = imfit_results["results"]["component0"]["shape"]["minoraxis"]["value"]
        fpa = imfit_results["results"]["component0"]["shape"]["positionangle"]["value"]
        point_source_fit_region = f"ellipse[[{fx0}pix, {fy0}pix], [{fmaj/2}arcsec, {fmin/2}arcsec], {fpa}deg]"
        if write_regions: write_region(region=point_source_fit_region, region_name=f"{subdir_path}/{image_name}.point_source_fit")
        
        # compare fit to beam: 
        # dArea is the fraction of 1 synthesized beam area the point source fit takes up
        # dPhi is the absolute offset of the position angles: fpa - bpa
        dArea, dPhi = compare_fit_to_beam(bmaj, bmin, bpa, fmaj, fmin, fpa)
    
        # find on-source RMS in case of non-detection
        region_nondetection = f"circle[[{x0}pix, {y0}pix], {(50**0.5)*beam_radius}arcsec]"
        nondetection_imstat_results = imstat(image_path, region=region_nondetection)
        on_source_rms = nondetection_imstat_results["rms"][0]

        # also fix the beam size and integrate over that
        '''
        fixed_x0 = x0
        fixed_y0 = y0
        fixed_maj = bmaj
        fixed_min = bmin
        fixed_pa = bpa
        init_guess_peak = flux_peak
    
        estfile = f"{directory}/estimates.txt"
        with open(estfile, "w") as f:
            #f.write(f"{init_guess_peak}, {fixed_x0}, {fixed_y0}, {fixed_maj}arcsec, {fixed_min}arcsec, {fixed_pa}deg, fxyabp")
            f.write("# x, y, peak, major, minor, pa, fit_x, fit_y, fit_peak, fit_major, fit_minor, fit_pa\n")
            f.write(f"{fixed_x0}, {fixed_y0}, {init_guess_peak}, {fixed_maj}arcsec, {fixed_min}arcsec, {fixed_pa}deg, abp\n")

        # results
        imfit_results_beam = imfit(imagename=image_path, rms=near_source_rms, estimates=estfile)
        flux_beam = imfit_results_beam["results"]["component0"]["flux"]["value"][0]
        flux_err_beam = imfit_results_beam["results"]["component0"]["flux"]["error"][0]
        flux_peak_beam = imfit_results_beam["results"]["component0"]["peak"]["value"]
        SNR_beam = flux_peak_beam/near_source_rms

        # shape: just make sure it's the same beam shape
        fmaj_beam = imfit_results_beam["results"]["component0"]["shape"]["majoraxis"]["value"]
        fmin_beam = imfit_results_beam["results"]["component0"]["shape"]["minoraxis"]["value"]
        fpa_beam = imfit_results_beam["results"]["component0"]["shape"]["positionangle"]["value"]
        dArea_beam, dPhi_beam = compare_fit_to_beam(bmaj, bmin, bpa, fmaj_beam, fmin_beam, fpa_beam)
        '''

    # if there is nothing there, imfit will fail. Reporting flux as 3x the RMS in a region over the source
    except:
        print(f"Imfit failed: check results carefully!")
        region_nondetection = f"circle[[{x0}pix, {y0}pix], {(50**0.5)*beam_radius}arcsec]"
        nondetection_imstat_results = imstat(image_path, region=region_nondetection)
        on_source_rms = nondetection_imstat_results["rms"][0]
        flux = 3*on_source_rms
        flux_peak = flux
        flux_err = 0
        SNR = 0
        point_source_fit_region = None
        dArea, dPhi = 0, 0
        result = "imfit fail: reporting flux as 3*on_source_RMS"

        flux_beam = 0
        flux_err_beam = 0
        flux_peak_beam = 0
        SNR_beam = 0
        dArea_beam, dPhi_beam = 0

    # store data
    results_dict = {"source": source, 
                    "Telescope": telescope, 
                    "Program": None, 
                    "Config": None,
                    "Observation Date": obs_date, 
                    "dT [days]": None, 
                    "Data source": "Jimmy", 
                    "Freq [Ghz]": freq, 
                    "Flux [mJy]": round(1000*flux, 4), 
                    "Flux_err [mJy]": round(1000*flux_err, 4), 
                    "dA [fit_area/beam]": round(dArea, 2), 
                    "dPhi [pa-pa_beam]": round(dPhi, 2),
                    "Flux_peak [mJy/beam]": round(1000*flux_peak, 4), 
                    "RMS [mJy/beam]": round(1000*near_source_rms, 4),
                    "SNR [-]": round(SNR, 1), 
                    "Detection": None, 
                    "Fit result": result, 
                    "Dynamic Range": dynamic_range, 
                    "Image process": None,
                    "Publish with?": None,
                    "Where is my image": None,
                    "Configuration": None}
    '''
                    "Flux_fixed_beam [mJy]": round(1000*flux_beam, 4),	
                    "Flux_err_fixed_beam [mJy]": round(1000*flux_err_beam, 4),	
                    "Flux_peak_fixed_beam [mJy/beam]": round(1000*flux_peak_beam, 4),
                    "SNR_fixed_beam [-]": round(SNR_beam, 1), 
                    "dArea_beam": round(dArea_beam, 2), 
                    "dPhi_beam": round(dPhi_beam, 2),
                    "RMS_whole_image": round(1000*rms_image, 4), 
                    "RMS_on_source": round(1000*on_source_rms, 4), }
    '''

    # print summary
    if print_results:
        print(f"{image_name}")
        max_len = max(len(k) for k in results_dict)
        for k, v in results_dict.items():
            if k in ["rms_image", "near_source_rms", "flux", "flux_err", "on_source_rms"]:
                print(f"{k:<{max_len}} : {round(v, 4)} mJy")
            else:
                if k in ["SNR", "dArea", "dPhi"]:
                    print(f"{k:<{max_len}} : {round(v,2)}")
                else:
                    print(f"{k:<{max_len}} : {v}")
        print("\n")

    if write_results:
        results_save = results_dict.copy()
        #results_save.pop("beam_region", None)
        #results_save.pop("point_source_fit_region", None)
        df_save = pd.DataFrame([results_save])
        df_save.to_csv(f"{subdir_path}/{image_name}.fit_results.csv")
    
    return results_dict

def run_fit_point_source(image_paths, source_free_region, print_results, write_results, write_regions, override_sfr_request, 
                         logfile_path, df_store):

    fit_results = []
    for image_path in image_paths:
        fit_result = fit_point_source(image_path+".image.tt0", 
                                       source_free_region=source_free_region, 
                                       print_results=print_results, 
                                       write_results=write_results,
                                       write_regions=write_regions,
                                       override_sfr_request=override_sfr_request)
    
        fit_results.append(fit_result)

    # write logfile
    if write_results: 
        df_results = pd.DataFrame(fit_results)
        df_save = pd.concat([df_store, df_results], axis=1)
        df_save.to_csv(logfile_path)

def convert_to_fits(image_path, image_prefix, crop_size=128):

    output_image_name = f"{image_prefix}.{crop_size}px.fits"

    # open image
    ia = image()
    ia.open(image_path)
    
    # get image coordinate system
    shape = ia.shape()  # [nx, ny, ...]
    nx, ny = shape[0], shape[1]
    
    # define pixel range for cropping
    half = crop_size // 2
    x1 = f"{pix_x - half}pix"
    y1 = f"{pix_y - half}pix"
    x2 = f"{pix_x + half - 1}pix"
    y2 = f"{pix_y + half - 1}pix"
    crop_region = f"box[[{x1}, {y1}], [{x2}, {y2}]]"
    
    # extract subimage
    subimage_name = 'cropped.image'
    if os.path.exists(subimage_name):
        os.system(f"rm -rf {subimage_name}")
    ia.subimage(outfile=subimage_name, region=crop_region)
    ia.close()
    
    # scale to be 99.95% 
    ia.open(subimage_name)
    data = ia.getchunk()
    ia.close()
    
    # get min and max for 99.95% clip
    #flat_data = data[0, 0, :, :].flatten()
    #low, high = np.percentile(flat_data, [0.05, 99.95])  # middle 99.9% range
    #scaled_data = np.clip(data[0, 0], low, high)

    # write to fits
    # make sure to flip CASA axes from [chan, stokes, y, x] to [y, x]
    output_image_path = f"{image_directory}/{output_image_name}"
    exportfits(subimage_name, fitsimage=output_image_path)
    print(f"Wrote {output_image_name} to fits")

    return output_image_name

'''
for image_path in image_paths:
    results_dict = fit_point_source(image_path, 
                                    source_free_region=source_free_region, 
                                    print_results=print_results, 
                                    write_results=write_results,
                                    write_regions=write_regions,
                                    override_sfr_request=override_sfr_request)
'''    

def get_theoretical_rms(listfile):
    return None

'''
def fit_point_source(image_name, estimates_file, region_flux, region_rms, region_non_detection, imfit_logfile):

    near_source_imstat_results = imstat(imagename=image_name, region=region_rms)
    near_source_rms = near_source_imstat_results["rms"][0]
    
    # if imfit works, it found something to fit to, but not necessarily your source
    # it's good to check manually by opening the image in CARTA
    try:
        rms = near_source_rms
        imfit_results = imfit(imagename=image_name, region=region_flux, rms=rms, logfile=imfit_logfile)
    
        flux = imfit_results["results"]["component0"]["flux"]["value"][0]
        flux_err = imfit_results["results"]["component0"]["flux"]["error"][0]
        flux_peak = imfit_results["results"]["component0"]["peak"]["value"][0]
        SNR = flux_peak/rms
    
        # just in case imfit does work but there isn't a 3 sigma detection
        if SNR >= 3:
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
        SNR = 0
        detection = False
        result = "imfit fail: flux is 3*RMS"

    # returning values in mJy
    return round(flux*1000, 4), round(flux_err*1000, 4), round(rms*1000, 4), round(SNR, 1), detection, result
'''