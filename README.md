# Image calibrated measurement sets from the VLA
Scripts to (allegedly) make life easier when imaging calibrated measurement sets from the VLA. 

## Getting started with GitHub
To obtain a local copy of this code + its supporting materials, clone the repository in Terminal via


    git clone https://github.com/jlynch2195/auto-image-VLA.git

You only have to do this once. To update your local copy to match the newest version on GitHub, make sure you're in the repository in your filesystem, then run

    git pull

You may get error/warning messages about needing to commit changes before pulling, meaning you've made edits to the existing files and git doesn't want to override those. If you want to keep your version, you can rename yours to avoid overwriting them. There's probably better practice; it's worth a Google.

## Repository contents:
1. run-auto-image.py: the main functionality; a script to create an image from any VLA observation, either single band or multi-frequency SEDs
2. run-rename-ms.py: renames your .ms file structure from the archive (eg. "program/pipeline.61139.014699073974/25A-060.sb48980598.eb49031908.60878.251126539355.ms") to "program/target/program.target.date/program.target.date.ms"
3. run-fit-point-source.py: also very useful; script to find flux of a point source, with both basic and advanced capabilities
4. vla-configuration-schedule.csv: table from https://science.nrao.edu/facilities/vla/proposing/configpropdeadlines
5. vla-resolution.csv: table from https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/resolution
6. config.example.yaml: an example file where you should define your imaging parameters. For how to use it, see note below.
7. .gitignore: ignore this
8. requirements.txt: required packages, can be installed via pip install -r requirements.txt

## run-rename-ms.py
The run-auto-image.py code relies on a very specific file structure at the moment. It expects that your .ms file lives in a directory with parents "~/program/target/program.target.date/program.target.date.ms". I do not feel like it is worth the time to generalize this at the moment. However, if your measurement set is fresh from the NRAO archive and has formatting similar to "~/program/pipeline.61139.014699073974/25A-060.sb48980598.eb49031908.60878.251126539355.ms", I at least wrote a script to turn that file structure to the required one.

Only run this once on your NRAO-formatted mspath by opening the "run-rename-ms.py" script, changing the "mspath" variable to the one you have (with the full path going back to ~/Users, etc), and running it in a casa terminal with command execfile("run-rename-ms.py"). Make sure that the correct "~/program/target/program.target.date/program.target.date.ms" got created in your filesystem, and now that directory is free to paste into the config.yaml file. The NRAO mspath can now be deleted to save space.

Some notes/errors:
* If the directory already exists, or if your mspath isn't directly from the archive, this will fail. At that point, it's probably time efficient to just manually rename your mspath with your target, observation date, and program, which is what we used to do before.
* If your target looks like a calibrator instead, that's because the target field is hard-coded as field = 2. To change it, create or open your listfile and find the correct target field ID. Then, in the rename_ms function inputs, change the field to the correct ID and rerun it.
* This is also an extra fail safe to make sure you actually know what you're imaging: the target, program, and obs_date come directly from the listfile, so the file directory that is setup reflects that. If the directories that the script creates are different from what you expect, double check you have the correct .ms file; it's easy to accidentally image the wrong one.

## Note on config.example.yaml
The run-auto-image.py script imports user inputs through the config.yaml file. To avoid pushing my local changes to config.yaml and overriding the example I wanted to keep here, I created an example file called config.example.yaml. Locally, you should make a copy of this file and rename it config.yaml via 

    cp config.example.yaml config.yaml

Then edit the config.yaml file with your inputs. The run-auto-image.py script does NOT look for your inputs in config.example.yaml!

## run-auto-image.py
The config.yaml file is where you specify inputs to the run-auto-image.py code. If the observation is a single-band, you must specify which band it is, and whether you would like an image for the entire band, two images for each of the half bands, or both. If the observation is multi-frequency, it will automatically determine the bands and image each of them, again at either the entire band, the half bands, or both based on user spec. For faint sources, I recommend using "both" so that if there are non-detections in each half band, the full band is automatically imaged as well and doesn't require a re-run. For multifrequency observations, the setup scans are assumed to be X-band and occupy spws 0 and 1, which are removed from X-band imaging. This is hard-coded in (for now) and may affect results if the setup scans are different.

The workflow is as follows. The script creates a listfile, scrapes it for spectral window ranges, and calculates the angular resolution based on the band and configuration. It then passes a dataframe of these values to a tclean helper that creates an image for each specified band/half-band. If try_point_source is toggled True, it will find the flux of your target source, assuming it is at or near the center of the image. More details on the fitting_procedure are below.

## run-fit-point-source.py
This script either takes in a list of images, or it will be automatically run on the images created during run-auto-image.py if try_point_source is toggled True in the config file.

If fitting_procedure is "basic", this script will attempt to fit a Gaussian the exact size and orientation of the observation's beam in a small region near the target source using the CASA task imfit. By restricting the fit size to exactly one beam, this accounts for only point source-like emission, which is typically of TDEs at our typical redshifts of >0.05-0.1. If imfit succeeds, it reports the integrated flux within the beam, and the error on the fit, which includes a contribution from the background RMS.

If there is no strong source, imfit may fail, in which case the flux is reported as 3x the RMS in a region ~100 beam areas at the target location. In some edge cases, there may be no strong source, but imfit may fit to a noise peak, returning a false positive. It is very important to open each image in CARTA to verify the fits. For this reason, the "detection" column is always left blank, and two sets of flux and flux errors are returned regardless of the actual result:
* Flux Det and Flux_err Det are the imfit results for exactly one beam at target location. These values are either real emission or from a noise peak, depending on how the image looks
* Flux Nondet is 3x the RMS at the target, and Flux_err Nondet is always zero
It is up to the user to determine if there is a detection, and hence which flux value to use: either a convincing detection with an error, or an upper limit!

If fitting_procedure is not "basic", a more advanced fitting procedure will take place. This allows for free-form fits (not restricted to beam shape or size) and allows the user to input source-free regions in the image to allow for a more accurate RMS reading. It can take x0 and y0 arguments to fit sources not at the center of the image. It will also write CARTA region files for the beam, the imfit Gaussian shape, and the region in which RMS is determined. This turned out to be overkill for our purposes and hence "basic" should be accurate enough in most cases.

## Example: single frequency observation
You have a calibrated measurement set called "main.ms", a C-band observation of SN2018cow, that you want to image for a quick look. You've cloned this repository and cd'ed into it on your local filesystem. You have the VLA pipeline version of CASA downloaded (https://casa.nrao.edu/casa_obtaining.shtml) and can start it by running "casa" in Terminal.

You first make a copy of config.example.yaml (see note below):

    cp config.example.yaml config.yaml

Then edit the config.yaml file to specify your parameters. You want a modest size image at a central frequency of 6 GHz (split: "whole) to maximize SNR. You specify EVLA_C as the band because the code needs that.

    measurement_set: "path/to/main.ms"
    source_name: "target"
    image_size: 512
    split: "whole"
    use_single_band: True
    single_band: "EVLA_C"

You want to get a flux measurement, but you don't need anything fancy. You just want the values printed in your terminal (print_results: True); no need for it to write a new Excel file with just a single line (write_results: False).

    try_point_source: True
    fitting_procedure: "basic"
    print_results: True
    write_results: False

You ignore the remaining parameters in config.yaml because "basic" fitting procedure doesn't utilize any of them. That's all the setup that's required! Then, in terminal, run these commands to execute the script:

    (base) cd place/where/this/code/lives/
    (base) casa
    (CASA) <1>: execfile("run-auto-image.py")

If successful, this should print a summary dataframe of the different images, their bands, spws, etc. Then it will notify you that it began imaging, and when it's finished, will print two lines (headers and values) that include flux values for what to report if a detection OR a non-detection. You open the image in CARTA, determine that there's a detection, and copy the values line into an Excel sheet, using a tab delimiter. You are happy that you didn't have to go to the VLA website to find the cell size or manually draw ellipses for fitting.

## Example: multi-frequency observation
You have a calibrated measurement set called "main-multi.ms", an LSCX-bands observation of SN2018cow, that you want to image to find information about the spectral energy distribution (SED). You've cloned this repository and cd'ed into it on your local filesystem. You have the VLA pipeline version of CASA downloaded (https://casa.nrao.edu/casa_obtaining.shtml) and can start it by running "casa" in Terminal.

You first make a copy of config.example.yaml (see note below):

    cp config.example.yaml config.yaml

Then edit the config.yaml file to specify your parameters. You want small images to speed up the process, and you want to split the full-band images into half bands, so that you can get two data points at different frequencies per band. You are pretty sure that it's bright enough to be detected at each half-band, but you toggle split: "both" to also image the full band: better to do this now just in case, for example, the 9 and 11 GHz images have faint detections that may form a much better detection when imaged at the full band! You are doing an SED, so you toggle use_single_band: False, and you leave single_band as whatever you want.

    measurement_set: "path/to/main.ms"
    source_name: "target"
    image_size: 256
    split: "both"
    use_single_band: False
    single_band: "EVLA_C"

You want to get a flux measurement for each image. This time, you want the results saved to an Excel file so that it's easy to copy over the entire block of fluxes and flux errors (write_results: True). This creates the file all_fit_results.csv in the same directory where all the images live. 

    try_point_source: True
    fitting_procedure: "basic"
    print_results: True
    write_results: True

You ignore the remaining parameters in config.yaml because "basic" fitting procedure doesn't utilize any of them. That's all the setup that's required! Then, in terminal, run these commands to execute the script:

    (base) cd place/where/this/code/lives/
    (base) casa
    (CASA) <1>: execfile("run-auto-image.py")

If successful, this should print a summary dataframe of the different images, their bands, spws, etc. which should have 12 lines (8 half bands, 4 full bands). Imaging should all take place first, and the fluxes are measured on each image only at the end. 

As an example, say you open the all_fit_results.csv file that shows quite large Flux Det values (~0.300 +/- 0.015 mJy) for L-band, but with a slight decline as frequency increases. In X-band, at 11 GHz, you find that the Flux Det value is 0.040 +/- 0.010 mJy, but since it is so faint, you open the image in CARTA just to be sure that your target is distinguishable from a noise peak. You can't convince yourself that it is, so you check 9 GHz, to the same conclusion: you do not have a detection at these frequencies. For these images, you report the flux value as the Flux Nondet, a few columns down in the Excel. For SED purposes, all is not lost: you imaged the full band (at 10 GHz), and inspection of the image yields a detection, which is due to increased SNR from "combining" the power of both the 11 and 9 GHz images. 
