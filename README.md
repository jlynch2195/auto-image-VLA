# Image calibrated measurement sets from the VLA
Scripts to (allegedly) make life easier when imaging calibrated measurement sets from the VLA. 

## Getting started with GitHub
To obtain a local copy of this code + its supporting materials, clone the repository in Terminal via


    git clone https://github.com/jlynch2195/auto-image-VLA.git

You only have to do this once. To update your local copy to match the newest version on GitHub, make sure you're in the repository in your filesystem, then run

    git pull

You may get error/warning messages about needing to commit changes before pulling, meaning you've made edits to the existing files and git doesn't want to override those. If you want to keep your version, you can rename yours to avoid overwriting them. There's probably better practice; it's worth a Google.

## Repository contents:
1. auto-image-singlefreq.py: script to create an image from a single-frequency VLA observation
2. (in progress) auto-image-multifreq.py: script to create a set of images from a multi-frequency VLA observation
3. vla-configuration-schedule.csv: table from https://science.nrao.edu/facilities/vla/proposing/configpropdeadlines
4. vla-resolution.csv: table from https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/resolution
5. config.example.yaml: an example file where you should define your imaging parameters. For how to use it, see note below.
6. .gitignore: ignore this
7. requirements.txt: required packages, can be installed via pip install -r requirements.txt

## Workflow of auto-image-singlefreq.py
1. User passes in path/to/measurement_set.ms, desired image size (default: 256 px), source name, VLA band, and any changes to CASA's tclean default values
2. Script creates a listfile and parses information for CASA's tclean task
3. Script creates an image using tclean
4. (Optional) Using CASA's imfit and imstat tasks, script fits for a point source flux value and RMS in circle centered at source coords. If no point source is detected, it reports 3*RMS within the circle as the upper limit.

## Workflow of auto-image-multifreq.py
(TODO): The multi-frequency script loops the single-frequency script to image several bands in one command. The user can specify if each band should be split into lower and upper spectral windows, as is good practice if you have strong detections. In the event that the imfit tasks finds non-detections in both lower and upper frequency windows for a single band, the script will re-image the measurement set for the full band to increase SNR.

## Example
You have a calibrated measurement set called "main.ms", a C-band observation of SN2018cow, that you want to image for a quick look. You've cloned this repository and cd'ed into it on your local filesystem. You have the VLA pipeline version of CASA downloaded (https://casa.nrao.edu/casa_obtaining.shtml) and can start it by running "casa" in Terminal.

You first make a copy of config.example.yaml (see note below):

    cp config.example.yaml config.yaml

Then edit the config.yaml file to specify your parameters:

    measurement_set: "full/path/to/meaurement_set/main.ms"
    source_name: "SN2018cow"
    band: "C"
    image_size: 256
    
    use_manual_spws: False
    manual_spws: "2~32"
    
    try_point_source: True

That's all the setup that's required! Then, in terminal, run these commands to execute the script:

    (base) cd place/where/this/code/lives/
    (base) casa
    (CASA) <1>: execfile("auto-image-singlefreq.py")

If successful, this should print the following output in Terminal:

    "Created listfile full/path/to/meaurement_set/listfile.txt"
    "Begin imaging full/path/to/meaurement_set/SN2018cow_6.0GHz_256px"
    "Finished imaging full/path/to/meaurement_set/SN2018cow_6.0GHz_256px in X mins"

    "Fitting point source in region 8.75arcsec centered at RA, DEC" 
    "Detection at 6.0GHz: X +/- X_err mJy"

Note: it's still highly recommended to open the image in CARTA to manually inspect for image artifacts and validate the flux results.


## Note on config.example.yaml
The imaging scripts import user inputs through the config.yaml file. To avoid pushing my local changes to config.yaml and overriding the example I wanted to keep here, I created an example file called config.example.yaml. Locally, you should make a copy of this file and rename it config.yaml via 

    cp config.example.yaml config.yaml

Then edit the config.yaml file with your inputs. The imaging script does NOT look for your inputs in config.example.yaml!
