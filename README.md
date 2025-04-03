# Image calibrated measurement sets from the VLA
Scripts to make life easier when imaging calibrated measurement sets from the VLA. 

## Repository contents:
1. auto-image-singlefreq.py: script to create an image from a single-frequency VLA observation
2. (in progress) auto-image-multifreq.py: script to create a set of images from a multi-frequency VLA observation
3. vla-configuration-schedule.csv: table from https://science.nrao.edu/facilities/vla/proposing/configpropdeadlines
4. vla-resolution.csv: table from https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/resolution
5. config.example.yaml: an example file where you should define your imaging parameters. For how to use it, see note below.
6. .gitignore: ignore this

## Workflow of auto-image-singlefreq.py
1. User passes in path/to/measurement_set.ms, desired image size (default: 256 px), source name, VLA band, and any changes to CASA's tclean default values
2. Script creates a listfile and parses information for CASA's tclean task
3. Script creates an image using tclean
4. (Optional) Using CASA's imfit and imstat tasks, script fits for a point source flux value and RMS in circle centered at source coords. If no point source is detected, it reports 3*RMS within the circle as the upper limit.

## Workflow of auto-image-multifreq.py
(TODO): The multi-frequency script loops the single-frequency script to image several bands in one command. The user can specify if each band should be split into lower and upper spectral window bands, as is the case in our science. In the event that the imfit tasks finds non-detections in both lower and upper frequency windows for a single band, the script will re-image the measurement set for the full band to increase SNR.

## Note on config.example.yaml
The imaging scripts import user inputs through the config.yaml file. To avoid pushing my local changes to config.yaml and overriding the example I wanted to keep here, I created an example file called config.example.yaml. When cloning this repo, this file will of course come along, and you should make a copy and rename it config.yaml via 

  cp config.example.yaml config.yaml

Then edit the config.yaml file with your inputs. The imaging script does not look for your inputs in config.example.yaml!
