# Image calibrated measurement sets from the VLA

Repository contains:
1. auto-image-singlefreq.py: script to create an image from a single-frequency VLA observation
2. (in progress) auto-image-multifreq.py: script to create a set of images from a multi-frequency VLA observation
3. vla-configuration-schedule.csv: table from https://science.nrao.edu/facilities/vla/proposing/configpropdeadlines
4. vla-resolution.csv: table from https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/resolution

The workflow for the imaging scripts are as follows:
1. User passes in path/to/measurement_set.ms, desired image size (default: 128 px), and any changes to CASA's tclean default values
2. Script creates a listfile and parses information for CASA's tclean task + creates an image using tclean
3. (optionally) Using CASA's imfit and imstat tasks, fits for a point source flux value and RMS in circle centered at source coords. If no point source is detected, it reports 3*RMS within the circle as the upper limit.

The multi-frequency script loops the single-frequency script to image several bands in one command. The user can specify if each band should be split into lower and upper spectral window bands, as is the case in our science. In the event that the imfit tasks finds non-detections in both lower and upper frequency windows for a single band, the script will re-image the measurement set for the full band to increase SNR.
