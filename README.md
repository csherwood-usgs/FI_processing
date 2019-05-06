# FI_processing
Matlab scripts to process Fire Island oceanographic data. This repo is in ..\proj\2019_FireIslandAsymmetry\FI_processing on my desktop.
---
### Overall goals of this repo
* Compare asymmetry deduced from bulk wave statistics via Ursell number with asymmetry measure with near-bottom instruments
* But first, demonstrate/test routines for reading/processing instrument data
---
### Asymmetry
* Measures of asymmetry vary. Some are best calculated for a singe wave, others are amenable to field measurements.
* TODO: summarize asymmetry formulae and code
---
### Progress so far
* The stats file has some missing depths...we skip these in processing, but have not investigated the underlying bursts. Some of the bursts where depth is missing are also missing some pressure measurements. These break the calcs.
* TODO: Add some basic QA/QC to the bursts to replace missing or bad data where possible, and to skip bursts that don't meet QA standards.
* TODO: Complete analysis of rotation into wave-propagation direction. I think `wds` returns the correct direction, and `puvq` returns a direction that is 180 deg. out...the direction waves are coming *from*.
* TODO: Complete calcs. of comparable statistics. Best version of *du/dt* for acceleration statistics? Hilbert transform (per Malarky and Davies, 2012) for acceleration asymmetry?
### Code so far:
* `puv_proc_FI`    - Main script to read in ADV data and process burst by burst
* `puvq`           - Function to calculate wave parameters from burst measurements of *u, v, p*
* `ubstatsr`       - Function to calculate a range of burst statistics
* `ruessink_asymm` - Function to calculate asymmetry parameters from Ursell number
* `urotate`        - Unfinished function to rotate velocities into direction of wave propagation
* `ustats`         - Unfinished function to calculate statistics on u returned by `urotate`
---
### Nortek scripts
We thank NortekUSA for several scripts written by Lee Gordon and distributed here:

http://www.nortekusa.com/usa/knowledge-center/table-of-contents/waves

These are:

`wds.m` - Directional wave spectra calcs

`hs.m`  - Summary statistics of output of `wds`

`wavek` and `logavg` - called by `wds`
