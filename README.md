# Fort Collins snowstorm maps
This repository includes the code required to use data from the ERA5 (ECWMF Reanalysis) API, along with [MetPy](https://unidata.github.io/MetPy/latest/index.html) and other python tools, to produce a series of maps for snowstorms that have occurred in Fort Collins, Colorado since 1979.  My images are on the web at [http://schumacher.atmos.colostate.edu/teaching/ats641/snowstorms](http://schumacher.atmos.colostate.edu/teaching/ats641/snowstorms).  Some instructions for building these are given below, and the code should be flexible enough that by providing your own list of cases that you'd like to look at, you could make similar maps for other cases of interest.

## First steps

This assumes that you have a working version of python on your computer and have a number of standard scientific packages installed.  If not, I suggest following the instructions used in the [Unidata python workshop](https://unidata.github.io/python-training/) to get started.

This code also requires you to get an account for ECMWF's "Climate Data Store" API.  Instructions for getting an account and setting up a few other things can be found at [this link](https://cds.climate.copernicus.eu/api-how-to).

## Which cases are we interested in?

I built this to analyze storms with at least 12" of total snowfall, where a 'storm' is defined as one or more consecutive days of measurable snowfall, at the Fort Collins weather station since 1979.  It plots maps starting four days before the first day of the storm, through four days after the first day of the storm; with maps twice a day at 00 and 12 UTC.  Here, these are listed in a file named `snowstorms_12inch_fortcollins_post1979.csv`, which also includes the total and daily snowfall for each storm.  *But here's the easy part: all you need is a text or csv file with a single column of dates to read in, and the code will make plots for your list of cases.*  And wait, there's more! This does not actually download ERA5 files (instead, it reads them remotely), so it will only take up minimal disk space. 

## The code

This is organized into one giant python program, `winterstorms_era5_loop.py` (or an updated version, `winterstorms_era5_metpyV1_loop.py`).  This may not be the finest programming practice but it will do everything in a single step.  If you wanted to modify this for your own purposes, there are just a few things that you might need to change:

1) The filename that is read in as `case_data_all` (currently `snowstorms_12inch_fortcollins_post1979.csv`)
2) The different fields to plot have 'True/False' settings for whether you want them plotted or not
3) The number of dates and times you want to plot.  Currently, the code finds the start date of the event from the csv file, then starts plotting four days before that date (`n_days_before=4`) and continues through four days after (`n_days_after=4`). If you wanted more or fewer days, just change the 4 to something else.  Right now it also is set to only get data at 00 and 12 UTC from ERA5, but ERA5 has files every hour, so in your calls to the ECWMF API you could add more times.
4) lat/lon areas for downloaded data and/or plots.  These you will need to adjust in various places in the code, both in the call to the ECWMF API, and also for individual plots (which have varying lat/lon bounds).
5) What fields are plotted or details of the plots. These are at your discretion to customize.

If you're working through the code and making modifications, etc., you may want to copy chunks into a jupyter notebook - this is how I wrote the code.  


## Running the code

Once everything is in place, you should be able to simply run the code with:
`python winterstorms_era5_loop.py` 
from the command line.  It will automatically create subdirectories for each storm and will write the png image files to those subdirectories.


## Acknowledgments

Huge thanks to [Unidata](https://unidata.github.io/python-training/) for their training materials and tools. 
