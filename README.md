# Fort Collins snowstorm maps
This repository includes the code required to use data from the ERA5 (ECWMF Reanalysis) API, along with [MetPy](https://unidata.github.io/MetPy/latest/index.html) and other python tools, to produce a series of maps for snowstorms that have occurred in Fort Collins, Colorado since 1979.  Some instructions are given below, and the code should be flexible enough that by providing your own list of cases that you'd like to look at, you could make similar maps for other cases of interest.

## First steps

This assumes that you have a working version of python on your computer and have a number of standard scientific packages installed.  If not, I suggest following the instructions used in the [Unidata python workshop](https://unidata.github.io/python-training/) to get started.

This code also requires you to get an account for ECMWF's "Climate Data Store" API.  Instructions for getting an account and setting up a few other things can be found at [this link](https://cds.climate.copernicus.eu/api-how-to).



