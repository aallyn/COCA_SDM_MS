# COCA_SDM_MS
Public repository for data used in the COCA SDM methods paper (Allyn et al.)

# Overview
The code and data provided here support our analysis to develop a novel method that blends quantitative distribution model projections with qualitative expert climate vulnerability assessment expectations to more accurately project marine species distribution and abundance under future climate scenarios in the Northeast U.S. Large Marine Ecosystem. To reach this overall goal, many interim analysis steps were completed, including (1) developing two-stage delta lognormal generalized additive models for a suite of marine species, (2) rigorously validating the predictive performance of the distribution models, (3) comparing quantitative distribution models projections using sea surface temperatures from an ensemble of climate models run under the RCP8.5 “business as usual” scenario to expectations from a expert climate vulnerability assessment, and, (4) proposing and evaluating a novel method for combining the quantitative distribution models with qualitative expert assessments. 

# Details
Before trying to execute the code, make sure to unzip all of the data files in the Data folder so that the only sub-folder within the Data folder is "BottomTrawlStrata" as in: Data/BottomTrawlStrata. The only other road block to running the code (I hope) is getting the OISST data. These are publicly available [here](https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html). Or, feel free to [send me an email](mailto:aallyn@gmri.org) and I'll be happy to pass along the netcdf file we created specifically for this study -- the file was too big to upload to GitHub.
