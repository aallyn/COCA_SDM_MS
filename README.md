# COCA_SDM_MS
Public repository for data used in the COCA SDM methods paper (Allyn et al.)

# Overview
The code and data provided here support our analysis to develop a novel method that blends quantitative distribution model projections with qualitative expert climate vulnerability assessment expectations to more accurately project marine species distribution and abundance under future climate scenarios in the Northeast U.S. Large Marine Ecosystem. To reach this overall goal, many interim analysis steps were completed, including (1) developing two-stage delta lognormal generalized additive models for a suite of marine species, (2) rigorously validating the predictive performance of the distribution models, (3) comparing quantitative distribution models projections using sea surface temperatures from an ensemble of climate models run under the RCP8.5 “business as usual” scenario to expectations from a expert climate vulnerability assessment, and, (4) proposing and evaluating a novel method for combining the quantitative distribution models with qualitative expert assessments. 

# Details
## Code files
You should only need to open and run the code "sdm_neva_merge_biomass_ind.R" to reproduce all of the analyses and corresponding figures.

## Data files
In the data folder there are the following files:
1. Assesmentfishspecies.csv - A csv file with common and scientific names for the species assessed by Hare et al. (2016) in their climate vulnerability assessment for the Northeast U.S. Large Marine Ecosystem.
2. BottomTrawlStrata.zip - A zipped file of the NOAA Northeast Fisheries Science Center bottom trawl survey strata polygon shapefile.
3. JHareDirectionalEffect.csv - A csv file with the directional effect voting results by species from Hare et al. (2016).
4. JHareQualitativeDataResult.csv - A csv file with the exposure and sensitivity voting results by species from Hare et al. (2016).
5. JHareSppFunctionalGroup.csv - A csv file with species name and then assigned functional group by Hare et al. (2016).
6. NELME_clipped_shp.zip - A zipped file of the Northeast U.S. Large Marine Ecosystem border polygon shapefile.
7. RegionalShapefiles.zip - A zipped file of the Gulf of Maine and Southern New England/Mid-Atlantic bight borders polygon shapefile. 
8. SST.CMIP5.1982-2099.anom.nc.zip - A zipped netcdf file with the CMIP5 RCP8.5 ensemble mean, 5th and 95th percentile projected temperature anomalies.
9. model.dat.rds.zip - A zipped Rdata folder with the NOAA Northeast Fisheries Science Survey spring and fall bottom trawl survey data and environmental characteristics of each trawl (i.e., depth and seasonal average sea surface temperature). 
10. projections.dat.rds.zip -  A zipped Rdata folder with the environmental data used to project species distribution and abundance in 2055.

## Running the code
Before trying to execute the code, make sure to unzip all of the data files in the Data folder so that the only sub-folder within the Data folder is "BottomTrawlStrata" as in: Data/BottomTrawlStrata. The only other road block to running the code (I hope) is getting the OISST data. These are publicly available [here](https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html). Or, feel free to [send me an email](mailto:aallyn@gmri.org) and I'll be happy to pass along the netcdf file we created specifically for this study -- the file was too big to upload to GitHub.
