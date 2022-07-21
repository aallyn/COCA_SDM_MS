#####
## Creating netcdf files for hosting on NOAA's DisMap portal
#####
library(tidyverse)
library(sf)
library(raster)
library(zoo)
library(rgeos)
library(mgcv)
library(tidync)

# New SF nonsense
sf_use_s2(FALSE)

proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19

## Interested in yearly predictions out until around mid-century time period. In the COCA work, we just did a snap shot out to 2055 and then to 2100. So, going to need to go back and grab the fitted models, get the right temperature data, and then make the projections.

## Directories
mods_dir <- "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/COCA15_ClimVuln/COCA-SDM/Results/NormalVoting_BiomassIncPresNoExposure_03152019"
oisst.dir<-  "/Users/aallyn/Library/CloudStorage/Box-Box/RES_Data/OISST/ThroughFeb2020.grd"
rcp85.mu.dir<- paste("/Users/aallyn/Library/CloudStorage/Box-Box/RES_Data/", "CMIP5_SST/ECW_RCP85_mu.grd", sep = "")
sp.in<- "/Users/aallyn/Library/CloudStorage/Box-Box/RES_Data/Shapefiles/" 

## Original fitted data -- going to need this for rescaling the future projections
# Fish assessment species
fish.spp<- read.csv(here::here("Data", "Assesmentfishspecies.csv"))

# Read it in, do some quick formatting to fit the GAM
dat<- readRDS(here::here("Data", "model.dat.rds"))%>% 
  filter(., SVSPP %in% fish.spp$SVSPP) %>%
  left_join(., fish.spp, by = "SVSPP") 

# Additional processing/filtering -- remove Atlantic sturgeon as an endangered species
dat<- dat[dat$COMNAME != "ATLANTIC STURGEON", ]

# Training vs. testing
train.start<- "1982-01-01"
train.end<- "2010-12-31"
test.start<- "2011-01-01"
test.end<- "2016-01-01"

dat$TRAIN.TEST<- ifelse(as.Date(dat$DATE) >= train.start & as.Date(dat$DATE) <= train.end, "TRAIN", 
                        ifelse(as.Date(dat$DATE) >= test.start & as.Date(dat$DATE) <= test.end, "TEST", "Neither"))

# Bottom trawl strata
bstrat<- st_read(here::here("Data/BottomTrawlStrata/", "BTS_Strata.shp"))

# Get names of strata
bstrat.names<- unique(bstrat$STRATA)

# Reduce dataset to use tows within the bottom trawl strata
dat<- dat[dat$STRATUM %in% bstrat.names,]

# Training data by season: fall and then spring, scaling variables before fitting models and defining "presence" response for presence model component
dat.train.f<- dat %>%
  filter(., TRAIN.TEST == "TRAIN" & SEASON == "FALL") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST))) 

dat.train.f$PRESENCE.BIOMASS<- ifelse(dat.train.f$BIOMASS > 0, 1, 0)

# Need to keep mean and sd from rescale to use when we predict or project to other time periods
base.depth.mean.f<- mean(abs(dat.train.f$DEPTH))
base.depth.sd.f<- sd(abs(dat.train.f$DEPTH))
base.shelf.mean.f<- mean(dat.train.f$SHELF_POS)
base.shelf.sd.f<- sd(dat.train.f$SHELF_POS)
base.temp.mean.f<- mean(dat.train.f$SEASONALMU.OISST, na.rm = T)
base.temp.sd.f<- sd(dat.train.f$SEASONALMU.OISST, na.rm = T)
fall.rescale.df<- data.frame(SEASON = "FALL", mean.t = base.temp.mean.f, sd.t = base.temp.sd.f, mean.depth = base.depth.mean.f, sd.depth = base.depth.sd.f, mean.shelf = base.shelf.mean.f, sd.shelf = base.shelf.sd.f)

# Now spring...
dat.train.s<- dat %>%
  filter(., TRAIN.TEST == "TRAIN" & SEASON == "SPRING") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST))) 

dat.train.s$PRESENCE.BIOMASS<- ifelse(dat.train.s$BIOMASS > 0, 1, 0)

# Temps to rescale other variables
base.depth.mean.sp<- mean(abs(dat.train.s$DEPTH))
base.depth.sd.sp<- sd(abs(dat.train.s$DEPTH))
base.shelf.mean.sp<- mean(dat.train.s$SHELF_POS)
base.shelf.sd.sp<- sd(dat.train.s$SHELF_POS)
base.temp.mean.sp<- mean(dat.train.s$SEASONALMU.OISST, na.rm = T)
base.temp.sd.sp<- sd(dat.train.s$SEASONALMU.OISST, na.rm = T)
spring.rescale.df<- data.frame(SEASON = "SPRING", mean.t = base.temp.mean.sp, sd.t = base.temp.sd.sp, mean.depth = base.depth.mean.sp, sd.depth = base.depth.sd.sp, mean.shelf = base.shelf.mean.sp, sd.shelf = base.shelf.sd.sp)

## Getting the future projected temperatures
seasons<- c("Spring", "Fall")
season <- seasons[2]

## Projections
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19
  
# Empty stack
pred.rast.stack<- stack()
  
# Add OISST
name.ind<- nlayers(pred.rast.stack)+1
stack0<- raster::stack(oisst.dir)
  
# Move to monthly?
oisst.min.date<- as.Date("1981-09-01")
oisst.max.date<- as.Date("2020-03-09")
  
# Calculate monthly mean temperature -- this would be compared to the sstclim data (monthly climate ensemble)
oisst.dates<- seq.Date(from = oisst.min.date, to = oisst.max.date, by = "day")
oisst.dat<- setZ(stack0, oisst.dates)
  
# Aggregate daily to monthly data
oisst.monthly <- zApply(oisst.dat, by = as.yearmon, mean)

# Basline
sst.stack<- stack()

years<- c("2014", "2015", "2016", "2017", "2018")

for(i in seq_along(years)){
    
    dates.use<- switch(season,
                         "Fall" = paste(c("Sep", "Oct", "Nov"), rep(years[i]), sep = "."),
                         "Spring" = paste(c("Mar", "Apr", "May"), rep(years[i]), sep = "."),
                         "Summer" = paste(c("Jun", "Jul", "Aug", "Sep"), rep(years[i]), sep = "."))
    
    sst.temp<- calc(oisst.monthly[[which(names(oisst.monthly) %in% dates.use)]], mean)
    
    sst.stack<- stack(sst.stack, sst.temp)
    print(years[i])
}

names(sst.stack)<- paste(season, years, sep = ".")

sst.basemeans <- calc(sst.stack[[c(1, 2, 3, 4, 5)]], mean)
pred.rast.stack <- stack(pred.rast.stack, sst.basemeans)
names(pred.rast.stack)[c(1)]<- c("Baseline")

# Climate
name.ind <- nlayers(pred.rast.stack) + 1
stack0<- stack(rcp85.mu.dir)
proj4string(stack0)<- proj.wgs84
stack0<- resample(stack0, oisst.monthly[[1]])
clim.dates <- seq.Date(from = as.Date("1980-01-16"), to = as.Date("2060-01-16"), by = "month")
clim.stack <- setZ(stack0, clim.dates)
clim.stack.zind<- getZ(clim.stack)

# Baseline stack, store seasonal means and then average them all
rcp85.mu.stack <- stack()

years <- seq(from = 1982, to = 2055, sep = "")

for(i in seq_along(years)){
    
    dates.use<- switch(season,
                          "Fall" = as.Date(paste(rep(years[i]), c("09-16", "10-16", "11-16"), sep = "-")),
                          "Spring" = as.Date(paste(rep(years[i]), c("03-16", "04-16", "05-16"), sep = "-")),
                         "Summer" = as.Date(paste(rep(years[i]), c("06-16", "07-16", "08-16", "09-16"), sep = "-")))
    clim.mu.temp<- calc(clim.stack[[which(clim.stack.zind %in% dates.use)]], mean)
     
    rcp85.mu.stack<- stack(rcp85.mu.stack, clim.mu.temp)
    print(years[i])
}

names(rcp85.mu.stack) <- paste(season, years, "rcp85.mu", sep = ".")

# Add it to the pred rast
pred.rast.stack<- stack(pred.rast.stack, rcp85.mu.stack)

# Other predictors
# Add depth
depth.temp<- raster(paste(sp.in, "NEShelf_Etopo1_bathy.tiff", sep = ""))
DEPTH<- resample(depth.temp, pred.rast.stack[[1]])
pred.rast.stack<- stack(pred.rast.stack, DEPTH)
    
# Add TRI
TRI.temp<- terrain(DEPTH, opt = "TRI")
TRI<- resample(TRI.temp, pred.rast.stack[[1]])
pred.rast.stack<- stack(pred.rast.stack, TRI)
names(pred.rast.stack)[76:77]<- c("DEPTH", "TRI")

pred.rast.stack<- projectRaster(pred.rast.stack, crs = proj.utm)
    
# Get these variables
# Mask out points outside of NELME
nelme.rast <- pred.rast.stack[[1]]
nelme.rast[] <- NA
nelme<- read_sf(paste(sp.in, "NELME_regions/NELME_sf.shp", sep = ""))
nelme.utm<- st_transform(nelme, proj.utm)
nelme.buff<- gBuffer(as_Spatial(nelme.utm), width = 40000)
nelme.rast<- rasterize(nelme.buff, nelme.rast)
pred.rast.stack.m<- raster::mask(pred.rast.stack, mask = nelme.rast, inverse = FALSE)
    
# Now values
pred.rast.stack.m.wgs84<- projectRaster(pred.rast.stack.m, crs = proj.wgs84)
points<- data.frame(coordinates(pred.rast.stack.m.wgs84[[1]]))
coordinates(points)<- ~x+y
proj4string(points) <- proj.wgs84

pred.df<- data.frame(raster::extract(pred.rast.stack.m.wgs84, points))

if(season == "Spring"){
    spring.rast.pred <- data.frame(points, pred.df)
    spring.rast.pred$SEASON <- "SPRING"
    saveRDS(spring.rast.pred, file = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/COCA15_ClimVuln/COCA-SDM/Data/spring.rast.preds.DisMap.rds")
} else if(season == "Fall"){
    fall.rast.pred <- data.frame(points, pred.df)
    fall.rast.pred$SEASON <- "FALL"
    saveRDS(fall.rast.pred, file = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/COCA15_ClimVuln/COCA-SDM/Data/fall.rast.preds.DisMap.rds")
}

## Getting the right scaling
temp.scale<- function(new.temp, base.temp.mean, base.temp.sd){
  if(is.na(new.temp)){
    temp.scaled<- NA
    return(temp.scaled)
  } else {
    temp.scaled<- (new.temp - base.temp.mean)/base.temp.sd
    return(temp.scaled)
  }
}

spring.preds<- readRDS("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/COCA15_ClimVuln/COCA-SDM/Data/spring.rast.preds.DisMap.rds")
fall.preds<- readRDS("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/COCA15_ClimVuln/COCA-SDM/Data/fall.rast.preds.DisMap.rds")

base.preds.sp<- spring.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH) %>%
  mutate(., "SEASON" = rep("SPRING", nrow(.)))

base.preds.sp<- base.preds.sp %>%
  left_join(., spring.rescale.df, by = "SEASON")
base.preds.sp$DEPTH.Scale<- mapply(temp.scale, abs(base.preds.sp$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
base.preds.sp$SEASONALMU.OISST.Scale<- mapply(temp.scale, base.preds.sp$Baseline, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

base.preds.sp<- base.preds.sp %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH) %>%
  mutate(., "SEASON" = rep("FALL", nrow(.))) 
base.preds.f$DEPTH.Scale<- mapply(temp.scale, abs(base.preds.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
base.preds.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, base.preds.f$Baseline, fall.rescale.df$mean.t, fall.rescale.df$sd.t)
base.preds.f<- base.preds.f %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds<- base.preds.f %>%
  bind_rows(., base.preds.sp)

## Future
fut.preds.sp <- spring.preds %>%
    mutate(., "SEASON" = rep("SPRING", nrow(.))) %>%
    left_join(., spring.rescale.df, by = "SEASON")
fut.preds.sp$DEPTH.Scale<- mapply(temp.scale, abs(fut.preds.sp$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)

scale_cols<- 4:77
for(i in seq_along(scale_cols)){
    fut.preds.sp[, scale_cols[i]]<- mapply(temp.scale, fut.preds.sp[, scale_cols[i]], spring.rescale.df$mean.t, spring.rescale.df$sd.t)
}

fut.preds.sp<- fut.preds.sp %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds.f<- fall.preds %>%
  mutate(., "SEASON" = rep("FALL", nrow(.))) %>%
  left_join(., fall.rescale.df, by = "SEASON")
fut.preds.f$DEPTH.Scale<- mapply(temp.scale, abs(fut.preds.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)

for(i in seq_along(scale_cols)){
    fut.preds.f[, scale_cols[i]]<- mapply(temp.scale, fut.preds.f[, scale_cols[i]], fall.rescale.df$mean.t, fall.rescale.df$sd.t)
}

fut.preds.f<- fut.preds.f %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds<- fut.preds.f %>%
  bind_rows(., fut.preds.sp)

## Get the model fits and make the projections
sp.sst.cols<- paste("Spring", seq(from = 1982, to = 2055), "rcp85.mu", sep = ".")
f.sst.cols<- paste("Fall", seq(from = 1982, to = 2055), "rcp85.mu", sep = ".")

res.files<- list.files(mods_dir, "mcmc_")
res_ind<- 1

for(i in seq_along(res.files)){
  spp<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][1])
  season<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][2])
  spp.season.match<- paste(spp, season, sep = ".")
  
  # Model fit -- presence and biomass SDM
  mod.fitted.p<- readRDS(paste(mods_dir, gsub("mcmc_", "gamfitpres", res.files[i]), sep = "/"))
  mod.fitted.b<- readRDS(paste(mods_dir, gsub("mcmc_", "gamfitbio", res.files[i]), sep = "/"))
  ilink<- family(mod.fitted.b)$linkinv
  gam.coef <- names(coef(mod.fitted.p))

  # Make predictions
  sst.cols.use <- switch(season,
      "SPRING" = sp.sst.cols,
      "FALL" = f.sst.cols
  )

  for(j in seq_along(sst.cols.use)){
      newdat.mu<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
      newdat.mu <- newdat.mu %>%
          unnest() %>%
          dplyr::select(., c("x", "y", "DEPTH.Scale", sst.cols.use[j]))
      names(newdat.mu)[4] <- "SEASONALMU.OISST.Scale"
      pred.p<- round(predict.gam(mod.fitted.p, newdata = newdat.mu, type = "response"), 2)
      pred.b<- round(as.numeric(pred.p) * exp(as.numeric(predict.gam(mod.fitted.b, newdata = newdat.mu, type = "response"))), 2)

      # Gather up
      pred.out.temp <- data.frame("COMNAME" = spp, "SEASON" = season, "Proj.Class" = sst.cols.use[j], "x" = newdat.mu$x, "y" = newdat.mu$y, "Projection" = pred.b) %>%
          group_by(COMNAME, SEASON, Proj.Class) %>%
          nest()
          
          if (res_ind == 1) {
              all_proj_out <- pred.out.temp
          } else {
              all_proj_out<- bind_rows(all_proj_out, pred.out.temp)
          }
      res_ind <- res_ind + 1
      }
      
      print(round(i / length(res.files)), 2)
  
}

# Source functions
source("https://raw.githubusercontent.com/GMRI-SEL/LabFunctionsandCode/master/GenerateSharedPathsScript.R")
source("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Functions/Raster_To_NetCDF_Functions.R")

# Read in SDM results and a little processing
results<- all_proj_out

mod.res<- read_csv(paste(mods_dir, "mod.results.csv", sep = "/")) %>%
  drop_na(RMSE.SDM.B, CorrCoeff.SDM.B, Bias.SDM.B, RMSE.NEVA.B, CorrCoeff.NEVA.B, Bias.NEVA.B)

# Exploring cut offs... AUC > 0.7 in both seasons
mod.spp.keep<- mod.res %>%
  filter(., AUC.SDM >= 0.7 & CorrCoeff.SDM.B >= 0) %>%
  group_by(., COMNAME) %>%
  summarize_at(vars(SEASON), n_distinct) %>%
  filter(., SEASON == 2)

mod.res<- mod.res %>%
  filter(., COMNAME %in% mod.spp.keep$COMNAME)

results.sub<- results %>%
  dplyr::filter(., COMNAME %in% unique(mod.res$COMNAME))

# Dataframes to raster stack 
results.stack.list<- results.sub %>%
  mutate(., "Raster" = map(data, tibb_to_rast))
results.stack<- raster::stack(results.stack.list$Raster)
names(results.stack)<- paste(as.character(results.sub$COMNAME), as.character(results.sub$SEASON), as.character(results.sub$Proj.Class), sep = "/")

# Species subset?
spp<- gsub(" ", ".", unique(as.character(results.sub$COMNAME)))
#spp<- c("AMERICAN.PLAICE", "ATLANTIC.COD", "ATLANTIC.HALIBUT", "POLLOCK", "ATLANTIC.WOLFFISH", "HADDOCK", "OCEAN.POUT", "WHITE.HAKE", "WINTER.FLOUNDER", "WITCH.FLOUNDER", "YELLOWTAIL.FLOUNDER")

# Output file stem
out.file.stem<- "~/Desktop/ProjectionsForDisMap/"
# Fill value
fillvalue.use = 1e32
# Coder
coder.use = "Andrew Allyn"

# For each species, output netcdf files
for(i in seq_along(spp)){

  # Get species, set netcdf file title and output
  spp.use<- spp[i]
  title.use<- gsub("[.]",  " ", str_to_sentence(spp.use))
  out.file.use<- paste(out.file.stem, title.use, ".nc", sep = "")

  # Get species raster layers from full stack
  raster.temp<- results.stack[[which(grepl(spp.use, names(results.stack)))]]
  na.value<- -99999
  raster.temp[is.na(raster.temp)]<- na.value
  raster.temp[is.infinite(raster.temp)]<- na.value
  raster.out <- raster.temp
  proj4string(raster.out) <- proj.wgs84

  names_list <- str_split(names(raster.temp), pattern = "[.]")
  names_use_temp <- do.call("rbind", lapply(names_list, "[", c(4, 5)))
  names_use <- paste(names_use_temp[, 1], names_use_temp[, 2], sep = "_")
  names(raster.out) <- names_use
  writeRaster(raster.out, out.file.use, overwrite = TRUE, format = "CDF", varname = "Biomass", varunit = "kg*tow-1", longname = "Predicted biomass (kg*tow-1)", xname = "Longitude", yname = "Latitude", zname = "Time (Season-Year)"
    )

  # print
  print(paste(spp.use, " is done", sep = ""))
}

# Quick check...
t <- tidync(out.file.use)
hyper_grids(t)





