################ 
## COCA Methods Manuscript
################
# Preliminaries -----------------------------------------------------------
## Key components: No longer sampling using MCMC update -- exhaustive evaulation of all potential GAM models. Adding in an option to allow for either using a fixed baseline OR using a new baseline and future for each of the potential GAM models

# Load required libraries
library(mgcv)
library(tidyverse)
library(mvtnorm) 
library(Smisc)
library(rgeos)
library(akima)
library(sf)
library(viridis)
library(cowplot)
library(corrr)
library(maptools)
library(ROCR)
library(hydroGOF)
library(sdm)
library(ggrepel)
library(ggmap)
library(raster)
library(tools)
library(randomForest)
library(broom)
library(snakecase)
library(here)

out.dir<- here::here("Results")
if(!file.exists(out.dir)){
  dir.create(out.dir)
}

# Core functions: log likelihood, posterior, bayes  ----------------------------------------
# NegativeLL function -- the likelihood of NEVA votes given the base and future, which comes out of a potential gam model. In this formulaion, we propose that the NEVA acts on the biomass model component. Not the presence component. 
loglike_bio_func<- function(gam.mod0.p, gam.mod0.b, season, cand.params.b, base.preds, fut.preds, nevaD, nevaS, nevaE, dir.brks, fix.params = NULL){
  
  if(FALSE){
    gam.mod0.p = dat.use$mod.fitted.p[[1]]
    gam.mod0.b = dat.use$mod.fitted.b[[1]]
    season = season.use
    cand.params.b = cand.params.b
    base.preds = base.preds
    fut.preds = fut.preds
    nevaD = round(c(dir.dat$Negative, dir.dat$Neutral, dir.dat$Positive)*dir.wt, 0)
    nevaS = round(c(sens.dat$Low, sens.dat$Moderate, sens.dat$High, sens.dat$Very.High)*sens.wt, 0)
    nevaE = NULL
    dir.brks<- c(0.33, 0.67)
    fix.params<- fix.params
  }
  
  ilink.p<- family(gam.mod0.p)$linkinv
  ilink.b<- family(gam.mod0.b)$linkinv
  
  # NEVA conversion components
  # New baseline and future for specific candidate parameters
  base.preds.use<- base.preds$Data[[match(season, base.preds$SEASON)]]
  fut.preds.use<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
  
  # Get predictions for presence model component -- no sampling
  base.p<- predict.gam(gam.mod0.p, newdata = base.preds.use, type = "link")
  fut.p<- predict.gam(gam.mod0.p, newdata = fut.preds.use, type = "link")
  
  # Get predictions for biomass model component, which does involve sampling
  lpmat.base.b<- predict.gam(gam.mod0.b, newdata = base.preds.use, type = "lpmatrix")
  base.bio<- as.numeric(lpmat.base.b %*% t(cand.params.b))
  lpmat.fut.b<- predict.gam(gam.mod0.b, newdata = fut.preds.use, type = "lpmatrix")
  fut.bio<- as.numeric(lpmat.fut.b %*% t(cand.params.b))
  
  # Combined predictions
  base.c<- as.numeric(base.p * base.bio)
  fut.c<- as.numeric(fut.p * fut.bio)
  #base.c<- as.numeric(base.bio)
  #fut.c<- as.numeric(fut.bio)
  
  # Likelihood of the GAM biomass model...uses normal distribution
  # For proposed GAM, get new baseline and future statistics
  bmn<- mean(base.c, na.rm = T)
  bsd<- sd(base.c, na.rm = T)
  fmn<- mean(fut.c, na.rm = T)
  fsd<- sd(fut.c, na.rm = T)
  
  # Calculate directional effect and vulnerability (sensitivity and exposure) breaks
  Dbrks<- quantile(na.omit(base.c), prob = dir.brks)
  Vbrks<- quantile(na.omit(base.c), prob = seq(from = 0, to = 1, length.out = 6))
  md<- length(Vbrks)/2+1
  
  ## Directional effect piece: likelihood of NEVA votes given our bmn, bsd, fmn, fsd
  if(!is.null(nevaD)){
    dtmp<- pnorm(Dbrks, fmn, fsd) # Cumulative distribution function, p(x) <= dbrks 0.33 quant and p(x) <= dbrks 0.67 quant from N(fmn, fsd) dist
    pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
    
    like.dir<- dmultinom(nevaD, prob = pd, log = TRUE) # Likelihood of observed directional effect NEVA votes given multinomial described by pd per bin probabilities
  }
  
  ## Vulnerability piece: sensitivity and exposure
  if(!is.null(nevaS)){
    sentmp<- pnorm(Vbrks, fmn, fsd) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
    sentmp<- c(sentmp[1], diff(sentmp), 1-sentmp[length(sentmp)])
    
    # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
    send<- rep(NA, md)
    send[1]<- sentmp[md]
    
    for(k in 1:(md-1)){
      send[k+1]<- sentmp[md+k]+sentmp[md-k]
    }
    
    like.sens<- try(dmultinom(nevaS, prob = send, log = TRUE), silent = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
  }
  
  if(!is.null(nevaE)){
    exptmp<- pnorm(Vbrks, fmn, fsd) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
    exptmp<- c(exptmp[1], diff(exptmp), 1-exptmp[length(exptmp)])
    
    # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
    expd<- rep(NA, md)
    expd[1]<- exptmp[md]
    
    for(k in 1:(md-1)){
      expd[k+1]<- exptmp[md+k]+exptmp[md-k]
    }
    
    like.exp<- try(dmultinom(nevaE, prob = expd, log = TRUE), silent = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
  }
  
  # Calculate and return relevant likelihood
  if(is.null(nevaD)){
    if(class(like.sens) == "try-error" | class(like.exp) == "try-error") {
      likelihood<- NA
    } else {
      likelihood<- like.sens + like.exp
    }
  } else if(is.null(nevaS) & is.null(nevaE)){
    if(class(like.dir) == "try-error"){
      likelihood<- NA
    } else {
      likelihood<- like.dir
    }
  } else if(!is.null(nevaS) & is.null(nevaE)){
    classes<- c(class(like.dir), class(like.sens))
    if(any(classes == "try-error")){
      likelihood<- NA
    } else {
      likelihood<- like.dir + like.sens
    }
  } else {
    classes<- c(class(like.dir), class(like.sens), class(like.exp))
    if(any(classes == "try-error")){
      likelihood<- NA
    } else {
      likelihood<- like.dir + like.sens + like.exp
    }
  }
  return(likelihood)
}

# Probability of the biomass model 
logprior_bio_func<- function(gam.mod0.b, cand.params.b){
  # Likelihood of the GAM baseline future (fitted curves) for presence component
  # Gather coefficiences and varcov matrix, needed to define the normal dist of each candidate parameter value
  coefs.mod.b<- coef(gam.mod0.b)
  vcov.mod.b<- vcov(gam.mod0.b)
  
  # Likelihood of these candidate parameter values, given coef (mean) and sigma (se)
  cand.params.logprior.b<- dmvnorm(cand.params.b, mean = coefs.mod.b, sigma = vcov.mod.b, log = TRUE)
  return(data.frame("Bio.Prior" =  cand.params.logprior.b))
}

# Posterior (loglike + logprior)
logposterior_bio_func<- function(gam.mod0.p, gam.mod0.b, season, cand.params.b, base.preds, fut.preds, nevaD, nevaS, nevaE, dir.brks, fix.params = NULL){
  # Likelihood of NEVA dir and vuln votes, given GAM baseline future (fitted curves)
  if(class(cand.params.b) == "numeric"){
    cand.params.b<- matrix(cand.params.b, nrow = 1, ncol = length(cand.params.b), byrow = T, dimnames = list(NULL, names(coef(gam.mod0.b))))
  }
  
  # Likelihood of NEVA votes given the biomass model
  loglike.out<- loglike_bio_func(gam.mod0.p, gam.mod0.b, season, cand.params.b, base.preds, fut.preds, nevaD, nevaS, nevaE, dir.brks, fix.params = NULL)
  
  # Probability of prior (gam biomass model)
  logprior.out<- logprior_bio_func(gam.mod0.b, cand.params.b)
  
  # Add em all together
  return(data.frame("Likelihood" = loglike.out, "Bio.Prior" = logprior.out$Bio.Prior, "Posterior" = loglike.out + logprior.out$Bio.Prior))
}

# Bayes Algorithm
sdm_neva_bayes_bio<- function(gam.mod0.p, gam.mod0.b, season, cand.params.b, base.preds, fut.preds, nevaD, nevaS, nevaE, dir.brks, fix.params = NULL){
  
  # Calculate likelihood, prior and posterior given candidate
  post<- logposterior_bio_func(gam.mod0.p, gam.mod0.b, season, cand.params.b, base.preds, fut.preds, nevaD, nevaS, nevaE, dir.brks, fix.params = NULL)
  likes<- cbind(post$Likelihood, post$Bio.Prior, post$Posterior)
  return(likes)
}

# Data prep - Alaways run -------------------------------------------------
# Function to get scaled temperatures
temp.scale<- function(new.temp, base.temp.mean, base.temp.sd){
  if(is.na(new.temp)){
    temp.scaled<- NA
    return(temp.scaled)
  } else {
    temp.scaled<- (new.temp - base.temp.mean)/base.temp.sd
    return(temp.scaled)
  }
}

## Loading in the data and some prep -- unzip any zipped folders and add just inner files to the Data folder. Zipping was done to reduce file size for hosting on GitHub.
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

## Testing dataframes and applying correct scale to match rescaled variables used in the model fitting process
dat.test.f<- dat %>%
  filter(., TRAIN.TEST == "TEST" & SEASON == "FALL") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM))) %>%
  left_join(., fall.rescale.df, by = "SEASON")
dat.test.f$DEPTH.Scale<- mapply(temp.scale, abs(dat.test.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
dat.test.f$SHELF_POS.Scale<- mapply(temp.scale, dat.test.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
dat.test.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, dat.test.f$SEASONALMU.OISST, fall.rescale.df$mean.t, fall.rescale.df$sd.t)

dat.test.f$PRESENCE.BIOMASS<- ifelse(dat.test.f$BIOMASS > 0, 1, 0)

## Testing dataframes
dat.test.s<- dat %>%
  filter(., TRAIN.TEST == "TEST" & SEASON == "SPRING") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM))) %>%
  left_join(., spring.rescale.df, by = "SEASON")
dat.test.s$DEPTH.Scale<- mapply(temp.scale, abs(dat.test.s$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
dat.test.s$SHELF_POS.Scale<- mapply(temp.scale, dat.test.s$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
dat.test.s$SEASONALMU.OISST.Scale<- mapply(temp.scale, dat.test.s$SEASONALMU.OISST, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

dat.test.s$PRESENCE.BIOMASS<- ifelse(dat.test.s$BIOMASS > 0, 1, 0)

# Create nested dataframes, one for testing, one for training
# Training
dat.train<- dat.train.f %>%
  bind_rows(., dat.train.s) %>%
  group_by(., COMNAME, SEASON) %>%
  nest(.key = "TRAIN.DATA") %>%
  arrange(COMNAME)

# Testing
dat.test<- dat.test.f %>%
  bind_rows(., dat.test.s) %>%
  group_by(., COMNAME, SEASON) %>%
  nest(.key = "TEST.DATA") %>%
  arrange(COMNAME) 

## We will also need correctly "rescaled" versions for the projection time periods
# Need base.preds, fut.preds, nevaD and nevaV
spring.preds<- readRDS(here::here("Data", "spring.rast.preds06122018.rds"))
fall.preds<- readRDS(here::here("Data", "fall.rast.preds06122018.rds"))

base.preds.sp<- spring.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("SPRING", nrow(.)))

base.preds.sp<- base.preds.sp %>%
  left_join(., spring.rescale.df, by = "SEASON")
base.preds.sp$DEPTH.Scale<- mapply(temp.scale, abs(base.preds.sp$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
base.preds.sp$SHELF_POS.Scale<- mapply(temp.scale, base.preds.sp$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
base.preds.sp$SEASONALMU.OISST.Scale<- mapply(temp.scale, base.preds.sp$Baseline, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

base.preds.sp<- base.preds.sp %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("FALL", nrow(.))) 
base.preds.f$DEPTH.Scale<- mapply(temp.scale, abs(base.preds.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
base.preds.f$SHELF_POS.Scale<- mapply(temp.scale, base.preds.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
base.preds.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, base.preds.f$Baseline, fall.rescale.df$mean.t, fall.rescale.df$sd.t)
base.preds.f<- base.preds.f %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds<- base.preds.f %>%
  bind_rows(., base.preds.sp)

## Future
fut.preds.sp<- spring.preds %>%
  dplyr::select(., x, y, Spring.2055.mu, Spring.2055.pct05, Spring.2055.pct95, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("SPRING", nrow(.))) %>%
  left_join(., spring.rescale.df, by = "SEASON")
fut.preds.sp$DEPTH.Scale<- mapply(temp.scale, abs(fut.preds.sp$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
fut.preds.sp$SHELF_POS.Scale<- mapply(temp.scale, fut.preds.sp$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
fut.preds.sp$SEASONALMU.2055.RCP85.OISST.Scale<- mapply(temp.scale, fut.preds.sp$Spring.2055.mu, spring.rescale.df$mean.t, spring.rescale.df$sd.t)
fut.preds.sp$SEASONAL05.2055.RCP85.OISST.Scale<- mapply(temp.scale, fut.preds.sp$Spring.2055.pct05, spring.rescale.df$mean.t, spring.rescale.df$sd.t)
fut.preds.sp$SEASONAL95.2055.RCP85.OISST.Scale<- mapply(temp.scale, fut.preds.sp$Spring.2055.pct95, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

fut.preds.sp<- fut.preds.sp %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Fall.2055.mu, Fall.2055.pct05, Fall.2055.pct95, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("FALL", nrow(.))) %>%
  left_join(., fall.rescale.df, by = "SEASON")
fut.preds.f$DEPTH.Scale<- mapply(temp.scale, abs(fut.preds.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
fut.preds.f$SHELF_POS.Scale<- mapply(temp.scale, fut.preds.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
fut.preds.f$SEASONALMU.2055.RCP85.OISST.Scale<- mapply(temp.scale, fut.preds.f$Fall.2055.mu, fall.rescale.df$mean.t, fall.rescale.df$sd.t)
fut.preds.f$SEASONAL05.2055.RCP85.OISST.Scale<- mapply(temp.scale, fut.preds.f$Fall.2055.pct05, fall.rescale.df$mean.t, fall.rescale.df$sd.t)
fut.preds.f$SEASONAL95.2055.RCP85.OISST.Scale<- mapply(temp.scale, fut.preds.f$Fall.2055.pct95, fall.rescale.df$mean.t, fall.rescale.df$sd.t)

fut.preds.f<- fut.preds.f %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds<- fut.preds.f %>%
  bind_rows(., fut.preds.sp)

## Also going to want these as prediction dataframes to plot GAM smooths
# Prediction dataframe
pred.dat.f<- with(base.preds$Data[[1]],
                  data.frame(SEASONALMU.OISST.Scale = c(seq(min(SEASONALMU.OISST.Scale, na.rm = T), max(SEASONALMU.OISST.Scale, na.rm = T), length = 500), rep(mean(SEASONALMU.OISST.Scale, na.rm = T), 1000)),
                             DEPTH.Scale = c(rep(mean(DEPTH.Scale, na.rm = T), 500), seq(min(DEPTH.Scale, na.rm = T), max(DEPTH.Scale, na.rm = T), length = 500), rep(mean(DEPTH.Scale, na.rm = T), 500)),
                             SHELF_POS.Scale = c(rep(mean(SHELF_POS.Scale, na.rm = T), 1000), seq(min(SHELF_POS.Scale, na.rm = T), max(SHELF_POS.Scale, na.rm = T), length = 500)), "SEASON" = rep("FALL", 1500)))
rescaled.dat.f<- data.frame("SST" = seq(min(dat.train.f$SEASONALMU.OISST, na.rm = T), max(dat.train.f$SEASONALMU.OISST, na.rm = T), length = 500),
                                 "Depth" = seq(min(abs(dat.train.f$DEPTH), na.rm = T), max(abs(dat.train.f$DEPTH), na.rm = T), length = 500),
                                 "Shelf_Pos" = seq(min(dat.train.f$SHELF_POS, na.rm = T), max(dat.train.f$SHELF_POS, na.rm = T), length = 500),
                                 "Season" = rep("FALL", 500))

pred.dat.s<- with(base.preds$Data[[2]],
                  data.frame(SEASONALMU.OISST.Scale = c(seq(min(SEASONALMU.OISST.Scale, na.rm = T), max(SEASONALMU.OISST.Scale, na.rm = T), length = 500), rep(mean(SEASONALMU.OISST.Scale, na.rm = T), 1000)),
                             DEPTH.Scale = c(rep(mean(DEPTH.Scale, na.rm = T), 500), seq(min(DEPTH.Scale, na.rm = T), max(DEPTH.Scale, na.rm = T), length = 500), rep(mean(DEPTH.Scale, na.rm = T), 500)),
                             SHELF_POS.Scale = c(rep(mean(SHELF_POS.Scale, na.rm = T), 1000), seq(min(SHELF_POS.Scale, na.rm = T), max(SHELF_POS.Scale, na.rm = T), length = 500)), "SEASON" = rep("SPRING", 1500)))
rescaled.dat.s<- data.frame("SST" = seq(min(dat.train.s$SEASONALMU.OISST, na.rm = T), max(dat.train.s$SEASONALMU.OISST, na.rm = T), length = 500),
                            "Depth" = seq(min(abs(dat.train.s$DEPTH), na.rm = T), max(abs(dat.train.s$DEPTH), na.rm = T), length = 500),
                            "Shelf_Pos" = seq(min(dat.train.s$SHELF_POS, na.rm = T), max(dat.train.s$SHELF_POS, na.rm = T), length = 500),
                            "Season" = rep("SPRING", 500))

pred.dat<- bind_rows(pred.dat.f, pred.dat.s)
rescaled.dat<- bind_rows(rescaled.dat.f, rescaled.dat.s)

# Bring in vulnerability assessment datasets
dir.dat<- read.csv(here::here("Data", "JHareDirectionalEffect.csv")) %>%
  mutate(., "COMNAME" = toupper(Species)) %>%
  dplyr::select(., -Species) %>%
  group_by(COMNAME) %>%
  nest(., .key = "Dir.Data")
vuln.dat<- read.csv(here::here("Data", "JHareQualitativeDataResults.csv")) %>%
  group_by(., Species, Attribute.Category) %>%
  summarise_at(., .vars = c("Low", "Moderate", "High", "Very.High"), mean) 
exp.dat<- vuln.dat %>%
  filter(., Attribute.Category == "Exposure.Factor") %>%
  mutate(., "COMNAME" = toupper(Species)) %>%
  data.frame(.) %>%
  dplyr::select(., -Species, -Attribute.Category) %>%
  group_by(COMNAME) %>%
  nest(., .key = "Exp.Data")
sens.dat<- vuln.dat %>%
  filter(., Attribute.Category == "Sensitivity.Attribute") %>%
  mutate(., "COMNAME" = toupper(Species)) %>%
  data.frame(.) %>%
  dplyr::select(., -Species, -Attribute.Category) %>% 
  group_by(COMNAME) %>%
  nest(., .key = "Sens.Data")
qual.dat<- dir.dat %>%
  left_join(., exp.dat, by = "COMNAME") %>%
  left_join(., sens.dat, by = "COMNAME")

# Bind to training data..
dat.full<- dat.train %>%
  left_join(., qual.dat, by = "COMNAME")

# Now to test data
dat.full<- dat.full %>%
  left_join(., dat.test, by = c("COMNAME", "SEASON"))

# Model fitting -----------------------------------------------------------
# Set up the loop
scenarios.type = "Dir.Sens"
mod.sims<- 1000
dir.brks<- c(0.33, 0.67)
out.dir<- out.dir
fix.params<- NULL
set.seed(13)

gam_fit_full_func<- function(df, response){
  if(response == "Presence"){
    gam.mod0<- gam(PRESENCE.BIOMASS ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = binomial(link = logit), select = TRUE)
    return(gam.mod0)
  } 
  
  if(response == "Biomass"){
    gam.mod0<- gam(BIOMASS.MOD ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = gaussian, select = TRUE)
    return(gam.mod0)
  }
}

gam_fit_red_func<- function(df, response){
  if(response == "Presence"){
    gam.mod0<- gam(PRESENCE.BIOMASS ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = binomial(link = logit), select = TRUE)
    return(gam.mod0)
  } 
  
  if(response == "Biomass"){
    gam.mod0<- gam(BIOMASS.MOD ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = gaussian, select = TRUE)
    return(gam.mod0)
  }
}

resids_map_func<- function(test.data, response, predicted, type) {
  # Preliminary spatial stuff
  # Spatial projections
  proj.wgs84<- "+init=epsg:4326" #WGS84
  proj.utm<- "+init=epsg:2960" #UTM 19
  
  # NELME domaine
  nelme<- st_read(here::here("Data", "NELME_clipped.shp"))
  st_crs(nelme)<- proj.wgs84
  
  #Bounds
  xlim.use<- c(-77, -65)
  ylim.use<- c(35, 45)
  
  states <- c("Maine", "New Hampshire", "Massachusetts", "Vermont", "New York", "Rhode Island", "Connecticut", "Delaware", "New Jersey", "Maryland", "Pennsylvania", "Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", "District of Columbia", "West Virgina")
  provinces <- c("Ontario", "Québec", "Nova Scotia", "New Brunswick")
  
  us <- raster::getData("GADM",country="USA",level=1)
  us.states <- us[us$NAME_1 %in% states,]
  us.states <- gSimplify(us.states, tol = 0.025, topologyPreserve = TRUE)
  canada <- raster::getData("GADM",country="CAN",level=1)
  ca.provinces <- canada[canada$NAME_1 %in% provinces,]
  ca.provinces <- gSimplify(ca.provinces, tol = 0.025, topologyPreserve = TRUE)
  
  us.states.f<- fortify(us.states, NAME_1)
  ca.provinces.f<- fortify(ca.provinces, NAME_1)
  
  if(response == "Presence" & type == "response"){
    resids<- predicted - test.data$PRESENCE.BIOMASS
    dat.resids<- data.frame("year" = test.data$EST_YEAR, "x" = test.data$DECDEG_BEGLON, "y" = test.data$DECDEG_BEGLAT, "resid" = resids)
    
    plots.out<- vector("list", length(unique(dat.resids$year)))
    
    zlim.use<- c(min(dat.resids$resid, na.rm = T), max(dat.resids$resid, na.rm = T))
    
    for(k in 1:length(unique(dat.resids$year))){
      data.use<- dat.resids[dat.resids$year == unique(dat.resids$year)[k],]
      pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$resid))
      pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
                              xo=seq(-87.99457, -57.4307, length = 115),
                              yo=seq(22.27352, 48.11657, length = 133))
      pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
      pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
      
      # Clip to nelme
      pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
      coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
      row.names(coords.keep)<- NULL
      pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
      names(pred.df.use)<- c("X", "Y", "z")
      
      # Plot
      plots.out[[k]]<- ggplot() + 
        geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
        scale_fill_gradient2(name = paste(unique(dat.resids$year)[k], " pred - obs", sep = ""), low = "blue", high = "red", mid = "white", midpoint = 0.0, na.value = "white", limits = c(-1,1)) +
        geom_map(data = us.states.f, map = us.states.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        geom_map(data = ca.provinces.f, map = ca.provinces.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        ylim(ylim.use) + ylab("Lat") +
        scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
        coord_fixed(1.3) + 
        theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=10), legend.title=element_text(size=10), plot.margin = unit(c(0, 0, 0, 0), "in"))
    }
    
    out<- plot_grid(plots.out[[1]], plots.out[[2]], plots.out[[3]], plots.out[[4]], plots.out[[5]], nrow = 2, ncol = 3, scale = 1)
  }
  
  if(response == "Biomass" & type == "response"){
    resids<- predicted - test.data$BIOMASS
    dat.resids<- data.frame("year" = test.data$EST_YEAR, "x" = test.data$DECDEG_BEGLON, "y" = test.data$DECDEG_BEGLAT, "resid" = resids)
    
    plots.out<- vector("list", length(unique(dat.resids$year)))
    
    zlim.use<- c(min(dat.resids$resid, na.rm = T), max(dat.resids$resid, na.rm = T))
    
    for(k in 1:length(unique(dat.resids$year))){
      data.use<- dat.resids[dat.resids$year == unique(dat.resids$year)[k],]
      pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$resid))
      pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
                              xo=seq(-87.99457, -57.4307, length = 115),
                              yo=seq(22.27352, 48.11657, length = 133))
      pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
      pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
      
      # Clip to nelme
      pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
      coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
      row.names(coords.keep)<- NULL
      pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
      names(pred.df.use)<- c("X", "Y", "z")
      
      # Plot
      plots.out[[k]]<- ggplot() + 
        geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
        scale_fill_gradient2(name = paste(unique(dat.resids$year)[k], " pred - obs", sep = ""), low = "blue", high = "red", mid = "white", midpoint = 0.0, na.value = "white") +
        geom_map(data = us.states.f, map = us.states.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        geom_map(data = ca.provinces.f, map = ca.provinces.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        ylim(ylim.use) + ylab("Lat") +
        scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
        coord_fixed(1.3) + 
        theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=10), legend.title=element_text(size=10), plot.margin = unit(c(0, 0, 0, 0), "in"))
    }
    
    out<- plot_grid(plots.out[[1]], plots.out[[2]], plots.out[[3]], plots.out[[4]], plots.out[[5]], nrow = 2, ncol = 3, scale = 1)
  }
  
  return(out)
}

predict_func<- function(mod.fitted.p, mod.fitted.b, response, percentile, test.data) {
  temp<- dplyr::select(test.data, one_of(c("DEPTH.Scale", percentile)))
  test.data<- data.frame(na.omit(temp))
  out.p<- round(as.numeric(predict.gam(mod.fitted.p, newdata = test.data, type = "response", se.fit = TRUE)$fit), 3)
  
  if(response == "Biomass"){
    out.b<- exp(round(as.numeric(predict.gam(mod.fitted.b, newdata = test.data, type = "response", se.fit = TRUE)$fit), 3))
    out.c<- out.p * out.b
    return(out.c)
  } else {
    return(out.p)
  }
}

neva_predict_bio_func<- function(mod.fitted.p, mod.fitted.b, test.data, neva.best.fit.b = best.fit.mat.b){
  temp<- dplyr::select(test.data, one_of(c("DECDEG_BEGLON", "DECDEG_BEGLAT", "DEPTH.Scale", "SEASONALMU.OISST.Scale")))
  temp.data<- data.frame(na.omit(temp))
  
  # Presence prediction
  neva.pred.p<- predict.gam(mod.fitted.p, newdata = temp.data, type = "response")
  
  # Biomass prediction
  lpmat<- predict.gam(mod.fitted.b, newdata = temp.data, type = "lpmatrix")
  neva.pred.b<- as.numeric(lpmat %*% t(neva.best.fit.b))
  ilink<- family(mod.fitted.b)$linkinv
  neva.pred.b<- exp(ilink(neva.pred.b)) 
  return(neva.pred.p*neva.pred.b)
}

pred_ranges_func<- function(predicted){
  pred.ranges<- data.frame("Min.Pred" = min(predicted, na.rm = T), "Max.Pred" = max(predicted, na.rm = T), "Mean.Pred" = mean(predicted, na.rm = T))
  return(pred.ranges)
}

auc_func<- function(test.data, predicted) {
  temp<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS")))
  test.data<- data.frame(na.omit(temp))
  if(all(test.data$PRESENCE.BIOMASS == 0)){
    return(NA)
  } else {
    col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
    dat<- prediction(predictions = predicted, labels = test.data[,col.ind])
    return(performance(dat, measure = "auc")@y.values[[1]])
  }
}

rmse_func<- function(test.data, predicted, response) {
  if(response == "Presence"){
    test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS"))) 
    test.data<- data.frame(na.omit(test.data))
    col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
    if(all(test.data$PRESENCE.BIOMASS == 0)){
      return(NA)
    } else {
      rmse.out<- rmse(sim = as.numeric(predicted), obs = test.data$PRESENCE.BIOMASS)
      rmse.out<- rmse.out/sd(test.data[,col.ind], na.rm = TRUE)
      return(rmse.out)
    }
  }
  
  if(response == "Biomass"){
    test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) 
    test.data<- data.frame(na.omit(test.data))
    col.ind<- which(colnames(test.data) == "BIOMASS")
    if(all(test.data$BIOMASS == 0)){
      return(NA)
    } else {
      rmse.out<- rmse(sim = as.numeric(predicted), obs = test.data$BIOMASS)
      rmse.out<- rmse.out/sd(test.data[,col.ind], na.rm = TRUE)
      return(rmse.out)
    }
  }
}

calib_stat_func<- function(test.data, predicted){
  test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS"))) 
  test.data<- data.frame(na.omit(test.data))
  col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
  if(all(test.data$PRESENCE.BIOMASS == 0)){
    return(NA)
  } else {
    calib.stat<- round(calibration(x = test.data[,col.ind], p = predicted)@statistic, 3)
    return(calib.stat)
  }
}

calib_plot_func<- function(test.data, predicted){
  test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS"))) 
  test.data<- data.frame(na.omit(test.data))
  col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
  if(all(test.data$PRESENCE.BIOMASS == 0)){
    return(NA)
  } else {
    calib.all<- calibration(x = test.data[,col.ind], p = predicted)
    calib.df<- calib.all@calibration
    calib.plot<- ggplot() +
      geom_point(data = calib.df, aes(x = pedicted_pobability_bin, y = observed_poportion)) +
      xlim(c(0, 1)) + 
      ylim(c(0, 1)) +
      ylab("Proportion of Observed Occurrences") +
      xlab("Predicted Probability of Occurrence") +
      geom_abline(intercept = 0, linetype = "dashed")
    return(calib.plot)
  }
}

corr_coeff_func<- function(test.data, predicted, response){
  if(response == "Presence"){
    test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS"))) 
    test.data<- data.frame(na.omit(test.data))
    col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
    if(all(test.data$PRESENCE.BIOMASS == 0)){
      return(NA)
    } else {
      mean.obs<- mean(test.data[,col.ind])
      mean.mod<- mean(predicted, na.rm = TRUE)
      sd.obs<- sd(test.data[,col.ind])
      sd.mod<- sd(predicted, na.rm = TRUE)
      samps<- nrow(test.data)
      corr.coeff<- ((1/samps)*(sum((predicted - mean.mod)*(test.data[,col.ind] - mean.obs))))/(sd.obs*sd.mod)
      return(corr.coeff)
    }
  }
  
  
  if(response == "Biomass"){
    test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) 
    test.data<- data.frame(na.omit(test.data))
    col.ind<- which(colnames(test.data) == "BIOMASS")
    if(all(test.data$BIOMASS == 0)){
      return(NA)
    } else {
      mean.obs<- mean(test.data[,col.ind])
      mean.mod<- mean(predicted, na.rm = TRUE)
      sd.obs<- sd(test.data[,col.ind])
      sd.mod<- sd(predicted, na.rm = TRUE)
      samps<- nrow(test.data)
      corr.coeff<- ((1/samps)*(sum((predicted - mean.mod)*(test.data[,col.ind] - mean.obs))))/(sd.obs*sd.mod)
      return(corr.coeff)
    }
  }
}

sdm_bias_func<- function(test.data, predicted, response){
  if(response == "Presence"){
    test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS"))) 
    test.data<- data.frame(na.omit(test.data))
    col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
    if(all(test.data$PRESENCE.BIOMASS == 0)){
      return(NA)
    } else {
      bias<- sd(predicted, na.rm = TRUE)/sd(test.data[,col.ind], na.rm = TRUE)
      return(bias)
    }
  }
  
  if(response == "Biomass"){
    test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) 
    test.data<- data.frame(na.omit(test.data))
    col.ind<- which(colnames(test.data) == "BIOMASS")
    if(all(test.data$BIOMASS == 0)){
      return(NA)
    } else {
      bias<- sd(predicted)/sd(test.data[,col.ind])
      return(bias)
    }
  }
}

mod.results<- data.frame(dat.full[,c(1:2)])
mod.results$DevExp.P<- rep(NA, nrow(mod.results))
mod.results$DevExp.B<- rep(NA, nrow(mod.results))
mod.results$Temp.DevExp.P<- rep(NA, nrow(mod.results))
mod.results$Temp.DevExp.B<- rep(NA, nrow(mod.results))
mod.results$AUC.SDM<- rep(NA, nrow(mod.results))
mod.results$RMSE.SDM.P<- rep(NA, nrow(mod.results))
mod.results$RMSE.SDM.B<- rep(NA, nrow(mod.results))
mod.results$CorrCoeff.SDM.P<- rep(NA, nrow(mod.results))
mod.results$CorrCoeff.SDM.B<- rep(NA, nrow(mod.results))
mod.results$Bias.SDM.P<- rep(NA, nrow(mod.results))
mod.results$Bias.SDM.B<- rep(NA, nrow(mod.results))
mod.results$Calib.SDM<- rep(NA, nrow(mod.results))
mod.results$Pred.Min.SDM<- rep(NA, nrow(mod.results))
mod.results$Pred.Max.SDM<- rep(NA, nrow(mod.results))
mod.results$Pred.BaseRate.SDM<- rep(NA, nrow(mod.results))
mod.results$RMSE.NEVA.B<- rep(NA, nrow(mod.results))
mod.results$CorrCoeff.NEVA.B<- rep(NA, nrow(mod.results))
mod.results$Bias.NEVA.B<- rep(NA, nrow(mod.results))

for(i in 1:nrow(dat.train)){
  
  # Get the data and fit the model
  dat.use<- dat.full[i,]
  season.use<- as.character(dat.use$SEASON[[1]])
  
  if(is.null(dat.use$Dir.Data[[1]])){
    print("No Vuln Assessment")
    next
  }
  
  dat.use<- dat.use %>%
    mutate(., "mod.fitted.p" = pmap(list(df = TRAIN.DATA, response = list("Presence")), possibly(gam_fit_full_func, NA)),
           "mod.fitted.b" = pmap(list(df = TRAIN.DATA, response = list("Biomass")), possibly(gam_fit_full_func, NA)))
  gam.mod0.p<- dat.use$mod.fitted.p[[1]]
  
  # No SDM data..
  if(is.infinite(summary(gam.mod0.p)$dev.expl) | summary(gam.mod0.p)$dev.expl > 0.97 | is.na(dat.use$mod.fitted.b)){
    print("Bad SDM fit")
    next()
  }
  
  # Else keep going
  saveRDS(dat.use$mod.fitted.p[[1]], file = paste(out.dir, "gamfitpres", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), ".rds", sep = ""))
  saveRDS(dat.use$mod.fitted.b[[1]], file = paste(out.dir, "gamfitbio", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), ".rds", sep = ""))
  
  # Store deviance explained
  mod.results$DevExp.P[i]<- round(summary(dat.use$mod.fitted.p[[1]])$dev.expl, 3)
  mod.results$DevExp.B[i]<- round(summary(dat.use$mod.fitted.b[[1]])$dev.expl, 3)
  
  # Compare with reduced model...
  gam.p.red<- gam(PRESENCE.BIOMASS ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.use$TRAIN.DATA[[1]], family = binomial(link = logit), select = TRUE)
  gam.b.red<- gam(BIOMASS.MOD ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.use$TRAIN.DATA[[1]], family = gaussian, select = TRUE)
  dev.p.temp<- (summary(dat.use$mod.fitted.p[[1]])$dev.expl - summary(gam.p.red)$dev.expl)/summary(dat.use$mod.fitted.p[[1]])$dev.expl
  dev.b.temp<- (summary(dat.use$mod.fitted.b[[1]])$dev.expl - summary(gam.b.red)$dev.expl)/summary(dat.use$mod.fitted.b[[1]])$dev.expl
  mod.results$Temp.DevExp.P[i]<- round(dev.p.temp, 3)
  mod.results$Temp.DevExp.B[i]<- round(dev.b.temp, 3)
  
  # Predictions and residuals
  dat.use<- dat.use %>%
    mutate(., "Predicted.SDM.P" = pmap(list(mod.fitted.p = mod.fitted.p, mod.fitted.b = mod.fitted.b, percentile = "SEASONALMU.OISST.Scale", test.data = TEST.DATA, response = "Presence"), possibly(predict_func, NA)),
           "Predicted.SDM.B" = pmap(list(mod.fitted = mod.fitted.p, mod.fitted.b = mod.fitted.b, percentile = "SEASONALMU.OISST.Scale", test.data = TEST.DATA, response = "Biomass"), possibly(predict_func, NA)),
           "Residuals.SDM.P" = pmap(list(test.data = TEST.DATA, response = list("Presence"), predicted = Predicted.SDM.P, type = "response"), possibly(resids_map_func, NA)),
           "Residuals.SDM.B" = pmap(list(test.data = TEST.DATA, response = list("Biomass"), predicted = Predicted.SDM.B, type = "response"), possibly(resids_map_func, NA)),
           "Pred.Ranges.SDM" = map(Predicted.SDM.P, possibly(pred_ranges_func, NA)),
           "AUC.SDM" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.P), possibly(auc_func, NA)),
           "RMSE.SDM.P" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.P, response = "Presence"), possibly(rmse_func, NA)),
           "RMSE.SDM.B" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.B, response = "Biomass"), possibly(rmse_func, NA)),
           "Calib.Stat.SDM" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.P), possibly(calib_stat_func, NA)),
           "Calib.Plot.SDM" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.P), possibly(calib_plot_func, NA)),
           "CorrCoeff.SDM.P" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.P, response = "Presence"), possibly(corr_coeff_func, NA)),
           "CorrCoeff.SDM.B" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.B, response = "Biomass"), possibly(corr_coeff_func, NA)), 
           "Bias.SDM.P" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.P, response = "Presence"), possibly(sdm_bias_func, NA)),
           "Bias.SDM.B" = pmap(list(test.data = TEST.DATA, predicted = Predicted.SDM.B, response = "Biomass"), possibly(sdm_bias_func, NA)))
  
  mod.results$AUC.SDM[i]<- round(as.numeric(dat.use$AUC.SDM[[1]]), 3)
  mod.results$RMSE.SDM.P[i]<- round(as.numeric(dat.use$RMSE.SDM.P[[1]]), 3)
  mod.results$RMSE.SDM.B[i]<- round(as.numeric(dat.use$RMSE.SDM.B[[1]]), 3)
  mod.results$Calib.SDM[i]<- round(as.numeric(dat.use$Calib.Stat.SDM[[1]]), 3)
  mod.results$Pred.Min.SDM[i]<- round(as.numeric(dat.use$Pred.Ranges.SDM[[1]]$Min.Pred), 3)
  mod.results$Pred.Max.SDM[i]<- round(as.numeric(dat.use$Pred.Ranges.SDM[[1]]$Max.Pred), 3)
  mod.results$Pred.BaseRate.SDM[i]<- round(as.numeric(dat.use$Pred.Ranges.SDM[[1]]$Mean.Pred), 3)
  mod.results$CorrCoeff.SDM.P[i]<- round(as.numeric(dat.use$CorrCoeff.SDM.P[[1]]), 3)
  mod.results$CorrCoeff.SDM.B[i]<- round(as.numeric(dat.use$CorrCoeff.SDM.B[[1]]), 3)
  mod.results$Bias.SDM.P[i]<- round(as.numeric(dat.use$Bias.SDM.P[[1]]), 3)
  mod.results$Bias.SDM.B[i]<- round(as.numeric(dat.use$Bias.SDM.B[[1]]), 3)

  # A few plots --
  # Spatial residual maps
  resids.plot<- dat.use$Residuals.SDM.P[[1]]
  ggsave(paste(out.dir, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_PresResidsSDMMap.jpg", sep = ""), resids.plot, width = 11, height = 8, dpi = 125, units = "in")
  
  resids.plot<- dat.use$Residuals.SDM.B[[1]]
  ggsave(paste(out.dir, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_BioResidsSDMMap.jpg", sep = ""), resids.plot, width = 11, height = 8, dpi = 125, units = "in")
  
  # Validation plots
  calib.plot<- dat.use$Calib.Plot.SDM[[1]]
  if(any(is.na(calib.plot))){
    next
  }
  
  calib.plot<- calib.plot +
    annotate("text", x = 0.85, y = 0.15, label = paste("Pred ranges = c(", round(mod.results$Pred.Min.SDM[i], 2), ", ", round(mod.results$Pred.Max.SDM[i], 2), ")", sep = "")) +
    annotate("text", x = 0.85, y = 0.1, label = paste("Base rate = ", round(mod.results$Pred.BaseRate.SDM[i], 2))) +
    annotate("text", x = 0.85, y = 0.05, label = paste("AUC = ", round(mod.results$AUC.SDM[i], 2))) +
    annotate("text", x = 0.85, y = 0.0, label = paste("Calib = ", round(mod.results$Calib.SDM[i], 2)))
  ggsave(paste(out.dir, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_CalibSDMPlot.jpg", sep = ""), calib.plot, width = 11, height = 8, dpi = 125, units = "in")
  
  # Draws of candidate models from the posterior for biomass component
  coef.mod0.b<- coef(dat.use$mod.fitted.b[[1]])
  vc.mod0.b<- vcov(dat.use$mod.fitted.b[[1]])
  cand.params.all.b<- rmvnorm(mod.sims, mean = coef.mod0.b, sigma = vc.mod0.b) 
  
  # Qualitative Data
  dir.dat<- dat.use$Dir.Data[[1]]
  sens.dat<- dat.use$Sens.Data[[1]]
  exp.dat<- dat.use$Exp.Data[[1]]
  
  # Adjustments
  sens.wt<- 5/5
  exp.wt<- 5/4
  dir.wt<- 5/3
  
  # Fix any of these?
  if(!is.null(fix.params)){
    if(fix.params == "All"){
      # Modify to match all column names of cand.params except intercept
      fix.params.use<- names(coef.mod0)[-which(names(coef.mod0) %in% "(Intercept)")]
      
      # Cand.params.all index...
      cand.params.all.ind<- which(colnames(cand.params.all) %in% fix.params.use, arr.ind = T)
      
      # Coef index
      coef.mod0.fix<- coef.mod0[names(coef.mod0) %in% fix.params.use]
      coef.mod0.fix.mat<- matrix(coef.mod0.fix, nrow = 1, ncol = length(coef.mod0.fix), byrow = T)
      
      # Move over values
      cand.params.all[,cand.params.all.ind]<- coef.mod0.fix.mat[rep(1:nrow(coef.mod0.fix.mat), times = nrow(cand.params.all)),]
      
      # Randomly generate intercept values...
      cand.params.all[,1]<- runif(mod.sims, -3, 3)
    } else {
      # Modify to match column names of cand.params
      fix.params.use<- paste(paste("s(", fix.params, ")", sep = ""), ".", seq(from = 1, to = 9, by = 1), sep = "")
      
      # Cand.params.all index...
      cand.params.all.ind<- which(colnames(cand.params.all) %in% fix.params.use, arr.ind = T)
      
      # Coef index
      coef.mod0.fix<- coef.mod0[names(coef.mod0) %in% fix.params.use]
      coef.mod0.fix.mat<- matrix(coef.mod0.fix, nrow = 1, ncol = length(coef.mod0.fix), byrow = T)
      
      # Move over values
      cand.params.all[,cand.params.all.ind]<- coef.mod0.fix.mat[rep(1:nrow(coef.mod0.fix.mat), times = nrow(cand.params.all)),]
    }
  }
  
  out.likes<- array(dim = c(mod.sims, 3))
  out.cands.b<- array(dim = c(mod.sims, ncol(cand.params.all.b)))
  
  for(k in 1:nrow(cand.params.all.b)){
    cand.params.b<- cand.params.all.b[k,]
    out.cands.b[k,]<- cand.params.b
    
    if(scenarios.type == "Dir") {
      out.likes[k,]<- sdm_neva_bayes_bio(gam.mod0.p = dat.use$mod.fitted.p[[1]], gam.mod0.b = dat.use$mod.fitted.b[[1]], season = season.use, cand.params.b = cand.params.b, base.preds = base.preds, fut.preds = fut.preds, nevaD = round(c(dir.dat$Negative, dir.dat$Neutral, dir.dat$Positive)*dir.wt, 0), nevaS = NULL, nevaE = NULL, dir.brks = dir.brks, fix.params = fix.params)
    } 
    
    if(scenarios.type == "Vuln"){
      out.likes[k,]<- sdm_neva_bayes_bio(gam.mod0.p = dat.use$mod.fitted.p[[1]], gam.mod0.b = dat.use$mod.fitted.b[[1]], season = season.use, cand.params.b = cand.params.b, base.preds = base.preds, fut.preds = fut.preds, nevaD = NULL, nevaS = round(c(sens.dat$Low, sens.dat$Moderate, sens.dat$High, sens.dat$Very.High)*sens.wt, 0), nevaE = round(c(exp.dat$Low, exp.dat$Moderate, exp.dat$High, exp.dat$Very.High)*exp.wt, 0), dir.brks = dir.brks, fix.params = fix.params)
    }
    
    if(scenarios.type == "Both"){
      out.likes[k,]<- sdm_neva_bayes_bio(gam.mod0.p = dat.use$mod.fitted.p[[1]], gam.mod0.b = dat.use$mod.fitted.b[[1]], season = season.use, cand.params.b = cand.params.b, base.preds = base.preds, fut.preds = fut.preds, nevaD = round(c(dir.dat$Negative, dir.dat$Neutral, dir.dat$Positive)*dir.wt, 0), nevaS = round(c(sens.dat$Low, sens.dat$Moderate, sens.dat$High, sens.dat$Very.High)*sens.wt, 0), nevaE = round(c(exp.dat$Low, exp.dat$Moderate, exp.dat$High, exp.dat$Very.High)*exp.wt, 0), dir.brks = dir.brks, fix.params = fix.params)
    }
    
    if(scenarios.type == "Dir.Sens"){
      out.likes[k,]<- sdm_neva_bayes_bio(gam.mod0.p = dat.use$mod.fitted.p[[1]], gam.mod0.b = dat.use$mod.fitted.b[[1]], season = season.use, cand.params.b = cand.params.b, base.preds = base.preds, fut.preds = fut.preds, nevaD = round(c(dir.dat$Negative, dir.dat$Neutral, dir.dat$Positive)*dir.wt, 0), nevaS = round(c(sens.dat$Low, sens.dat$Moderate, sens.dat$High, sens.dat$Very.High)*sens.wt, 0), nevaE = NULL, dir.brks = dir.brks, fix.params = fix.params)
    }
  }
  
  # Bring together results in a list
  out.likes[is.infinite(out.likes)]<- NA
  mcmc.results<- list(out.likes, out.cands.b)
  
  # Processing and results plots
  likes.df<- data.frame(mcmc.results[[1]])
  likes.nsamps<- nrow(likes.df)
  
  posts.df.b<- data.frame(mcmc.results[[2]])
  posts.nsamps<- nrow(posts.df.b)
  
  # Likes dataframe and plots of mixing
  names(likes.df)<- c("Likelihood (NEVA|SDM)", "Prior P(SDM.B)", "Posterior")
  likes.df$Iteration<- rep(seq(from = 1, to = likes.nsamps, by = 1))
  likes.df<- likes.df %>%
    gather(., Sample, Value, -Iteration)
  likes.df$Sample<- factor(likes.df$Sample, levels = c("Likelihood (NEVA|SDM)", "Prior P(SDM.B)", "Posterior"))
  
  out.plot<- ggplot(likes.df) +
    geom_line(aes(x = Iteration, y = Value), alpha = 0.25) +
    #scale_color_manual(values = c('#4daf4a', '#377eb8', '#e41a1c')) + 
    ylab("Log Likelihood") + 
    facet_wrap(~Sample, scales = "free") + 
    theme_bw()
  ggsave(paste(out.dir, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_LikePriorPost.jpg", sep = ""), out.plot, width = 11, height = 8, dpi = 125, units = "in")
  
  # Predictions from best model 
  like.posts.b<- likes.df[likes.df$Sample == "Prior P(SDM.B)",]
  posts.samps.b<- posts.df.b
  
  # Get maximimum value and extract row from posts.samps
  mod.ind.b<- like.posts.b$Iteration[which.max(like.posts.b$Value)]
  best.fit.b<- posts.samps.b[mod.ind.b,]
  best.fit.mat.b<- matrix(as.numeric(best.fit.b), nrow = 1, ncol = length(best.fit.b), byrow = T, dimnames = list(NULL, names(coef(dat.use$mod.fitted.b))))
  
  # Make predictions with these values
  dat.use<- dat.use %>%
    mutate(., "Predicted.NEVA.B" = pmap(list(mod.fitted.p = mod.fitted.p, mod.fitted.b = mod.fitted.b, test.data = TEST.DATA, neva.best.fit.b = list(best.fit.mat.b)), possibly(neva_predict_bio_func, NA)),
           "Residuals.NEVA.B" = pmap(list(test.data = TEST.DATA, predicted = Predicted.NEVA.B, response = "Biomass", type = "response"), possibly(resids_map_func, NA)),
           "RMSE.NEVA.B" = pmap(list(test.data = TEST.DATA, predicted = Predicted.NEVA.B, response = "Biomass"), possibly(rmse_func, NA)),
           "CorrCoeff.NEVA.B" = pmap(list(test.data = TEST.DATA, predicted = Predicted.NEVA.B, response = "Biomass"), possibly(corr_coeff_func, NA)), 
           "Bias.NEVA.B" = pmap(list(test.data = TEST.DATA, predicted = Predicted.NEVA.B, response = "Biomass"), possibly(sdm_bias_func, NA)))
  
  # Store key results
  mod.results$RMSE.NEVA.B[i]<- round(as.numeric(dat.use$RMSE.NEVA.B[[1]]), 3)
  mod.results$CorrCoeff.NEVA.B[i]<- round(as.numeric(dat.use$CorrCoeff.NEVA.B[[1]]), 3)
  mod.results$Bias.NEVA.B[i]<- round(as.numeric(dat.use$Bias.NEVA.B[[1]]), 3)
  
  # A few plots --
  # Spatial residual maps
  resids.plot<- dat.use$Residuals.NEVA.B[[1]]
  ggsave(paste(out.dir, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_BioResidsNEVAMap.jpg", sep = ""), resids.plot, width = 11, height = 8, dpi = 125, units = "in")
  
  # Prediction Maps
  source(here("Code", "sdm_predictionmaps_func.R"))
  maps<- plot_func(response = "Presence", mod.fitted.p = dat.use$mod.fitted.p[[1]], mod.fitted.b = dat.use$mod.fitted.b[[1]], season = season.use, base.preds = base.preds, fut.preds = fut.preds, like.posts = like.posts.b, posts.samps = posts.df.b)
  ggsave(paste(out.dir, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_PresenceMaps.jpg", sep = ""), maps, width = 10, height = 11, dpi = 125, units = "in")
  rm(maps)
  maps<- plot_func(response = "Biomass", mod.fitted.p = dat.use$mod.fitted.p[[1]], mod.fitted.b = dat.use$mod.fitted.b[[1]], season = season.use, base.preds = base.preds, fut.preds = fut.preds, like.posts = like.posts.b, posts.samps = posts.df.b)
  ggsave(paste(out.dir, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_BiomassMaps.jpg", sep = ""), maps, width = 10, height = 11, dpi = 125, units = "in")
  rm(maps)
  
  # All possible prediction curves, the original smooth, and then the one we have selected...
  # Some prep
  ilink<- family(dat.use$mod.fitted.b[[1]])$linkinv

  pred.dat.use<- pred.dat[pred.dat$SEASON == season.use, ]
  rescaled.dat.use<- rescaled.dat[rescaled.dat$Season == season.use, ]
  
  # Predictor matrix
  pred.mat<- predict(dat.use$mod.fitted.b[[1]], newdata = pred.dat.use, type = "lpmatrix")

  # Original predictions
  pred0<- exp(ilink(pred.mat %*% coef(dat.use$mod.fitted.b[[1]])))
  
  # Best prediction
  like.posts = likes.df[likes.df$Sample == "Posterior",]
  mod.ind<- like.posts$Iteration[which.max(like.posts$Value)]
  pred.best.vec<- as.numeric(posts.df.b[mod.ind,])
  names(pred.best.vec)<- names(coef(dat.use$mod.fitted.b[[1]]))
  pred.best<- exp(ilink(pred.mat %*% pred.best.vec))
  
  # All predictions
  pred.all<- apply(cand.params.all.b, 1, function(x) exp(ilink(pred.mat %*% x)))
  
  # Plotting
  # SST
  want<- 1:500
  sst.all<- data.frame(pred.all[want,])
  colnames(sst.all)<- paste("Mod.", seq(from = 1, to = ncol(sst.all)), sep = "")
  sst.all<- sst.all %>%
    gather(., "Model", "Pred")
  sst.all$Value<- rep(rescaled.dat.use$SST, length(unique(sst.all$Model)))
  sst.all$Parameter<- rep("SST", nrow(sst.all))
  
  sst.base<- ggplot(sst.all, aes(x = Value, y = Pred, group = Model)) + 
    geom_line(show.legend = F, color = "#bdbdbd") 
  
  sst.dat<- data.frame("Parameter" = rep("SST", 1000), "Value" = rep(rescaled.dat.use$SST, 2), "Pred" = c(pred0[want], pred.best[want]), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
  
  sst.out<- sst.base +
    geom_line(data = sst.dat, aes(x = Value, y = Pred, group = Model, color = Model)) +
    scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
    xlab("SST") +
    theme_bw() 
  
  # Depth
  want<- 501:1000
  depth.all<- data.frame(pred.all[want,])
  colnames(depth.all)<- paste("Mod.", seq(from = 1, to = ncol(depth.all)), sep = "")
  depth.all<- depth.all %>%
    gather(., "Model", "Pred")
  depth.all$Value<- rep(rescaled.dat.use$Depth, length(unique(depth.all$Model)))
  depth.all$Parameter<- rep("Depth", nrow(depth.all))
  
  depth.base<- ggplot(depth.all, aes(x = Value, y = Pred, group = Model)) + 
    geom_line(show.legend = F, color = "#bdbdbd")
  
  depth.dat<- data.frame("Parameter" = rep("Depth", 1000), "Value" = rep(rescaled.dat.use$Depth, 2), "Pred" = c(pred0[want], pred.best[want]), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
  
  depth.out<- depth.base +
    geom_line(data = depth.dat, aes(x = Value, y = Pred, group = Model, color = Model)) +
    scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
    xlab("Depth") +
    theme_bw() 
  
  out<- plot_grid(sst.out + theme(legend.position="none"), depth.out + theme(legend.position="none"), nrow = 1, ncol = 2, align = "hv", scale = 1)
  legend<- get_legend(sst.out)
  out<- plot_grid(out, legend, rel_widths = c(3, 0.5))
  ggsave(paste(out.dir, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_PredCurv.jpg", sep = ""), out, width = 11, height = 8, dpi = 125, units = "in")
  
  # Update
  print(paste(tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), " is done!", sep = ""))
  
  # Save MCMC results
  saveRDS(mcmc.results, file = paste(out.dir, "mcmc_", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), ".rds", sep = ""))
}

write.csv(mod.results, file = paste(out.dir, "mod.results.csv", sep = ""))


# Model predictions -------------------------------------------------------
res.files<- list.files(out.dir, "mcmc_")
mod.results<- read.csv(paste(out.dir, "mod.results.csv", sep = ""))

for(i in seq_along(res.files)){
  spp<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][1])
  season<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][2])
  spp.season.match<- paste(spp, season, sep = ".")
  
  # Model fit -- presence and biomass SDM
  mod.fitted.p<- readRDS(paste(out.dir, gsub("mcmc_", "gamfitpres", res.files[i]), sep = ""))
  mod.fitted.b<- readRDS(paste(out.dir, gsub("mcmc_", "gamfitbio", res.files[i]), sep = ""))
  ilink<- family(mod.fitted.b)$linkinv
  gam.coef<- names(coef(mod.fitted.p))
  
  # Model fit -- NEVA
  likes.temp<- data.frame(read_rds(paste(out.dir, res.files[i], sep = ""))[[1]])
  names(likes.temp)<- c("Likelihood", "Prior", "Posterior")
  likes.temp$Iteration<- rep(seq(from = 1, to = nrow(likes.temp), by = 1))
  mod.ind.b<- likes.temp$Iteration[which.max(likes.temp$Posterior)]
  
  posts.temp<- data.frame(read_rds(paste(out.dir, res.files[i], sep = ""))[[2]])
  best.fit<- posts.temp[mod.ind.b,]
  best.fit.mat<- matrix(as.numeric(best.fit), nrow = 1, ncol = length(best.fit), byrow = T, dimnames = list(NULL, gam.coef))
  
  # Make predictions
  # SDM
  sdm.base<- base.preds$Data[[match(season, base.preds$SEASON)]]
  sdm.base<- sdm.base %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))
  sdm.base.p<- round(predict.gam(mod.fitted.p, newdata = sdm.base, type = "response"), 2)
  sdm.base.b<- round(as.numeric(sdm.base.p) * exp(as.numeric(predict.gam(mod.fitted.b, newdata = sdm.base, type = "response"))), 2)
  
  # Mean climate model SST
  newdat.mu<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
  newdat.2055.rcp85.mu<- newdat.mu %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.2055.RCP85.OISST.Scale"))
  names(newdat.2055.rcp85.mu)[4]<- "SEASONALMU.OISST.Scale"
  sdm.2055.rcp85.mu.p<- round(predict.gam(mod.fitted.p, newdata = newdat.2055.rcp85.mu, type = "response"), 2)
  sdm.2055.rcp85.mu.b<- round(as.numeric(sdm.2055.rcp85.mu.p) * exp(as.numeric(predict.gam(mod.fitted.b, newdata = newdat.2055.rcp85.mu, type = "response"))), 2)
  
  # pct05 climate model
  newdat.05<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
  newdat.2055.rcp85.05<- newdat.05 %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL05.2055.RCP85.OISST.Scale"))
  names(newdat.2055.rcp85.05)[4]<- "SEASONALMU.OISST.Scale"
  sdm.2055.rcp85.pct05.p<- round(predict.gam(mod.fitted.p, newdata = newdat.2055.rcp85.05, type = "response"), 2)
  sdm.2055.rcp85.pct05.b<- round(as.numeric(sdm.2055.rcp85.pct05.p) * exp(as.numeric(predict.gam(mod.fitted.b, newdata = newdat.2055.rcp85.05, type = "response"))), 2)
  
  # pct95 climate model
  newdat.95<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
  newdat.2055.rcp85.95<- newdat.95 %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL95.2055.RCP85.OISST.Scale"))
  names(newdat.2055.rcp85.95)[4]<- "SEASONALMU.OISST.Scale"
  sdm.2055.rcp85.pct95.p<- round(predict.gam(mod.fitted.p, newdata = newdat.2055.rcp85.95, type = "response"), 2)
  sdm.2055.rcp85.pct95.b<- round(as.numeric(sdm.2055.rcp85.pct95.p) * exp(as.numeric(predict.gam(mod.fitted.b, newdata = newdat.2055.rcp85.95, type = "response"))), 2)
  
  # Presence only dataframes
  sdm.map.base.p<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y, "pred" = sdm.base.p)
  sdm.map.2055.rcp85.mu.p<- data.frame("x" = newdat.2055.rcp85.mu$x, "y" = newdat.2055.rcp85.mu$y, "pred" = sdm.2055.rcp85.mu.p)
  sdm.map.2055.rcp85.pct05.p<- data.frame("x" = newdat.2055.rcp85.05$x, "y" = newdat.2055.rcp85.05$y, "pred.05" = sdm.2055.rcp85.pct05.p)
  sdm.map.2055.rcp85.pct95.p<- data.frame("x" = newdat.2055.rcp85.95$x, "y" = newdat.2055.rcp85.95$y, "pred.95" = sdm.2055.rcp85.pct95.p)
  
  # Biomass
  sdm.map.base.b<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y, "pred" = sdm.base.b)
  sdm.map.2055.rcp85.mu.b<- data.frame("x" = newdat.2055.rcp85.mu$x, "y" = newdat.2055.rcp85.mu$y, "pred" = sdm.2055.rcp85.mu.b)
  sdm.map.2055.rcp85.pct05.b<- data.frame("x" = newdat.2055.rcp85.05$x, "y" = newdat.2055.rcp85.05$y, "pred.05" = sdm.2055.rcp85.pct05.b)
  sdm.map.2055.rcp85.pct95.b<- data.frame("x" = newdat.2055.rcp85.95$x, "y" = newdat.2055.rcp85.95$y, "pred.95" = sdm.2055.rcp85.pct95.b)

  # Differences
  sdm.diff.2055.rcp85.p<- data.frame("x" = sdm.map.2055.rcp85.mu.p$x, "y" = sdm.map.2055.rcp85.mu.p$y, "pred" = sdm.map.2055.rcp85.mu.p$pred - sdm.map.base.p$pred)
  sdm.lwr.diff.2055.rcp85.p<- data.frame("x" = sdm.map.2055.rcp85.pct05.p$x, "y" = sdm.map.2055.rcp85.pct05.p$y, "pred" = sdm.map.2055.rcp85.pct05.p$pred - sdm.map.base.p$pred)
  sdm.upr.diff.2055.rcp85.p<- data.frame("x" = sdm.map.2055.rcp85.pct95.p$x, "y" = sdm.map.2055.rcp85.pct95.p$y, "pred" = sdm.map.2055.rcp85.pct95.p$pred - sdm.map.base.p$pred)
 
  # Percent Differences
  sdm.percdiff.2055.rcp85.p<- data.frame("x" = sdm.map.2055.rcp85.mu.p$x, "y" = sdm.map.2055.rcp85.mu.p$y, "pred" = 100*((sdm.map.2055.rcp85.mu.p$pred - sdm.map.base.p$pred)/sdm.map.base.p$pred))
  sdm.lwr.percdiff.2055.rcp85.p<- data.frame("x" = sdm.map.2055.rcp85.pct05.p$x, "y" = sdm.map.2055.rcp85.pct05.p$y, "pred" = 100*((sdm.map.2055.rcp85.pct05.p$pred - sdm.map.base.p$pred)/sdm.map.base.p$pred))
  sdm.upr.percdiff.2055.rcp85.p<- data.frame("x" = sdm.map.2055.rcp85.pct95.p$x, "y" = sdm.map.2055.rcp85.pct95.p$y, "pred" = 100*((sdm.map.2055.rcp85.pct95.p$pred - sdm.map.base.p$pred)/sdm.map.base.p$pred))
 
  names(sdm.map.base.p)[3]<- "Baseline.sdm.p"
  names(sdm.map.2055.rcp85.mu.p)[3]<- "Future_2055_rcp85_mean.sdm.p"
  names(sdm.map.2055.rcp85.pct05.p)[3]<- "Future_2055_rcp85_pct05.sdm.p"
  names(sdm.map.2055.rcp85.pct95.p)[3]<- "Future_2055_rcp85_pct95.sdm.p"
  names(sdm.diff.2055.rcp85.p)[3]<- "Future_2055_rcp85_mean_diff.sdm.p"
  names(sdm.lwr.diff.2055.rcp85.p)[3]<- "Future_2055_rcp85_pct05_diff.sdm.p"
  names(sdm.upr.diff.2055.rcp85.p)[3]<- "Future_2055_rcp85_pct95_diff.sdm.p"
  names(sdm.percdiff.2055.rcp85.p)[3]<- "Future_2055_rcp85_mean_percdiff.sdm.p"
  names(sdm.lwr.percdiff.2055.rcp85.p)[3]<- "Future_2055_rcp85_pct05_percdiff.sdm.p"
  names(sdm.upr.percdiff.2055.rcp85.p)[3]<- "Future_2055_rcp85_pct95_percdiff.sdm.p"
  
 
  # Biomass
  sdm.diff.2055.rcp85.b<- data.frame("x" = sdm.map.2055.rcp85.mu.b$x, "y" = sdm.map.2055.rcp85.mu.b$y, "pred" = sdm.map.2055.rcp85.mu.b$pred - sdm.map.base.b$pred)
  sdm.lwr.diff.2055.rcp85.b<- data.frame("x" = sdm.map.2055.rcp85.pct05.b$x, "y" = sdm.map.2055.rcp85.pct05.b$y, "pred" = sdm.map.2055.rcp85.pct05.b$pred - sdm.map.base.b$pred)
  sdm.upr.diff.2055.rcp85.b<- data.frame("x" = sdm.map.2055.rcp85.pct95.b$x, "y" = sdm.map.2055.rcp85.pct95.b$y, "pred" = sdm.map.2055.rcp85.pct95.b$pred - sdm.map.base.b$pred)
 
  # Percent Differences
  sdm.percdiff.2055.rcp85.b<- data.frame("x" = sdm.map.2055.rcp85.mu.b$x, "y" = sdm.map.2055.rcp85.mu.b$y, "pred" = 100*((sdm.map.2055.rcp85.mu.b$pred - sdm.map.base.b$pred)/sdm.map.base.b$pred))
  sdm.lwr.percdiff.2055.rcp85.b<- data.frame("x" = sdm.map.2055.rcp85.pct05.b$x, "y" = sdm.map.2055.rcp85.pct05.b$y, "pred" = 100*((sdm.map.2055.rcp85.pct05.b$pred - sdm.map.base.b$pred)/sdm.map.base.b$pred))
  sdm.upr.percdiff.2055.rcp85.b<- data.frame("x" = sdm.map.2055.rcp85.pct95.b$x, "y" = sdm.map.2055.rcp85.pct95.b$y, "pred" = 100*((sdm.map.2055.rcp85.pct95.b$pred - sdm.map.base.b$pred)/sdm.map.base.b$pred))
 
  names(sdm.map.base.b)[3]<- "Baseline.sdm.b"
  names(sdm.map.2055.rcp85.mu.b)[3]<- "Future_2055_rcp85_mean.sdm.b"
  names(sdm.map.2055.rcp85.pct05.b)[3]<- "Future_2055_rcp85_pct05.sdm.b"
  names(sdm.map.2055.rcp85.pct95.b)[3]<- "Future_2055_rcp85_pct95.sdm.b"
  names(sdm.diff.2055.rcp85.b)[3]<- "Future_2055_rcp85_mean_diff.sdm.b"
  names(sdm.lwr.diff.2055.rcp85.b)[3]<- "Future_2055_rcp85_pct05_diff.sdm.b"
  names(sdm.upr.diff.2055.rcp85.b)[3]<- "Future_2055_rcp85_pct95_diff.sdm.b"
  names(sdm.percdiff.2055.rcp85.b)[3]<- "Future_2055_rcp85_mean_percdiff.sdm.b"
  names(sdm.lwr.percdiff.2055.rcp85.b)[3]<- "Future_2055_rcp85_pct05_percdiff.sdm.b"
  names(sdm.upr.percdiff.2055.rcp85.b)[3]<- "Future_2055_rcp85_pct95_percdiff.sdm.b"
  
  # Getting combined biomass predictions
  # Make predictions with these values
  lpmat.base<- predict.gam(mod.fitted.b, newdata = base.preds$Data[[match(season, base.preds$SEASON)]], type = "lpmatrix")
  combo.map.base<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y, "pred" = sdm.base.p * exp(ilink(as.numeric(lpmat.base %*% t(best.fit.mat)))))
  
  # Mean climate model SST
  newdat.mu<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
  newdat.2055.rcp85.mu<- newdat.mu %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONALMU.2055.RCP85.OISST.Scale"))
  names(newdat.2055.rcp85.mu)[4]<- "SEASONALMU.OISST.Scale"
  lpmat.2055.rcp85.mu<- predict.gam(mod.fitted.b, newdata = newdat.2055.rcp85.mu, type = "lpmatrix")
  combo.map.2055.rcp85.mu<- data.frame("x" = newdat.2055.rcp85.mu$x, "y" = newdat.2055.rcp85.mu$y, "pred" = round(sdm.2055.rcp85.mu.p * exp(ilink(as.numeric(lpmat.2055.rcp85.mu %*% t(best.fit.mat)))), 2))
  
  # pct05 climate model
  newdat.05<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
  newdat.2055.rcp85.05<- newdat.05 %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL05.2055.RCP85.OISST.Scale"))
  names(newdat.2055.rcp85.05)[4]<- "SEASONALMU.OISST.Scale"
  lpmat.2055.rcp85.pct05<- predict.gam(mod.fitted.b, newdata = newdat.2055.rcp85.05, type = "lpmatrix")
  combo.map.2055.rcp85.pct05<- data.frame("x" = newdat.2055.rcp85.05$x, "y" = newdat.2055.rcp85.05$y, "pred.05" = round(sdm.2055.rcp85.pct05.p * exp(ilink(as.numeric(lpmat.2055.rcp85.pct05 %*% t(best.fit.mat)))), 2))
  
  # pct95 climate model
  newdat.95<- fut.preds$Data[[match(season, fut.preds$SEASON)]] 
  newdat.2055.rcp85.95<- newdat.95 %>%
    unnest() %>%
    dplyr::select(., c("x", "y", "DEPTH.Scale", "SEASONAL95.2055.RCP85.OISST.Scale"))
  names(newdat.2055.rcp85.95)[4]<- "SEASONALMU.OISST.Scale"
  lpmat.2055.rcp85.pct95<- predict.gam(mod.fitted.b, newdata = newdat.2055.rcp85.95, type = "lpmatrix")
  combo.map.2055.rcp85.pct95<- data.frame("x" = newdat.2055.rcp85.95$x, "y" = newdat.2055.rcp85.95$y, "pred.95" = round(sdm.2055.rcp85.pct95.p * exp(ilink(as.numeric(lpmat.2055.rcp85.pct95 %*% t(best.fit.mat)))), 2))
  
  # Differences
  combo.diff.2055.rcp85<- data.frame("x" = combo.map.2055.rcp85.mu$x, "y" = combo.map.2055.rcp85.mu$y, "pred" = combo.map.2055.rcp85.mu$pred - combo.map.base$pred)
  combo.lwr.diff.2055.rcp85<- data.frame("x" = combo.map.2055.rcp85.pct05$x, "y" = combo.map.2055.rcp85.pct05$y, "pred" = combo.map.2055.rcp85.pct05$pred - combo.map.base$pred)
  combo.upr.diff.2055.rcp85<- data.frame("x" = combo.map.2055.rcp85.pct95$x, "y" = combo.map.2055.rcp85.pct95$y, "pred" = combo.map.2055.rcp85.pct95$pred - combo.map.base$pred)
  
  # Calculate fish availability for the different datasets: combo.map.base, combo.map.fut.mu, combo.map.fut.pct05, combo.map.fut.pct95, combo.diff, combo.lwr.diff, combo.upr.diff
  names(combo.map.base)[3]<- "Baseline.combo.b"
  names(combo.map.2055.rcp85.mu)[3]<- "Future_2055_rcp85_mean.combo.b"
  names(combo.map.2055.rcp85.pct05)[3]<- "Future_2055_rcp85_pct05.combo.b"
  names(combo.map.2055.rcp85.pct95)[3]<- "Future_2055_rcp85_pct95.combo.b"
  names(combo.diff.2055.rcp85)[3]<- "Future_mean_2055_rcp85_diff.combo.b"
  names(combo.lwr.diff.2055.rcp85)[3]<- "Future_pct05_2055_rcp85_diff.combo.b"
  names(combo.upr.diff.2055.rcp85)[3]<- "Future_pct95_2055_rcp85_diff.combo.b"
  names(combo.percdiff.2055.rcp85)[3]<- "Future_mean_2055_rcp85_percdiff.combo.b"
  names(combo.lwr.percdiff.2055.rcp85)[3]<- "Future_pct05_2055_rcp85_percdiff.combo.b"
  names(combo.upr.percdiff.2055.rcp85)[3]<- "Future_pct95_2055_rcp85_percdiff.combo.b"
  
  projections.dat<- sdm.map.base.p %>%
    left_join(., sdm.map.2055.rcp85.mu.p, by = c("x", "y")) %>%
    left_join(., sdm.map.2055.rcp85.pct05.p, by = c("x", "y")) %>%
    left_join(., sdm.map.2055.rcp85.pct95.p, by = c("x", "y")) %>%
    left_join(., sdm.diff.2055.rcp85.p, by = c("x", "y")) %>%
    left_join(., sdm.lwr.diff.2055.rcp85.p, by = c("x", "y")) %>%
    left_join(., sdm.upr.diff.2055.rcp85.p, by = c("x", "y")) %>%
    left_join(., sdm.percdiff.2055.rcp85.p, by = c("x", "y")) %>%
    left_join(., sdm.lwr.percdiff.2055.rcp85.p, by = c("x", "y")) %>%
    left_join(., sdm.upr.percdiff.2055.rcp85.p, by = c("x", "y")) %>%
    left_join(., sdm.map.base.b, by = c("x", "y")) %>%
    left_join(., sdm.map.2055.rcp85.mu.b, by = c("x", "y")) %>%
    left_join(., sdm.map.2055.rcp85.pct05.b, by = c("x", "y")) %>%
    left_join(., sdm.map.2055.rcp85.pct95.b, by = c("x", "y")) %>%
    left_join(., sdm.diff.2055.rcp85.b, by = c("x", "y")) %>%
    left_join(., sdm.lwr.diff.2055.rcp85.b, by = c("x", "y")) %>%
    left_join(., sdm.upr.diff.2055.rcp85.b, by = c("x", "y")) %>%
    left_join(., sdm.percdiff.2055.rcp85.b, by = c("x", "y")) %>%
    left_join(., sdm.lwr.percdiff.2055.rcp85.b, by = c("x", "y")) %>%
    left_join(., sdm.upr.percdiff.2055.rcp85.b, by = c("x", "y")) %>%
    left_join(., combo.map.base, by = c("x", "y")) %>%
    left_join(., combo.map.2055.rcp85.mu, by = c("x", "y")) %>%
    left_join(., combo.map.2055.rcp85.pct05, by = c("x", "y")) %>%
    left_join(., combo.map.2055.rcp85.pct95, by = c("x", "y")) %>%
    left_join(., combo.diff.2055.rcp85, by = c("x", "y")) %>%
    left_join(., combo.lwr.diff.2055.rcp85, by = c("x", "y")) %>%
    left_join(., combo.upr.diff.2055.rcp85, by = c("x", "y")) %>%
    left_join(., combo.percdiff.2055.rcp85, by = c("x", "y")) %>%
    left_join(., combo.lwr.percdiff.2055.rcp85, by = c("x", "y")) %>%
    left_join(., combo.upr.percdiff.2055.rcp85, by = c("x", "y")) 
  projections.dat$COMNAME<- spp
  projections.dat$SEASON<- season
  projections.dat<- projections.dat %>%
    gather(., "Proj.Class", "Projection", -x, -y, -COMNAME, -SEASON) %>%
    group_by(., COMNAME, SEASON, Proj.Class) %>%
    nest(.key = "Projections")
    
  if(i == 1){
    result<- projections.dat
    print(paste(spp, season, "is done!", sep = " "))
  } else {
    result<- bind_rows(result, projections.dat)
    print(paste(spp, season, "is done!", sep = " "))
  }
}

saveRDS(result, file = paste(out.dir, "SDMPredictions.rds", sep = ""))

# Results — Filtering, always run -----------------------------------------------
# Real data
mod.res<- read_csv(paste(out.dir, "/mod.results.csv", sep = "")) %>%
  drop_na(RMSE.SDM.B, CorrCoeff.SDM.B, Bias.SDM.B, RMSE.NEVA.B, CorrCoeff.NEVA.B, Bias.NEVA.B)

# Exploring cut offs... AUC > 0.7 in both seasons
mod.spp.keep<- mod.res %>%
  filter(., AUC.SDM >= 0.7 & CorrCoeff.SDM.B >= 0) %>%
  group_by(., COMNAME) %>%
  summarize_at(vars(SEASON), n_distinct) %>%
  filter(., SEASON == 2)

mod.res<- mod.res %>%
  filter(., COMNAME %in% mod.spp.keep$COMNAME)

# Results — Taylor diagrams -----------------------------------------------
# Getting maxSD for plotting
maxsd<- max(mod.res$Bias.SDM.B, 1)
sd.r<- 1

# Empty plot first
# Creating empty plot first
plot.base<- ggplot() + 
  scale_x_continuous(name = "Standard deviation (normalized)", limits = c(0, maxsd), breaks = seq(from = 0, to = maxsd, by = 0.5)) +
  scale_y_continuous(name = "Standard deviation (normalized)", limits = c(0, maxsd), breaks = seq(from = 0, to = maxsd, by = 0.5)) +
  theme_classic()

# Coeff D rays 
grad.corr.lines = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
for(i in 1:length(grad.corr.lines)){
  x.vec<- c(0, maxsd*grad.corr.lines[i])
  y.vec<- c(0, maxsd*sqrt(1 - grad.corr.lines[i]^2))
  
  if(i ==1){
    coeffd.rays.df<- data.frame("Ray" = rep(1, length(x.vec)), "x" = x.vec, "y" = y.vec)
  } else {
    temp<- data.frame("Ray" = rep(i, length(x.vec)), "x" = x.vec, "y" = y.vec)
    coeffd.rays.df<- bind_rows(coeffd.rays.df, temp)
  }
}

# Add rays
plot.coeffd<- plot.base +
  geom_line(data = coeffd.rays.df, aes(x = x, y = y, group = Ray), lty = "longdash", col = "lightgray")

coeffd.labs<- coeffd.rays.df %>%
  group_by(Ray) %>%
  summarize(., 
            "x" = max(x, na.rm = TRUE), 
            "y" = max(y, na.rm = TRUE)) %>%
  data.frame()

coeffd.labs$Label<- grad.corr.lines

plot.coeffd<- plot.coeffd +
  geom_label(data = coeffd.labs, aes(x = x, y = y, label = Label), color = "gray", fill = "white", label.size = NA)

# SD arcs
# Need to add in SD arcs
sd.arcs<- seq(from = 0, to = maxsd, by = 0.5)

for(i in 1:length(sd.arcs)){
  x.vec<- sd.arcs[i]*cos(seq(0, pi/2, by = 0.03))
  y.vec<- sd.arcs[i]*sin(seq(0, pi/2, by = 0.03))
  
  if(i ==1){
    sd.arcs.df<- data.frame("Arc" = rep(sd.arcs[1], length(x.vec)), "x" = x.vec, "y" = y.vec)
  } else {
    temp<- data.frame("Arc" = rep(sd.arcs[i], length(x.vec)), "x" = x.vec, "y" = y.vec)
    sd.arcs.df<- bind_rows(sd.arcs.df, temp)
  }
}

# Add arcs to plot.base
plot.sd<- plot.coeffd +
  geom_line(data = sd.arcs.df, aes(x = x, y = y, group = Arc), lty = "dotted", color = "lightgray") 

# Now gamma? -- Standard deviation arcs around the reference point
gamma<- pretty(c(0, maxsd), n = 4)[-1]
gamma<- gamma[-length(gamma)]
labelpos<- seq(45, 70, length.out = length(gamma))

for(gindex in 1:length(gamma)) {
  xcurve <- cos(seq(0, pi, by = 0.03)) * gamma[gindex] + sd.r
  endcurve <- which(xcurve < 0)
  endcurve <- ifelse(length(endcurve), min(endcurve) - 1, 105)
  ycurve <- sin(seq(0, pi, by = 0.03)) * gamma[gindex]
  maxcurve <- xcurve * xcurve + ycurve * ycurve
  startcurve <- which(maxcurve > maxsd * maxsd)
  startcurve <- ifelse(length(startcurve), max(startcurve) + 1, 0)
  x.vec<- xcurve[startcurve:endcurve]
  y.vec<- ycurve[startcurve:endcurve]
  
  if(gindex ==1){
    gamma.df<- data.frame("Gamma" = rep(gamma[1], length(x.vec)), "x" = x.vec, "y" = y.vec)
  } else {
    temp<- data.frame("Gamma" = rep(gamma[gindex], length(x.vec)), "x" = x.vec, "y" = y.vec)
    gamma.df<- bind_rows(gamma.df, temp)
  }
}

gamma.df$Gamma<- factor(gamma.df$Gamma, levels = unique(gamma.df$Gamma))

# Add em
plot.gamma<- plot.sd +
  geom_line(data = gamma.df, aes(x = x, y = y, group = Gamma), lty = "solid", col = "lightgray")

# Label...
gamma.labs<- gamma.df %>%
  group_by(Gamma) %>%
  summarize("x" = mean(x, na.rm = TRUE), 
            "y" = median(y, na.rm = TRUE))

plot.gamma<- plot.gamma +
  geom_label(data = gamma.labs, aes(x = x, y = y, label = Gamma), color = "gray", fill = "white", label.size = NA)

# Add in reference point
plot.all<- plot.gamma +
  geom_point(aes(x = sd.r, y = 0), color = "black", size = 2.75)

# Add in species points
mod.results.td<- mod.res %>%
  mutate(., "TD.X" = Bias.SDM.B * CorrCoeff.SDM.B,
         "TD.Y" = Bias.SDM.B * sin(acos(CorrCoeff.SDM.B)))

# Join with functional groups
func.groups<- read.csv(here::here("Data", "JHareSppFunctionalGroup.csv"))
func.groups$COMNAME<- toupper(func.groups$COMNAME)

mod.results.td<- mod.results.td %>%
  left_join(., func.groups)
mod.results.td$COMNAME<- to_sentence_case(mod.results.td$COMNAME)

# Faceting isn't working and not sure why...loop?
func.groups<- c("Groundfish", "Pelagic", "Coastal", "Invertebrates", "Diadromous", "Elasmobranch")
plots.out.f<- vector("list", length(func.groups))
plots.out.s<- vector("list", length(func.groups))

# Colors stuff -- what is the highest number we will need?
colors.max<- mod.results.td %>%
  group_by(SEASON, Functional.Group) %>%
  summarize(., "Rows" = n())
colors.n<- max(colors.max$Rows)
colors<- distinctColorPalette(colors.n)

for(i in seq_along(func.groups)){
  group.use<- func.groups[i]
  dat.use.f<- mod.results.td %>%
    filter(., SEASON == "FALL" & as.character(Functional.Group) == group.use)
  dat.use.s<- mod.results.td %>%
    filter(., SEASON == "SPRING" & as.character(Functional.Group) == group.use)
  
  if(length(unique(dat.use.f$COMNAME)) == length(colors)){
    colors.use<- colors
  } else {
    colors.use<- sample(colors, length(unique(dat.use.f$COMNAME)))
  }
 
  plots.out.f[[i]]<- plot.all +
    geom_point(data = dat.use.f, aes(x = TD.X, y = TD.Y, fill = COMNAME), pch = 21, alpha = 0.85, size = 2.75) +
    scale_fill_manual(name = "Species", values = colors.use) +
    geom_text(aes(label = "Correlation coefficient", x = 0.75, y = 0.75), angle = -45)
  plots.out.s[[i]]<- plot.all +
    geom_point(data = dat.use.s, aes(x = TD.X, y = TD.Y, fill = COMNAME), pch = 21, alpha = 0.85, size = 2.75) +
    scale_fill_manual(name = "Species", values = colors.use) +
    geom_text(aes(label = "Correlation coefficient", x = 0.75, y = 0.75), angle = -45)
  print(i)
}

# Arrange in a grid
out.f<- plot_grid(plots.out.f[[1]], plots.out.f[[2]], plots.out.f[[3]], plots.out.f[[4]], plots.out.f[[5]], plots.out.f[[6]], nrow = 2, labels = func.groups, label_x = 0.25, label_y = 1)
ggplot2::ggsave(filename = paste(out.dir, "FallTaylorDiagram", ".jpg", sep = ""), plot = out.f, width = 18, height = 10, units = "in")
out.s<- plot_grid(plots.out.s[[1]], plots.out.s[[2]], plots.out.s[[3]], plots.out.s[[4]], plots.out.s[[5]], plots.out.s[[6]], nrow = 2, labels = func.groups, label_x = 0.25, label_y = 1)
ggplot2::ggsave(filename = paste(out.dir, "SpringTaylorDiagram", ".jpg", sep = ""), plot = out.s, width = 18, height = 10, units = "in")

# Results — SDM shelfwide and regional changes ----------------------------
# Read in projections
results<- read_rds(paste(out.dir, "SDMPredictions.rds", sep = "")) # This should have everything we need. 

# Filter species 
dat.sub<- results %>%
  filter(., COMNAME %in% mod.res$COMNAME)
dat.full<- dat.sub

# Spatial projections
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19

# NELME, GoM and Southern New England-Mid Atlantic Bight Regions
nelme<- as(st_read(here::here("Data", "NELME_sf.shp"), "Spatial"))
gom<- as(st_read(here::here("Data", "GoM_sf.shp"), "Spatial"))
south<- as(st_read(here::here("Data", "SNEandMAB_sf.shp"), "Spatial"))

# For these changes, don't want to use the cell by cell differences or percent differences...
dat.sub<- dat.full[-which(grepl("diff", dat.full$Proj.Class)),] 

# Overlay func
overlay_func<- function(df, region, proj.use = proj4string(nelme)){
  dat.use<- data.frame(df)
  pts.temp<- dat.use
  coordinates(pts.temp)<- ~x+y
  proj4string(pts.temp)<- proj.use
  
  region.avg<- switch(region,
         NELME = mean(dat.use[,3], na.rm = T),
         GOM = mean(data.frame(pts.temp[!is.na(over(pts.temp, as(gom, "SpatialPolygons"))),])[,3], na.rm = T),
         South = mean(data.frame(pts.temp[!is.na(over(pts.temp, as(south, "SpatialPolygons"))),])[,3], na.rm = T))
  return(region.avg)
}

preds.df.sub<- dat.sub %>%
  mutate(., "NELME.Mean" = purrr::map2(Projections, list("NELME"), possibly(overlay_func, NA)),
         "GOM.Mean" = purrr::map2(Projections, list("GOM"), possibly(overlay_func, NA)),
         "South.Mean" = purrr::map2(Projections, list("South"), possibly(overlay_func, NA)))


# Now, we want regional differences (raw and percentages)...
preds.df.sub<- preds.df.sub %>%
  ungroup() %>%
  dplyr::select(., COMNAME, SEASON, Proj.Class, NELME.Mean, GOM.Mean, South.Mean) %>%
  gather(., Region, Projections, -COMNAME, -SEASON, -Proj.Class) %>%
  mutate(., Proj.ClassandRegion = paste(Proj.Class, Region, sep = "_")) %>%
  dplyr::select(., COMNAME, SEASON, Proj.ClassandRegion, Projections) %>%
  spread(., Proj.ClassandRegion, Projections) %>%
  mutate_if(., is.list, as.numeric) %>%
  mutate(., "NELME.2055.rcp85.Mean.Change.mu.b" = Future_mean.sdm.b_NELME.Mean - Baseline.sdm.b_NELME.Mean,
         "NELME.2055.rcp85.Mean.Change.pct95.b" = Future_pct95.sdm.b_NELME.Mean - Baseline.sdm.b_NELME.Mean,
         "NELME.2055.rcp85.Mean.Change.pct05.b" = Future_pct05.sdm.b_NELME.Mean - Baseline.sdm.b_NELME.Mean,
         "NELME.2055.rcp85.Mean.Perc.Change.mu.b" = 100*(NELME.2055.rcp85.Mean.Change.mu.b/Baseline.sdm.b_NELME.Mean),
         "NELME.2055.rcp85.Mean.Perc.Change.pct95.b" = 100*(NELME.2055.rcp85.Mean.Change.pct95.b/Baseline.sdm.b_NELME.Mean),
         "NELME.2055.rcp85.Mean.Perc.Change.pct05.b" = 100*(NELME.2055.rcp85.Mean.Change.pct05.b/Baseline.sdm.b_NELME.Mean),
         "GOM.2055.rcp85.Mean.Change.mu.b" = Future_mean.sdm.b_GOM.Mean - Baseline.sdm.b_GOM.Mean,
         "GOM.2055.rcp85.Mean.Change.pct95.b" = Future_pct95.sdm.b_GOM.Mean - Baseline.sdm.b_GOM.Mean,
         "GOM.2055.rcp85.Mean.Change.pct05.b" = Future_pct05.sdm.b_GOM.Mean - Baseline.sdm.b_GOM.Mean,
         "GOM.2055.rcp85.Mean.Perc.Change.mu.b" = 100*(GOM.2055.rcp85.Mean.Change.mu.b/Baseline.sdm.b_GOM.Mean),
         "GOM.2055.rcp85.Mean.Perc.Change.pct95.b" = 100*(GOM.2055.rcp85.Mean.Change.pct95.b/Baseline.sdm.b_GOM.Mean),
         "GOM.2055.rcp85.Mean.Perc.Change.pct05.b" = 100*(GOM.2055.rcp85.Mean.Change.pct05.b/Baseline.sdm.b_GOM.Mean),
         "South.2055.rcp85.Mean.Change.mu.b" = Future_mean.sdm.b_South.Mean - Baseline.sdm.b_South.Mean,
         "South.2055.rcp85.Mean.Change.pct95.b" = Future_pct95.sdm.b_South.Mean - Baseline.sdm.b_South.Mean,
         "South.2055.rcp85.Mean.Change.pct05.b" = Future_pct05.sdm.b_South.Mean - Baseline.sdm.b_South.Mean,
         "South.2055.rcp85.Mean.Perc.Change.mu.b" = 100*(South.2055.rcp85.Mean.Change.mu.b/Baseline.sdm.b_South.Mean),
         "South.2055.rcp85.Mean.Perc.Change.pct95.b" = 100*(South.2055.rcp85.Mean.Change.pct95.b/Baseline.sdm.b_South.Mean),
         "South.2055.rcp85.Mean.Perc.Change.pct05.b" = 100*(South.2055.rcp85.Mean.Change.pct05.b/Baseline.sdm.b_South.Mean))

preds.df.sub.perc.plot<- preds.df.sub %>%
  dplyr::select(., COMNAME, SEASON, NELME.2055.rcp85.Mean.Perc.Change.mu.b, NELME.2055.rcp85.Mean.Perc.Change.pct95.b, NELME.2055.rcp85.Mean.Perc.Change.pct05.b, GOM.2055.rcp85.Mean.Perc.Change.mu.b, GOM.2055.rcp85.Mean.Perc.Change.pct95.b, GOM.2055.rcp85.Mean.Perc.Change.pct05.b, South.2055.rcp85.Mean.Perc.Change.mu.b, South.2055.rcp85.Mean.Perc.Change.pct95.b, South.2055.rcp85.Mean.Perc.Change.pct05.b) %>%
  gather(., "Region_Scenario", "Change", -COMNAME, -SEASON)
preds.df.sub.perc.plot$Region_Only<- unlist(lapply(strsplit(preds.df.sub.perc.plot$Region_Scenario, "[.]"), "[", 1))
preds.df.sub.perc.plot$Scenario_Only<- unlist(lapply(strsplit(sub("[.]", "*", preds.df.sub.perc.plot$Region_Scenario), "[*]"), "[", 2))
preds.df.sub.perc.plot$Climate_Only<- ifelse(grepl("mu", preds.df.sub.perc.plot$Region_Scenario), "Mean", 
                                             ifelse(grepl("pct95", preds.df.sub.perc.plot$Region_Scenario), "95th percentile", "5th percentile"))

## Alright, we are now after a species - scenario - season - region - mean change dataframe...
res<- preds.df.sub.perc.plot 
res$Change[is.infinite(res$Change)]<- NA
res$Change[res$Change >= 500]<- 500

# Biomass shelfwide pct95 and pct05 ...
res.plot<- res

# Lets add a functional group column...
# Merge with functional groups....
func.groups<- read.csv(here:here("Data", "JHareSppFunctionalGroup.csv"))
func.groups$COMNAME<- toupper(func.groups$COMNAME)

res.plot<- res.plot %>%
  left_join(., func.groups, by = "COMNAME")
res.plot<- res.plot[!is.na(res.plot$Functional.Group),]
res.plot$Functional.Group<- factor(res.plot$Functional.Group, levels = c("Groundfish", "Pelagic", "Coastal", "Invertebrates", "Diadromous", "Elasmobranch"))
res.plot<- res.plot %>%
  dplyr::select(., -Region_Scenario, -Scenario_Only) %>%
  dplyr::arrange(., COMNAME, Functional.Group, SEASON, Region_Only, Climate_Only, Change)
res.plot$COMNAME<- factor(res.plot$COMNAME, levels = unique(res.plot$COMNAME))
res.plot$SEASON<- factor(res.plot$SEASON, levels = c("FALL", "SPRING"))
res.plot$Climate_Only<- factor(res.plot$Climate_Only, levels = c("Mean", "95th percentile", "5th percentile"))
res.plot$Region_Only<- factor(res.plot$Region_Only, levels = c("NELME", "GOM", "South"))

res.plot.nelme<- res.plot %>%
  dplyr::filter(., Region_Only == "NELME")
res.plot.gom<- res.plot %>%
  dplyr::filter(., Region_Only == "GOM")
res.plot.south<- res.plot %>%
  dplyr::filter(., Region_Only == "South")

res.plot.all<- list(res.plot.nelme, res.plot.gom, res.plot.south)
names(res.plot.all)<- c("NELME", "GoM", "South")

## Climate variability across the shelf by season
df<- data.frame(res.plot.all[[1]])
df.null<- cbind(expand.grid(COMNAME = levels(res.plot.nelme$COMNAME), SEASON = unique(res.plot.nelme$SEASON), Climate_Only = unique(res.plot.nelme$Climate_Only), Region_Only = levels(res.plot.nelme$Region_Only), Change = NA))
df.null<- df.null %>%
  left_join(., func.groups, by = "COMNAME")
df<- rbind(df[,], df.null)
df$duplicated<- paste(df$COMNAME, df$SEASON, df$Climate_Only, df$Functional.Group)
df<- df[!duplicated(df$duplicated),] %>%
  arrange(., COMNAME, SEASON)

for(j in seq_along(levels(df$SEASON))){
  dat.use<- df %>%
    dplyr::filter(., SEASON == levels(df$SEASON)[j])
  
  dodge <- position_dodge(width = 1)
  
  dat.use.df<- dat.use %>%
    tidyr::complete(COMNAME, SEASON)
  dat.use.df$COMNAME<- str_to_title(dat.use.df$COMNAME)
  dat.use.df$COMNAME<- factor(dat.use.df$COMNAME, levels = rev(unique(dat.use.df$COMNAME)))
  
  dat.use.df<- dat.use.df %>%
    drop_na(Change)
  dat.use.df$COMNAME.Plot<- factor(dat.use.df$COMNAME, levels = rev(unique(dat.use.df$COMNAME)), labels = rev(unique(to_sentence_case(as.character(dat.use.df$COMNAME)))))
  dat.use.df$Scenario<- ifelse(dat.use.df$Climate_Only == "5th percentile", "5th percentile",
                               ifelse(dat.use.df$Climate_Only == "95th percentile", "95th percentile", "Mean"))
  dat.use.df$Scenario<- factor(dat.use.df$Scenario, levels = c("5th percentile", "Mean", "95th percentile"))
  
  plot.out<- ggplot(data = dat.use.df, aes(x = COMNAME.Plot, y = Change, color = Scenario)) + 
    geom_hline(yintercept = 0, color = "#bdbdbd") +
    geom_point(alpha = 0.7, size = 2.5) +
    scale_color_manual("Climate Scenario", values = c("#3182bd", "#636363", "#de2d26")) +
    ylab("Percent change in relative biomass") + 
    xlab("Species") +
    ylim(c(-100, 500)) +
    theme_bw() +
    theme(text = element_text(size = 12),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black")) +
    coord_flip() +
    facet_wrap(~Functional.Group, scales = "free_y") +
    ggtitle(paste(names(res.plot.all)[1], levels(df$SEASON)[j], sep = " "))
  
  ggplot2::ggsave(filename = paste(out.dir, names(res.plot.all)[1], levels(df$SEASON)[j], ".jpg", sep = ""), plot = plot.out, width = 11, height = 8, units = "in")
}

### Season, differences by region
df<- data.frame(rbind(res.plot.all[[2]], res.plot.all[[3]]))
df<- df %>%
  filter(., as.character(Climate_Only) == "Mean")
df$Region_Only<- factor(df$Region_Only, levels = c("South", "GOM"), labels = c("Southern New England/Mid Atlantic Bight", "Gulf of Maine"))
df.null<- cbind(expand.grid(COMNAME = levels(df$COMNAME), SEASON = unique(df$SEASON), Climate_Only = unique(df$Climate_Only), Region_Only = levels(df$Region_Only), Change = NA))
df.null<- df.null %>%
  left_join(., func.groups, by = "COMNAME")
df<- rbind(df[,], df.null)
df$duplicated<- paste(df$COMNAME, df$SEASON, df$Region_Only, df$Functional.Group)
df<- df[!duplicated(df$duplicated),] %>%
  arrange(., COMNAME, SEASON) 
spp.keep<- df %>%
  group_by(COMNAME) %>%
  summarize_at(vars(Change), n_distinct, na.rm = T) %>%
  filter(., Change == 4) %>%
  dplyr::select(., COMNAME)
df<- df %>%
  filter(., as.character(COMNAME) %in% as.character(spp.keep$COMNAME))

# One plot per season, regional differences
seasons<- c("FALL", "SPRING")

for(i in seq_along(seasons)){
  dat.use<- df %>%
    filter(., SEASON == seasons[i])
  dodge <- position_dodge(width = 1)
  
  dat.use$COMNAME<- str_to_title(dat.use$COMNAME)
  dat.use$COMNAME<- factor(dat.use$COMNAME, levels = rev(unique(dat.use$COMNAME)))
  
  dat.use.df<- dat.use %>%
    drop_na(Change)
 
  dat.use.df$COMNAME.Plot<- factor(dat.use.df$COMNAME, levels = rev(unique(dat.use.df$COMNAME)), labels = rev(unique(to_sentence_case(as.character(dat.use.df$COMNAME)))))
  
  plot.means<- ggplot(data = dat.use.df, aes(x = COMNAME.Plot, y = Change, fill = Region_Only)) + 
    geom_bar(stat = "identity", width = 0.6, position = position_dodge(width = 0.6)) +
    scale_fill_manual(name = "Region", values  = c("#ff7f00", "#377eb8")) +
    ylab("Percent change in relative biomass") + 
    xlab("Species") +
    geom_hline(yintercept = 0) +
    theme_bw() +
    theme(text = element_text(size = 12),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black")) +
    coord_flip() +
    facet_wrap(~Functional.Group, scales = "free_y") +
    guides(fill = guide_legend(reverse = TRUE)) +
    ggtitle(paste(seasons[i], sep = " "))
  
  ggplot2::ggsave(filename = paste(out.dir, seasons[i], ".jpg", sep = ""), plot = plot.means, width = 11, height = 8, units = "in")
}

# Results - SDM vs. NEVA comparison ---------------------------------------
results<- read_rds(paste(out.dir, "SDMPredictions.rds", sep = "")) # This should have everything we need. 

# Filter species 
dat.sub<- results %>%
  filter(., COMNAME %in% mod.res$COMNAME)
dat.full<- dat.sub

# Functional groups
func.groups<- read.csv(here::here("Data", "JHareSppFunctionalGroup.csv"))
func.groups$COMNAME<- toupper(func.groups$COMNAME)

dat.full<- dat.full %>%
  left_join(., func.groups)

# Jon Certainty Data
jon.df<- read_csv(here::here("Data", "JHareQualitativeResults.csv"))

# Filter to species with Moderate, High, Very high certainty
cert.keep<- c("Moderate", "High", "Very high")
dat.full<- dat.full %>%
  left_join(., jon.df, by = "COMNAME") %>%
  filter(., SENSITIVITY.CERTAINTY %in% cert.keep & DIRECTIONAL.EFFECT.CERTAINTY %in% cert.keep)

# Crossman et al. 2012. A scalar sensitivity weight measuring the ratio of change in species distribution relative to the extent of a species distribution in the future. 
proj.class.keep<- c("Baseline.sdm.b", "Future_mean.sdm.b", "Future_pct95.sdm.b", "Future_pct05.sdm.b")

## Sensitivity
dat.sens<- dat.full %>%
  filter(., Proj.Class %in% proj.class.keep) %>%
  dplyr::select(., COMNAME, SEASON, Proj.Class, Projections) %>%
  spread(., Proj.Class, Projections)

sens_weight_func<- function(base, future){
  # Get individual cell absolutel changes in species distribution
  absdiff.df<- data.frame("x" = base$x, "y" = base$y, "AbsChange" = abs(base$Projection - future$Projection))
  ext<- sum(future$Projection, na.rm = TRUE)
  weight<- (sum(absdiff.df$AbsChange, na.rm = TRUE))/ext
  return(as.numeric(weight))
}

dat.sens<- dat.sens %>%
  mutate(., "Mean.SensitivityWeight" = as.numeric(pmap(list(base = Baseline.sdm.b, future = Future_mean.sdm.b), possibly(sens_weight_func, NA))),
         "Pct95.SensitivityWeight" = as.numeric(pmap(list(base = Baseline.sdm.b, future = Future_pct95.sdm.b), possibly(sens_weight_func, NA))),
         "Pct05.SensitivityWeight" = as.numeric(pmap(list(base = Baseline.sdm.b, future = Future_pct05.sdm.b), possibly(sens_weight_func, NA))))

dat.sens.keep<- dat.sens %>%
  dplyr::select(., COMNAME, SEASON, Mean.SensitivityWeight, Pct95.SensitivityWeight, Pct05.SensitivityWeight) %>%
  gather(., Impact.Scenario, Data, -COMNAME, -SEASON)

# Add in Jon's sensitivity data...
sens.cert.keep<- c("Moderate", "High", "Very high")
dat.sens.keep<- dat.sens.keep %>%
  left_join(., jon.df) %>%
  dplyr::filter(., SENSITIVITY.CERTAINTY %in% sens.cert.keep)

dat.sens.keep$SENSITIVITY<- factor(dat.sens.keep$SENSITIVITY, levels = c("Low", "Moderate", "High", "Very high"))
dat.sens.keep$SENSITIVITY.CERTAINTY<- factor(dat.sens.keep$SENSITIVITY.CERTAINTY, levels = c("Moderate", "High", "Very high"))

# Plot....
# First, slight change so that "Infinte" snesitivity weights, which means there is NO habitat left, receive the highest sensitivity weight
temp<- which(is.infinite(unlist(dat.sens.keep$Data)))
stats<- dat.sens.keep %>%
  filter(., !is.infinite(Data)) %>%
  group_by(SEASON, Impact.Scenario) %>%
  summarize(.,
            "SD" = sd(Data, na.rm = TRUE),
            "Max" = max(Data, na.rm = TRUE))

dat.sens.keep$Data[temp]<- 158+35.5

# Add standard deviation?
dat.sens.avg<- dat.sens.keep %>%
  group_by(., Impact.Scenario, SENSITIVITY) %>%
  summarise_at(., "Data", c("mean", "sd"), na.rm = T) %>%
  drop_na(SENSITIVITY)
dat.sens.avg$Impact.Scenario<- factor(dat.sens.avg$Impact.Scenario, levels = c("Pct05.SensitivityWeight", "Mean.SensitivityWeight", "Pct95.SensitivityWeight"))

plot.out.sens<- ggplot(data = dat.sens.avg, aes(x = SENSITIVITY, y = mean)) +
  geom_point() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  ylab("SDM Potential Impact") +
  xlab("NEVA Sensitivity Ranking") +
  facet_wrap(~Impact.Scenario, nrow = 3, scales = "free_y") +
  geom_smooth(data = dat.sens.avg, aes(x = as.numeric(SENSITIVITY), y = mean), method = "lm", se = FALSE) +
  theme_bw() 
ggsave(paste(out.dir, "Sensitivity_vs_Impact.jpg", sep = ""), plot.out.sens, width = 11, height = 8, units = "in")

## Directional effect -- could look at this both on raw change as well as percent change. Everything we need should be in the dat.full object
dat.dir<- dat.full %>%
  dplyr::select(., COMNAME, SEASON, Proj.Class, Projections, DIRECTIONAL.EFFECT) %>%
  spread(., Proj.Class, Projections)

dir_eff_func<- function(base, future, scale){
  if(scale == "Raw"){
    diff<- data.frame("x" = base$x, "y" = base$y, "Raw.Diff" = future$Projection - base$Projection)
    return(diff)
  }
  
  if(scale == "Percent"){
    base.adj<- ifelse(base$Projection == 0, 0.0001, base$Projection)
    diff<- data.frame("x" = base$x, "y" = base$y, "Perc.Diff" = (100*((future$Projection - base.adj)/base.adj)))
    return(diff)
  }
  
  if(scale == "Percent.New"){
    base.mean<- mean(base$Projection, na.rm = TRUE)
    fut.mean<- mean(future$Projection, na.rm = TRUE)
    diff<- data.frame("x" = 1, "y" = 1, "Perc.Diff" = (100*((fut.mean - base.mean)/base.mean)))
    return(diff)
  }
  
}

dat.dir<- dat.dir %>%
  mutate(., "Mean.Raw.Diff" = pmap(list(base = Baseline.sdm.b, future = Future_mean.sdm.b, scale = list("Raw")), possibly(dir_eff_func, NA)),
         "Pct05.Raw.Diff" = pmap(list(base = Baseline.sdm.b, future = Future_pct05.sdm.b, scale = list("Raw")), possibly(dir_eff_func, NA)),
         "Pct95.Raw.Diff" = pmap(list(base = Baseline.sdm.b, future = Future_pct95.sdm.b, scale = list("Raw")), possibly(dir_eff_func, NA)),
         "Mean.Perc.Diff" = pmap(list(base = Baseline.sdm.b, future = Future_mean.sdm.b, scale = list("Percent")), possibly(dir_eff_func, NA)),
         "Pct05.Perc.Diff" = pmap(list(base = Baseline.sdm.b, future = Future_pct05.sdm.b, scale = list("Percent")), possibly(dir_eff_func, NA)),
         "Pct95.Perc.Diff" = pmap(list(base = Baseline.sdm.b, future = Future_pct95.sdm.b, scale = list("Percent")), possibly(dir_eff_func, NA)),
         "Mean.PercNew.Diff" = pmap(list(base = Baseline.sdm.b, future = Future_mean.sdm.b, scale = list("Percent.New")), possibly(dir_eff_func, NA)),
         "Pct05.PercNew.Diff" = pmap(list(base = Baseline.sdm.b, future = Future_pct05.sdm.b, scale = list("Percent.New")), possibly(dir_eff_func, NA)),
         "Pct95.PercNew.Diff" = pmap(list(base = Baseline.sdm.b, future = Future_pct95.sdm.b, scale = list("Percent.New")), possibly(dir_eff_func, NA)))

dat.dir<- dat.dir %>%
  dplyr::select(., COMNAME, SEASON, DIRECTIONAL.EFFECT, Mean.Raw.Diff, Pct05.Raw.Diff, Pct95.Raw.Diff, Mean.Perc.Diff, Pct05.Perc.Diff, Pct95.Perc.Diff, Mean.PercNew.Diff, Pct05.PercNew.Diff, Pct95.PercNew.Diff) %>%
  gather(., Proj.Class, Projections, -COMNAME, -SEASON, -DIRECTIONAL.EFFECT)

# Get average
means_func<- function(df){
  temp<- data.frame(df)
  mean(temp[,3], na.rm = T)
}

dat.dir<- dat.dir %>%
  mutate(., "Mean.Diff" = as.numeric(map(Projections, possibly(means_func, NA))))

dat.dir$DIRECTIONAL.EFFECT<- factor(dat.dir$DIRECTIONAL.EFFECT, levels = c("Negative", "Neutral", "Positive"))

dat.dir.avg<- dat.dir %>%
  group_by(., Proj.Class, DIRECTIONAL.EFFECT) %>%
  summarise_at(., "Mean.Diff", c("mean", "sd"), na.rm = T) %>%
  drop_na(DIRECTIONAL.EFFECT)
dat.dir.avg$Proj.Class<- factor(dat.dir.avg$Proj.Class, levels = c("Pct05.Raw.Diff", "Mean.Raw.Diff", "Pct95.Raw.Diff", "Pct05.Perc.Diff", "Mean.Perc.Diff", "Pct95.Perc.Diff", "Pct05.PercNew.Diff", "Mean.PercNew.Diff", "Pct95.PercNew.Diff"))
dat.dir.avg$DIRECTIONAL.EFFECT<- factor(dat.dir.avg$DIRECTIONAL.EFFECT, levels = c("Negative", "Neutral", "Positive"))

plot.out.dir.raw<- ggplot(data = subset(dat.dir.avg, Proj.Class %in% c("Pct05.Raw.Diff", "Mean.Raw.Diff", "Pct95.Raw.Diff")), aes(x = DIRECTIONAL.EFFECT, y = mean)) +
  geom_point() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  ylab("SDM Mean Projected Change") +
  xlab("NEVA Directional Effect Classification") +
  facet_wrap(~Proj.Class, nrow = 1, scales = "free_y") +
  geom_smooth(data = subset(dat.dir.avg, Proj.Class %in% c("Pct05.Raw.Diff", "Mean.Raw.Diff", "Pct95.Raw.Diff")), aes(x = as.numeric(DIRECTIONAL.EFFECT), y = mean), method = "lm", se = FALSE) +
  theme_bw() 
ggsave(paste(out.dir, "DirectionaEffect_vs_RawChange.jpg", sep = ""), plot.out.dir.raw, width = 11, height = 8, units = "in")

plot.out.dir.perc<- ggplot(data = subset(dat.dir.avg, Proj.Class %in% c("Pct05.Perc.Diff", "Mean.Perc.Diff", "Pct95.Perc.Diff")), aes(x = DIRECTIONAL.EFFECT, y = mean)) +
  geom_point() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  ylab("SDM Mean Projected Percent Change") +
  xlab("NEVA Directional Effect Classification") +
  facet_wrap(~Proj.Class, nrow = 1, scales = "free_y") +
  geom_smooth(data = subset(dat.dir.avg, Proj.Class %in% c("Pct05.Perc.Diff", "Mean.Perc.Diff", "Pct95.Perc.Diff")), aes(x = as.numeric(DIRECTIONAL.EFFECT), y = mean), method = "lm", se = FALSE) +
  theme_bw() 
ggsave(paste(out.dir, "DirectionaEffect_vs_PercChange.jpg", sep = ""), plot.out.dir.perc, width = 11, height = 8, units = "in")

plot.out.dir.percnew<- ggplot(data = subset(dat.dir.avg, Proj.Class %in% c("Pct05.PercNew.Diff", "Mean.PercNew.Diff", "Pct95.PercNew.Diff")), aes(x = DIRECTIONAL.EFFECT, y = mean)) +
  geom_point() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  ylab("SDM Mean Projected Percent Change") +
  xlab("NEVA Directional Effect Classification") +
  facet_wrap(~Proj.Class, nrow = 1, scales = "free_y") +
  geom_smooth(data = subset(dat.dir.avg, Proj.Class %in% c("Pct05.PercNew.Diff", "Mean.PercNew.Diff", "Pct95.PercNew.Diff")), aes(x = as.numeric(DIRECTIONAL.EFFECT), y = mean), method = "lm", se = FALSE) +
  theme_bw() 
ggsave(paste(out.dir, "DirectionaEffect_vs_PercNewChange.jpg", sep = ""), plot.out.dir.percnew, width = 11, height = 8, units = "in")

## What about looking into species that follow vs. deviate from expectations?
# Directional effect
dir.cuts<- list(c(-5000, -25, 25, 5000))

for(i in seq_along(dir.cuts)){
  dir.cuts.use<- dir.cuts[[i]]
  
  dat.dir.comp<- dat.dir %>%
    filter(., Proj.Class == "Mean.PercNew.Diff") %>%
    group_by(., COMNAME, DIRECTIONAL.EFFECT, Proj.Class) %>%
    summarize_at(., "Mean.Diff", mean, na.rm = T)
  
  dat.dir.comp$Binned<- cut(dat.dir.comp$Mean.Diff, breaks = dir.cuts.use, labels = c("Negative", "Neutral", "Positive"))
  dat.dir.comp$Binned<- factor(dat.dir.comp$Binned, levels = c("Negative", "Neutral", "Positive"))
  dat.dir.comp$DIRECTIONAL.EFFECT<- factor(dat.dir.comp$DIRECTIONAL.EFFECT, levels = c("Negative", "Neutral", "Positive"))
  dat.dir.comp$Match<- paste(dat.dir.comp$DIRECTIONAL.EFFECT, dat.dir.comp$Binned, sep = "_")
  
  dat.dir.comp.table.hare<- dat.dir.comp %>%
    group_by(., DIRECTIONAL.EFFECT) %>%
    dplyr::summarize(., "Count.Hare" = n())
  
  dat.dir.comp.table.us<- dat.dir.comp %>%
    group_by(., Binned) %>%
    dplyr::summarize(., "Count.Us" = n())
  
  dat.dir.comp.count<- dat.dir.comp %>%
    drop_na(., DIRECTIONAL.EFFECT, Binned) %>%
    group_by(., Match) %>%
    dplyr::summarize(., "Count" = n()) %>%
    mutate(., "Percent" = round(100*(Count/sum(Count)), 2)) %>%
    separate(Match, c("DIRECTIONAL.EFFECT", "Binned"), sep = "_")
  #dat.dir.comp.count$Title<- paste(dat.dir.comp.count$Count, "\n", "(", dat.dir.comp.count$Percent, ")", sep = "")
  
  dat.dir.comp.tile<- data.frame("DIRECTIONAL.EFFECT" = c("Negative\n(n = 18)", "Neutral\n(n = 9)", "Positive\n(n = 6)"), "Binned" = c("Negative\n(n = 9)", "Neutral\n(n = 7)", "Positive\n(n = 17)"))
  
  dat.dir.comp.count$DIRECTIONAL.EFFECT<- factor(dat.dir.comp.count$DIRECTIONAL.EFFECT, levels = c("Negative", "Neutral", "Positive"), labels = c("Negative\n(n = 18)", "Neutral\n(n = 9)", "Positive\n(n = 6)"))
  dat.dir.comp.count$Binned<- factor(dat.dir.comp.count$Binned, levels = c("Negative", "Neutral", "Positive"), labels = c("Negative\n(n = 9)", "Neutral\n(n = 7)", "Positive\n(n = 17)"))
  
  dat.dir.comp.plot<- ggplot() +
    geom_text(data = dat.dir.comp.count, aes(x = DIRECTIONAL.EFFECT, y = Binned, label = NA), fontface = "bold") +
    geom_tile(data = dat.dir.comp.tile, aes(x = DIRECTIONAL.EFFECT, y = Binned), fill = NA, color = "black") +
    xlab("NEVA Directional Effect") +
    ylab("SDM Directional Effect") +
    theme(panel.background = element_blank())
  
  dat.dir.labels<- dat.dir.comp %>%
    drop_na(., DIRECTIONAL.EFFECT, Binned) %>%
    left_join(., func.groups)
  
  dat.dir.labels$Functional.Group<- factor(dat.dir.labels$Functional.Group, levels = c("Groundfish", "Pelagic", "Coastal", "Invertebrates", "Diadromous", "Elasmobranch"))
  
  dat.dir.labels$DIRECTIONAL.EFFECT<- factor(dat.dir.labels$DIRECTIONAL.EFFECT, levels = c("Negative", "Neutral", "Positive"), labels = c("Negative\n(n = 18)", "Neutral\n(n = 9)", "Positive\n(n = 6)"))
  dat.dir.labels$Binned<- factor(dat.dir.labels$Binned, levels = c("Negative", "Neutral", "Positive"), labels = c("Negative\n(n = 9)", "Neutral\n(n = 7)", "Positive\n(n = 17)"))
  
  set.seed(13331)
  dat.dir.comp.plot2<- dat.dir.comp.plot +
    geom_text_repel(data = dat.dir.labels, aes(x = DIRECTIONAL.EFFECT, y = Binned, label = toTitleCase(tolower(COMNAME)), color = Functional.Group), size = 4, segment.color = NA, show.legend = FALSE) +
    geom_point(data = dat.dir.labels, aes(x = DIRECTIONAL.EFFECT, y = Binned, color = Functional.Group), alpha = 0.00001) +
    scale_color_manual(name = "Functional Group", values = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')) +
    ggtitle(paste(dir.cuts.use[2], ":", dir.cuts.use[3], " used for neutral bin", sep = "")) +
    guides(color = guide_legend(override.aes = list(size = 3, color = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'), alpha = 1))) +
    theme(panel.background = element_blank())
  ggsave(paste(out.dir, dir.cuts.use[2], "to", dir.cuts.use[3], "DirectionaEffectMatching.jpg", sep = ""), dat.dir.comp.plot2, width = 15, height = 10, units = "in")
}

# Sensitivity...
dat.sens.comp<- dat.sens.keep %>%
  dplyr::filter(., Impact.Scenario == "Mean.SensitivityWeight") %>%
  group_by(., COMNAME, Impact.Scenario, SENSITIVITY, SENSITIVITY.CERTAINTY, VULNERABILITY, VULNERABILITY.CERTAINTY, DIRECTIONAL.EFFECT, DIRECTIONAL.EFFECT.CERTAINTY) %>%
  summarize_at(., "Data", mean, na.rm = T)
sens.quants<- quantile(dat.sens.comp$Data, na.rm = TRUE)
dat.sens.comp$Binned<- cut(dat.sens.comp$Data, breaks = sens.quants, labels = c("Low", "Moderate", "High", "Very high"), include.lowest = TRUE)
dat.sens.comp$Binned<- factor(dat.sens.comp$Binned, levels = c("Low", "Moderate", "High", "Very high"))
dat.sens.comp$SENSITIVITY<- factor(dat.sens.comp$SENSITIVITY, levels = c("Low", "Moderate", "High", "Very high"))
dat.sens.comp$Match<- paste(dat.sens.comp$SENSITIVITY, dat.sens.comp$Binned, sep = "_")

dat.sens.comp.table.hare<- dat.sens.comp %>%
  group_by(., SENSITIVITY) %>%
  dplyr::summarize(., "Count.Hare" = n())

dat.sens.comp.table.us<- dat.sens.comp %>%
  group_by(., Binned) %>%
  dplyr::summarize(., "Count.Us" = n())

dat.sens.comp.count<- dat.sens.comp %>%
  drop_na(., SENSITIVITY, Binned) %>%
  group_by(., Match) %>%
  dplyr::summarize(., "Count" = n()) %>%
  mutate(., "Percent" = round(100*(Count/sum(Count)), 2)) %>%
  separate(Match, c("SENSITIVITY", "Binned"), sep = "_")
dat.sens.comp.count$Title<- paste(dat.sens.comp.count$Count, "\n", "(", dat.sens.comp.count$Percent, ")", sep = "")
dat.sens.comp.count$SENSITIVITY<- factor(dat.sens.comp.count$SENSITIVITY, levels = c("Low", "Moderate", "High", "Very high"), labels = c("Low\n(n = 12)", "Moderate\n(n = 8)", "High\n(n = 13)", "Very high\n(n = 0)"))
dat.sens.comp.count$Binned<- factor(dat.sens.comp.count$Binned, levels = c("Low", "Moderate", "High", "Very high"), labels = c("Low\n(n = 9)", "Moderate\n(n = 8)", "High\n(n = 8)", "Very high\n(n = 8)"))

dat.sens.comp.tile<- data.frame("SENSITIVITY" = c("Low\n(n = 12)", "Moderate\n(n = 8)", "High\n(n = 13)", "Very high\n(n = 0)"), "Binned" = c("Low\n(n = 9)", "Moderate\n(n = 8)", "High\n(n = 8)", "Very high\n(n = 8)"))

dat.sens.comp.tile$SENSITIVITY<- factor(dat.sens.comp.tile$SENSITIVITY, levels = c("Low\n(n = 12)", "Moderate\n(n = 8)", "High\n(n = 13)", "Very high\n(n = 0)"), labels = c("Low\n(n = 12)", "Moderate\n(n = 8)", "High\n(n = 13)", "Very high\n(n = 0)"))
dat.sens.comp.tile$Binned<- factor(dat.sens.comp.tile$Binned, levels = c("Low\n(n = 9)", "Moderate\n(n = 8)", "High\n(n = 8)", "Very high\n(n = 8)"), labels = c("Low\n(n = 9)", "Moderate\n(n = 8)", "High\n(n = 8)", "Very high\n(n = 8)"))

dat.sens.comp.plot<- ggplot() +
  geom_tile(data = dat.sens.comp.tile, aes(x = SENSITIVITY, y = Binned), fill = NA, color = "black") +
  geom_text(data = dat.sens.comp.count, aes(x = SENSITIVITY, y = Binned, label = NA), fontface = "bold") +
  xlab("NEVA Sensitivity") +
  ylab("SDM Sensitivity") +
  theme(panel.background = element_blank())

dat.sens.labels<- dat.sens.comp %>%
  drop_na(., SENSITIVITY, Binned) %>%
  left_join(., func.groups)

dat.sens.labels$Functional.Group<- factor(dat.sens.labels$Functional.Group, levels = c("Groundfish", "Pelagic", "Coastal", "Invertebrates", "Diadromous", "Elasmobranch"))

dat.sens.labels$SENSITIVITY<- factor(dat.sens.labels$SENSITIVITY, levels = c("Low", "Moderate", "High", "Very high"), labels = c("Low\n(n = 12)", "Moderate\n(n = 8)", "High\n(n = 13)", "Very high\n(n = 0)"))
dat.sens.labels$Binned<- factor(dat.sens.labels$Binned, levels = c("Low", "Moderate", "High", "Very high"), labels = c("Low\n(n = 9)", "Moderate\n(n = 8)", "High\n(n = 8)", "Very high\n(n = 8)"))

dat.sens.comp.plot2<- dat.sens.comp.plot +
  geom_text_repel(data = dat.sens.labels, aes(x = SENSITIVITY, y = Binned, label = toTitleCase(tolower(COMNAME)), color = Functional.Group), size = 4, segment.color = NA) +
  geom_point(data = dat.sens.labels, aes(x = SENSITIVITY, y = Binned, color = Functional.Group), alpha = 0.00001) +
  scale_color_manual(name = "Functional Group", values = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')) +
  guides(color = guide_legend(override.aes = list(size = 3, color = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'), alpha = 1))) +
  theme(panel.background = element_blank())
ggsave(paste(out.dir, "SensitivityMatching.jpg", sep = ""), dat.sens.comp.plot2, width = 17, height = 10, units = "in")

### Random Forest model to test match vs. mismatching
## Random forest for sensitivity
dat.sens.mv<- dat.sens.comp %>%
  ungroup() %>%
  dplyr::select(., COMNAME, Match)

dat.sens.mv$Match.Factor<- ifelse(dat.sens.mv$Match == "Low_Low" | dat.sens.mv$Match == "Moderate_Moderate"| dat.sens.mv$Match == "High_High" | dat.sens.mv$Match == "Very high_Very high", "Match",
                                  ifelse(dat.sens.mv$Match == "Moderate_Low" | dat.sens.mv$Match == "High_Moderate" | dat.sens.mv$Match == "Very high_High", "Miss One NEVA", 
                                         ifelse(dat.sens.mv$Match == "Low_Moderate" |  dat.sens.mv$Match == "Moderate_High" |  dat.sens.mv$Match == "High_Very high", "Miss One SDM",
                                                ifelse(dat.sens.mv$Match == "High_Low" | dat.sens.mv$Match == "Very high_Moderate", "Miss Two NEVA",
                                                       ifelse(dat.sens.mv$Match == "Low_High" | dat.sens.mv$Match == "Moderate_Very high", "Miss Two SDM",
                                                              ifelse(dat.sens.mv$Match == "Very high_Low", "Miss Three NEVA",
                                                                     ifelse(dat.sens.mv$Match == "Low_Very high", "Miss Three SDM", NA)))))))

dat.sens.mv$Match.Factor<- factor(dat.sens.mv$Match.Factor, levels = c("Match", "Miss One NEVA", "Miss One SDM", "Miss Two NEVA", "Miss Two SDM", "Miss Three NEVA", "Miss Three SDM"))
dat.sens.mv$Match.Group<- ifelse(grepl("NEVA", dat.sens.mv$Match.Factor), "NEVA",
                                 ifelse(grepl("SDM", dat.sens.mv$Match.Factor), "SDM", "None"))
dat.sens.mv$Match.Group<- factor(dat.sens.mv$Match.Group, levels = c("None", "SDM", "NEVA"))
dat.sens.mv$Match.LevelOnly<- str_replace_all(dat.sens.mv$Match.Factor, c(" NEVA" = "", " SDM" = ""))
dat.sens.mv$Match.LevelOnly<- factor(dat.sens.mv$Match.LevelOnly, levels = c("Match", "Miss One", "Miss Two", "Miss Three"))

# Add in attributes...
dir.eff.dat<- read.csv(here::here("Data", "JHareDirectionalEffect.csv"))
qual.dat<- read.csv(here::here("Data", "JHareQualitativeDataResults.csv"))

# Need to get one directional effect for each species and one rank for each species attribute for sensitivity and exposure...
# Wide to long format
qual.dat.l<- gather(qual.dat, "Score", "Votes", Low, Moderate, High, Very.High)
qual.dat.l<- arrange(qual.dat.l, Species, Functional.Group, Attribute, Score)

qual.dat.l$Score<- factor(qual.dat.l$Score, levels = c("Low", "Moderate", "High", "Very.High"))
qual.dat.l$Score.Numeric<- as.numeric(qual.dat.l$Score)

# Calculate weighted mean by species-attribute-attribute.category, then characterize these into low, moderate, high and very high. This will allow us to count them and apply the logic rule used by Hare et al.
qual.dat.wt.mean<- qual.dat.l %>%
  group_by(., Species, Attribute, Attribute.Category) %>%
  dplyr::summarise("Weighted.Mean" = weighted.mean(Score.Numeric, Votes)) %>%
  dplyr::mutate(., "Logic.Rule" = cut(Weighted.Mean, breaks = c(0, 2.5, 3.0, 3.5, 100), labels = c("Low", "Moderate", "High", "Very.High")))

# Now need to get to overall vulnerability...
# Get the counts (number of votes in each species - attribute category - vulneraibity) to apply logic rule
dat.scores<- qual.dat.wt.mean %>%
  group_by(., Species, Attribute.Category, Logic.Rule) %>%
  dplyr::summarise("Number" = n()) %>%
  data.frame(.)

# Apply logic rule to translate number of votes to a low-moderatre-high-very high rank for each species - attribute category
spp.ranks.split<- dat.scores %>%
  dplyr::mutate(Rank.Temp = as.numeric(ifelse(Logic.Rule == "Very.High" & Number >= 3, 4, 
                                              ifelse(Logic.Rule == "High" | Logic.Rule == "Very.High" & Number >= 2, 3, 
                                                     ifelse(Logic.Rule == "High" | Logic.Rule == "Very.High" | Logic.Rule == "Moderate" & Number >= 2, 2, 1))))) %>%
  dplyr::group_by(., Species, Attribute.Category) %>%
  dplyr::summarise(Rank = max(Rank.Temp)) 

# Calculate one overall vulnerability rank for each species, which is exposure factor rank * sensitivity attribute rank (ranges 1-16)
spp.ranks.overall<- spp.ranks.split %>%
  spread(., Attribute.Category, Rank) %>% 
  dplyr::group_by(., Species) %>%
  dplyr::summarise(Overall.Rank = Exposure.Factor*Sensitivity.Attribute) %>% 
  dplyr::mutate(., "Overall.Rank.Code" = cut(Overall.Rank, breaks = c(0, 3, 6, 9, 16), labels = c("Low", "Moderate", "High", "Very.High")))

# Directional effect
dir.eff.dat.l<- gather(dir.eff.dat, "Score", "Votes", Negative, Neutral, Positive)
dir.eff.dat.l<- arrange(dir.eff.dat.l, Species, Functional.Group, Score)

# Convert directional effect to number
dir.eff.dat.l$Score.Numeric<- ifelse(dir.eff.dat.l$Score == "Negative", -1,
                                     ifelse(dir.eff.dat.l$Score == "Neutral", 0, 1))

# Calculate weighted mean by species, then characterize these into negative. neutral, positive. 
dir.eff.dat.wt.mean<- dir.eff.dat.l %>%
  group_by(., Species) %>%
  dplyr::summarise("Weighted.Mean" = weighted.mean(Score.Numeric, Votes)) %>%
  dplyr::mutate(., "Directional.Effect" = cut(Weighted.Mean, breaks = c(-100, -0.333, 0.333, 100), labels = c("Negative", "Neutral", "Positive")))

# Okay...now what? Let's add in the response columns, which is the SDM and NEVA rank to each of the qual.dat.wt.mean and dir.eff.dat.wt.mean dataframes
qual.dat.wt.mean$Attribute.Full<- paste(qual.dat.wt.mean$Attribute.Category, "_", qual.dat.wt.mean$Attribute, sep = "")
qual.dat.mv<- qual.dat.wt.mean %>%
  ungroup() %>%
  dplyr::select(., -Weighted.Mean, -Attribute, -Attribute.Category) %>%
  spread(., Attribute.Full, Logic.Rule) 
qual.dat.mv$COMNAME<- toupper(as.character(qual.dat.mv$Species))

dat.sens.mod<- dat.sens.mv %>%
  left_join(., qual.dat.mv, by = "COMNAME") %>%
  drop_na(Match.Factor, Match.LevelOnly)

dat.sens.mod<- data.frame(dat.sens.mod[,c(5, 19:30)])

# Can't have any 0 obs for response classes...
dat.mod.full<- data.frame(dat.sens.mod[,])

rf.f<- randomForest(Match.LevelOnly ~ ., data = dat.mod.full, importance = TRUE)
rf.f
varImpPlot(rf.f)

# Not great, reduce down to just matching vs. mismatching
dat.sens.mod$Match.LevelRed<- ifelse(grepl("Miss", dat.sens.mod$Match.LevelOnly), "Miss", "Match")
dat.sens.mod$Match.LevelRed<- factor(dat.sens.mod$Match.LevelRed, levels = c("Match", "Miss"))

dat.mod.LevelRed<- dat.sens.mod[,c(14, 2:13)]
rf.f.red<- randomForest(Match.LevelRed ~ ., data = dat.mod.LevelRed, importance = TRUE)
rf.f.red
varImpPlot(rf.f.red)

# Partial variable importance
pop.gr<- partialPlot(rf.f.red, dat.mod.LevelRed, Sensitivity.Attribute_Population.Growth.Rate, "Miss")
adult.mob<- partialPlot(rf.f.red, dat.mod.LevelRed, Sensitivity.Attribute_Adult.Mobility, "Miss")

part.df<- data.frame("Category" = c(rep("Population.Growth.Rate", 4), rep("Adult.Mobility", 4)), "Rank" = c(pop.gr[[1]], adult.mob[[1]]), "Probability" = c(pop.gr[[2]], adult.mob[[2]]))
part.df<- part.df %>%
  complete(., Category, Rank)
part.df$Rank<- factor(part.df$Rank, levels = c("Low", "Moderate", "High", "Very.High"))
part.df$Probability<- ifelse(part.df$Probability > 1, 1, part.df$Probability)

part.plot<- ggplot() +
  geom_bar(data = part.df, aes(x = Rank, y = Probability, fill = Category), stat = "identity", position = "dodge") +
  scale_fill_manual(name = "Trait/Factor", values = c('#66a61e','#e6ab02')) +
  ylab("Marginal probability of being classified as a 'Miss'") +
  xlab("NEVA Sensitivity")
ggsave(paste(out.dir, "RFPartImpPlot.jpg", sep = ""), part.plot, width = 11, height = 8)


## Random forest for directional effect...
dat.dir.mv<- dat.dir.comp %>%
  ungroup() %>%
  dplyr::select(., COMNAME, Match)

dat.dir.mv$Match.Factor<- ifelse(dat.dir.mv$Match == "Negative_Negative" | dat.dir.mv$Match == "Neutral_Neutral"| dat.dir.mv$Match == "Positive_Positivie", "Match",
                                 ifelse(dat.dir.mv$Match == "Negative_Neutral" | dat.dir.mv$Match == "Positive_Neutral", "Miss One NEVA", 
                                        ifelse(dat.dir.mv$Match == "Neutral_Negative" |  dat.dir.mv$Match == "Neutral_Positive", "Miss One SDM",
                                               ifelse(dat.dir.mv$Match == "Positive_Negative", "Miss Two NEVA",
                                                      ifelse(dat.dir.mv$Match == "Negative_Positive", "Miss Two SDM", NA)))))

dat.dir.mv$Match.Factor<- factor(dat.dir.mv$Match.Factor, levels = c("Match", "Miss One NEVA", "Miss One SDM", "Miss Two NEVA", "Miss Two SDM"))
dat.dir.mv$Match.Group<- ifelse(grepl("NEVA", dat.dir.mv$Match.Factor), "NEVA",
                                ifelse(grepl("SDM", dat.dir.mv$Match.Factor), "SDM", "None"))
dat.dir.mv$Match.Group<- factor(dat.dir.mv$Match.Group, levels = c("None", "SDM", "NEVA"))
dat.dir.mv$Match.LevelOnly<- str_replace_all(dat.dir.mv$Match.Factor, c(" NEVA" = "", " SDM" = ""))
dat.dir.mv$Match.LevelOnly<- factor(dat.dir.mv$Match.LevelOnly, levels = c("Match", "Miss One", "Miss Two"))

dat.dir.mod<- dat.dir.mv %>%
  left_join(., qual.dat.mv, by = "COMNAME") %>%
  drop_na(Match.Factor, Match.LevelOnly)

dat.dir.mod<- data.frame(dat.dir.mod[,c(5, 19:30)])

# Can't have any 0 obs for response classes...
summary(dat.dir.mod)
dat.dir.mod<- droplevels(dat.dir.mod)

dat.mod.full<- data.frame(dat.dir.mod[,c(1, 2:13)])

rf.f<- randomForest(Match.LevelOnly ~ ., data = dat.mod.full, importance = TRUE)
rf.f
varImpPlot(rf.f)

# Not great, reduce down to just matching vs. mismatching
dat.dir.mod$Match.LevelRed<- ifelse(grepl("Miss", dat.dir.mod$Match.LevelOnly), "Miss", "Match") 
dat.dir.mod$Match.LevelRed<- factor(dat.dir.mod$Match.LevelRed, levels = c("Match", "Miss"))

dat.mod.LevelRed<- data.frame(dat.dir.mod[,c(14, 2:13)])
rf.f.red<- randomForest(Match.LevelRed ~ ., data = dat.mod.LevelRed, importance = TRUE)
rf.f.red
varImpPlot(rf.f.red)

# Partial variable importance
other.stress<- partialPlot(rf.f.red, dat.mod.LevelRed, Sensitivity.Attribute_Other.Stressors, "Miss")
stock.size<- partialPlot(rf.f.red, dat.mod.LevelRed, Sensitivity.Attribute_Stock.Size.Status, "Miss")
pop.gr<- partialPlot(rf.f.red, dat.mod.LevelRed, Sensitivity.Attribute_Population.Growth.Rate, "Miss")
complex<- partialPlot(rf.f.red, dat.mod.LevelRed, Sensitivity.Attribute_Complexity.in.Reproductive.Strategy, "Miss")

part.df<- data.frame("Category" = c(rep("Other.Stressors", 4), rep("Stock.Size.Status", 4), rep("Population.Growth.Rate", 4), rep("Complexity.In.Reproductive.Strategy", 4)), "Rank" = rep(c("Low", "Moderate", "High", "Very.High"), 4), "Probability" = c(c(other.stress[[2]], 0), stock.size[[2]], pop.gr[[2]], c(complex[[2]], 0)))
part.df<- part.df %>%
  complete(., Category, Rank)
part.df$Rank<- factor(part.df$Rank, levels = c("Low", "Moderate", "High", "Very.High"))
part.df$Probability<- ifelse(part.df$Probability > 1, 1, part.df$Probability)

part.plot<- ggplot() +
  geom_bar(data = part.df, aes(x = Rank, y = Probability, fill = Category), stat = "identity", position = "dodge") +
  scale_fill_manual(name = "Trait/Factor", values = c('#1b9e77','#d95f02','#7570b3','#e7298a')) +
  ylab("Marginal probability of being classified as a 'Miss'") +
  xlab("NEVA Sensitivity")
ggsave(paste(out.dir, "RFPartImpPlotNew.jpg", sep = ""), part.plot, width = 11, height = 8)

## Random forest searching.....
# Establish a list of possible values for mtry, nodesize and sampsize
mtry<- seq(2, ncol(dat.mod.f) * 0.8, 2)
nodesize<- seq(2, 8, 2)
sampsize<- nrow(dat.mod.f) * c(0.7, 0.8)

# Create a data frame containing all combinations 
hyper.grid<- expand.grid(mtry = mtry, nodesize = nodesize, sampsize = sampsize)

# Create an empty vector to store OOB error values
oob_err<- c()

# Write a loop over the rows of hyper_grid to train the grid of models
for(i in 1:nrow(hyper.grid)) {
  
  # Train a Random Forest model
  model<- randomForest(formula = Match.Bin ~ ., 
                       data = dat.mod.f,
                       mtry = hyper.grid$mtry[i],
                       nodesize = hyper.grid$nodesize[i],
                       sampsize = hyper.grid$sampsize[i])
  
  # Store OOB error for the model                      
  oob_err[i]<- model$err.rate[nrow(model$err.rate), "OOB"]
}

# Identify optimal set of hyperparmeters based on OOB error
opt_i <- which.min(oob_err)
print(hyper.grid[opt_i,])

# Nothing new, okay to keep rf.f.bin.
# Results — Combo vs. SDM predictive ability to baseline ------------------
# Calculate difference in key statistics between the two models
# Maybe depending on certainty of SDM and on Jon's Assessment?
# Jon Certainty Data
jon.df<- read_csv(here::here("Data", "JHareQualitativeResults.csv"))

# Join to model results
mod.res.all<- mod.res %>%
  left_join(., jon.df, by = "COMNAME") 

# Overall SDM certainty?
auc.splits<- seq(0.7, 1.0, length.out = 5)
mod.res.all$SDM.Cert<- cut(mod.res.all$AUC.SDM, auc.splits, labels = c("Low", "Moderate", "High", "Very high"), include.lowest = TRUE)

# Overall NEVA certainty
cert.combos<- expand.grid(c("Low", "Moderate", "High", "Very high"), c("Low", "Moderate", "High", "Very high"))
cert.scores<- expand.grid(c(1, 2, 3, 4), c(1, 2, 3, 4))
ovrall.cert.df<- data.frame("Both.Certs" = paste(cert.combos$Var1, cert.combos$Var2, sep = "_"), "Cert.Score" = cert.scores$Var1 + cert.scores$Var2)
cert.splits<- seq(2, 8, length.out = 5)
ovrall.cert.df$Cert.Score.Rank<- cut(ovrall.cert.df$Cert.Score, cert.splits, labels = c("Low", "Moderate", "High", "Very high"), include.lowest = TRUE)

# Join?
mod.res.all$Both.Certs<- paste(mod.res.all$DIRECTIONAL.EFFECT.CERTAINTY, mod.res.all$SENSITIVITY.CERTAINTY, sep = "_")
mod.res.all<- mod.res.all %>%
  left_join(., ovrall.cert.df)

mod.diffs<- mod.res.all %>%
  mutate(., 
            "CorrCoeff.Diff" = CorrCoeff.NEVA.B - CorrCoeff.SDM.B,
            "Bias.Diff" = Bias.NEVA.B - Bias.SDM.B,
            "RMSE.Diff" = RMSE.NEVA.B - RMSE.SDM.B) %>%
  group_by(SDM.Cert, Cert.Score.Rank) %>%
  summarize(., 
            "MeanCorrCoeff.Diff" = round(mean(CorrCoeff.Diff), 3), 
            "MeanBias.Diff" = round(mean(Bias.Diff), 3),
            "MeanRMSE.Diff" = round(mean(RMSE.Diff), 3)) %>%
  drop_na()
write_csv(mod.diffs, paste(out.dir, "NEVAminusSDMonlyModelPredictiveAbility", sep = ""))
# Results — Combo vs. SDM prediction intervals ----------------------------
## The overall idea here is to see if NEVA models improve the predictions (accuracy or uncertainty)
rmse_func<- function(predicted, test.data) {
  if(all(test.data$BIOMASS == 0)){
    return(NA)
  } else {
    return(nrmse(sim = as.numeric(as.data.frame(predicted)[,1]), obs = test.data$BIOMASS))
  }
}

# Going to need access to the likelihood, prior and posterior information. This should be all in the mcmc results file for each species...
res.files<- list.files(out.dir, "mcmc_")

# Results empty storage file
neva.vs.sdm<- data.frame("COMNAME" = rep(NA, length(res.files)), "SEASON" = rep(NA, length(res.files)), "SDM.SD.mu" = rep(NA, length(res.files)), "NEVA.SD.mu" = rep(NA, length(res.files)), "NEVA.RMSE.mu" = rep(NA, length(res.files)), "NEVA.RMSE.sd" = rep(NA, length(res.files)), "SDM.RMSE.mu" = rep(NA, length(res.files)), "SDM.RMSE.sd" = rep(NA, length(res.files)))

# Loop over each file -- species season model fit information
for(i in seq_along(res.files)){
  
  print(res.files[i])
  
  # SDM Model fit -- presence and biomass
  gam.p.temp<- readRDS(paste(out.dir, gsub("mcmc_", "gamfitpres", res.files[i]), sep = ""))
  gam.b.temp<- readRDS(paste(out.dir, gsub("mcmc_", "gamfitbio", res.files[i]), sep = ""))
  ilink<- family(gam.b.temp)$linkinv
  gam.coef<- names(coef(gam.p.temp))
  
  # Ranks of candidate draws...
  # Get ranks for SDM only and for NEVA
  likes.temp<- data.frame(read_rds(paste(out.dir, res.files[i], sep = ""))[[1]])
  names(likes.temp)<- c("Likelihood", "Prior", "Posterior")
  likes.temp$Iteration<- rep(seq(from = 1, to = nrow(likes.temp), by = 1))
  likes.temp$SDM.Rank<- order(likes.temp$Prior)
  likes.temp$NEVA.Rank<- order(likes.temp$Posterior)
  
  # Need to make predictions from each of these iterations
  # Getting model coefficients for each of the candidate models
  mods.temp<- data.frame(read_rds(paste(out.dir, res.files[i], sep = ""))[[2]])
  colnames(mods.temp)<- gam.coef
  
  # Determining species-season which will then be used to get the testing data to validate predictions
  spp.match<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][1])
  season.match<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][2])
  col.check<- paste(spp.match, season.match, sep = ".")
  
  if(grepl("sturgeon", res.files[i])){
    neva.vs.sdm$COMNAME[i]<- spp.match
    neva.vs.sdm$SEASON[i]<- season.match
    neva.vs.sdm$SDM.SD.mu[i]<- NA
    neva.vs.sdm$NEVA.SD.mu[i]<- NA
    neva.vs.sdm$NEVA.RMSE.mu[i]<- NA
    neva.vs.sdm$NEVA.RMSE.sd[i]<- NA
    neva.vs.sdm$SDM.RMSE.mu[i]<- NA
    neva.vs.sdm$SDM.RMSE.sd[i]<- NA
    next()
  }
  
  # Collecting testing data
  test.data<- dat.full %>%
    filter(., COMNAME == spp.match & SEASON == season.match) %>%
    dplyr::select(TEST.DATA) %>%
    unnest() %>%
    data.frame()
  temp<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS")))
  test.data<- data.frame(na.omit(temp))
  
  # Predicting presence/absence, does NOT change with different candidate models, which only influence biomass curves
  pred.p<- predict.gam(gam.p.temp, newdata = test.data, type = "response")
  
  # Predicting biomass component for each of the 1000 candidate models
  # Prediction storage
  preds.out<- data.frame("Pt" = seq(from = 1, nrow(test.data)), "Pred" = rep(NA, nrow(test.data)))
  
  for(j in 1:nrow(mods.temp)){
    fit.mat<- matrix(as.numeric(mods.temp[j,]), nrow = 1, ncol = length(mods.temp[j,]), byrow = T, dimnames = list(NULL, gam.coef))
    
    # Make predictions with these values
    lpmat.pred<- predict.gam(gam.b.temp, newdata = test.data, type = "lpmatrix")
    preds.out[,j+1]<- as.numeric(exp(ilink(as.numeric(lpmat.pred %*% t(fit.mat))))*pred.p)
    names(preds.out)[j+1]<- paste("Iteration.", j, sep = "")
  }
  
  # Get the SDM predictions...
  sdm.best<- likes.temp[likes.temp$SDM.Rank<= 100, ]
  sdm.preds.best<- preds.out[,sdm.best$Iteration+1] %>%
    gather(., "Iteration", "Prediction") %>%
    separate(., Iteration, c("Remove", "Iteration")) %>%
    dplyr::select(., -Remove) %>%
    group_by(Iteration) %>%
    nest(., .key = "Pred")
  
  # Get the NEVA predictions
  neva.best<- likes.temp[likes.temp$NEVA.Rank<= 100, ]
  neva.preds.best<- preds.out[,neva.best$Iteration+1] %>%
    gather(., "Iteration", "Prediction") %>%
    separate(., Iteration, c("Remove", "Iteration")) %>%
    dplyr::select(., -Remove) %>%
    group_by(Iteration) %>%
    nest(., .key = "Pred") 
  
  # Prediction uncertainty -- for the top 85% of each method calculate the average per cell variability
  sdm.preds.best.sd<- preds.out[,sdm.best$Iteration+1] %>%
    mutate(., "Cell" = seq(from = 1, to = nrow(.))) %>%
    gather(., "Iteration", "Prediction", -Cell) %>%
    group_by(Cell) %>%
    dplyr::summarize(., "Cell.SD" = sd(Prediction, na.rm = T)) 
  sdm.preds.best.sd<- mean(sdm.preds.best.sd$Cell.SD)

  neva.preds.best.sd<- preds.out[,neva.best$Iteration+1] %>%
    mutate(., "Cell" = seq(from = 1, to = nrow(.))) %>%
    gather(., "Iteration", "Prediction", -Cell) %>%
    group_by(Cell) %>%
    dplyr::summarize(., "Cell.SD" = sd(Prediction, na.rm = T)) 
  neva.preds.best.sd<- mean(neva.preds.best.sd$Cell.SD)
  
  # For the top 85% of each method, validate each to get a RMSE and AUC value
  test.data<- dat.full %>%
    filter(., COMNAME == spp.match & SEASON == season.match) %>%
    dplyr::select(TEST.DATA) %>%
    unnest() %>%
    data.frame()
  temp<- dplyr::select(test.data, one_of(c("BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale")))
  test.data<- data.frame(na.omit(temp))
  
  sdm.summary<- sdm.preds.best %>%
    mutate(., "RMSE" = map2(Pred, list(test.data), possibly(rmse_func, NA))) %>%
    dplyr::select(., Iteration, RMSE) %>%
    unnest() %>%
    summarize_at(., vars("RMSE"), funs(mean, sd), na.rm = T)
  
  neva.summary<- neva.preds.best %>%
    mutate(., "RMSE" = map2(Pred, list(test.data), possibly(rmse_func, NA))) %>%
    dplyr::select(., Iteration, RMSE) %>%
    unnest() %>%
    summarize_at(., vars("RMSE"), funs(mean, sd), na.rm = T)
  
  neva.vs.sdm$COMNAME[i]<- spp.match
  neva.vs.sdm$SEASON[i]<- season.match
  neva.vs.sdm$SDM.SD.mu[i]<- sdm.preds.best.sd
  neva.vs.sdm$NEVA.SD.mu[i]<- neva.preds.best.sd
  neva.vs.sdm$NEVA.RMSE.mu[i]<- neva.summary$mean
  neva.vs.sdm$NEVA.RMSE.sd[i]<- neva.summary$sd
  neva.vs.sdm$SDM.RMSE.mu[i]<- sdm.summary$mean
  neva.vs.sdm$SDM.RMSE.sd[i]<- sdm.summary$sd
}

# Improvement by reducing uncertainty?
neva.vs.sdm.sd<- neva.vs.sdm$SDM.SD.mu - neva.vs.sdm$NEVA.SD.mu

# Improvement by increasing prediction accuracy?
rmse.cert.dat.mu<- neva.vs.sdm %>%
  dplyr::select(., COMNAME, SEASON, NEVA.RMSE.mu, SDM.RMSE.mu) %>%
  mutate(., "RMSE.mu.diff" = round(NEVA.RMSE.mu - SDM.RMSE.mu), 3) %>%
  dplyr::select(., COMNAME, SEASON, RMSE.mu.diff) 
rmse.cert.dat.sd<- neva.vs.sdm %>%
  dplyr::select(., COMNAME, SEASON, NEVA.RMSE.sd, SDM.RMSE.sd) %>%
  mutate(., "RMSE.sd.diff" = round(NEVA.RMSE.sd - SDM.RMSE.sd), 3) %>%
  dplyr::select(., COMNAME, SEASON, RMSE.sd.diff) 
rmse.cert.dat<- data.frame(rmse.cert.dat.mu, "RMSE.sd.diff" = rmse.cert.dat.sd$RMSE.sd.diff)

# Join with overall uncertainty...
certainty.levels<- data.frame("SENSITIVITY.CERTAINTY" = c("Low", "Moderate", "High", "Very high"), "DIRECTIONAL.EFFECT.CERTAINTY" = c("Low", "Moderate", "High", "Very high"), "Score" = c(1, 2, 3, 4))

# Jon's data
jon.df<- read_csv(here::here("Data", "JHareQualitativeResults.csv"))

rmse.cert.dat<- rmse.cert.dat %>%
  left_join(., jon.df, by = "COMNAME")

rmse.cert.dat<- rmse.cert.dat %>%
  mutate(., "Certainty" = paste(SENSITIVITY.CERTAINTY, "_", DIRECTIONAL.EFFECT.CERTAINTY, sep = ""))

rmse.cert.dat<- rmse.cert.dat %>%
  left_join(., dplyr::select(certainty.levels, c(SENSITIVITY.CERTAINTY, Score)), by = "SENSITIVITY.CERTAINTY")
names(rmse.cert.dat)[12]<- "SENSITIVITY.Score"

rmse.cert.dat<- rmse.cert.dat %>%
  left_join(., dplyr::select(certainty.levels, c(DIRECTIONAL.EFFECT.CERTAINTY, Score)), by = "DIRECTIONAL.EFFECT.CERTAINTY")
names(rmse.cert.dat)[13]<- "DirectionalEffect.Score"

rmse.cert.dat$Certainty.Score<- rmse.cert.dat$SENSITIVITY.Score + rmse.cert.dat$DirectionalEffect.Score
rmse.cert.dat$Certainty.Score<- factor(rmse.cert.dat$Certainty.Score, levels = c(2, 3, 4, 5, 6, 7, 8), labels = c("Low", "Low_Moderate", "Moderate", "Moderate_High", "High", "High_Very high", "Very high"))
rmse.cert.dat$RMSE.mu.diff.code<- ifelse(rmse.cert.dat$RMSE.mu.diff > 0, "Better", 
                                         ifelse(rmse.cert.dat$RMSE.mu.diff == 0, "No change", "Worse"))
rmse.cert.dat$RMSE.mu.diff.code<- factor(rmse.cert.dat$RMSE.mu.diff.code, levels = c("Better", "No change", "Worse"))
rmse.cert.dat$RMSE.sd.diff.code<- ifelse(rmse.cert.dat$RMSE.sd.diff > 0, "Better", 
                                         ifelse(rmse.cert.dat$RMSE.sd.diff == 0, "No change", "Worse"))
rmse.cert.dat$RMSE.sd.diff.code<- factor(rmse.cert.dat$RMSE.sd.diff.code, levels = c("Better", "No change", "Worse"))

rmse.cert.dat$SENSITIVITY.CERTAINTY<- factor(rmse.cert.dat$SENSITIVITY.CERTAINTY, levels = c("Low", "Moderate", "High", "Very high"))
rmse.cert.dat$DIRECTIONAL.EFFECT.CERTAINTY<- factor(rmse.cert.dat$DIRECTIONAL.EFFECT.CERTAINTY, levels = c("Low", "Moderate", "High", "Very high"))

rmse.cert.plot<- ggplot(na.omit(rmse.cert.dat), aes(x = SENSITIVITY.CERTAINTY, y = DIRECTIONAL.EFFECT.CERTAINTY, color = RMSE.mu.diff.code, label = str_to_title(COMNAME))) + 
  scale_color_manual(name = "Effect of including\nNEVA on RMSE", values = c("#4daf4a", "gray", "#e41a1c")) + 
  geom_text_repel(segment.color = NA, show.legend = F) +
  xlab("Sensitivity certainty") +
  ylab("Directional effect certainty") +
  facet_wrap(~SEASON, scales = "free", nrow = 2)
ggsave(paste(out.dir, "rmse.plot.jpg", sep = ""), height = 8, width = 11, units = "in")


# Results — Combo vs SDM parameter intervals ------------------------------
# Going to need access to the likelihood, prior and posterior information. This should be all in the mcmc results file for each species...
res.files<- list.files(out.dir, "mcmc_")

# Results
for(i in seq_along(res.files)){
  
  spp<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][1])
  season<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][2])
  
  # Model fit -- presence and biomass
  gam.b.temp<- readRDS(paste(out.dir, gsub("mcmc_", "gamfitbio", res.files[i]), sep = ""))
  ilink<- family(gam.b.temp)$linkinv
  gam.coef<- names(coef(gam.b.temp))
  
  # Ranks of candidate draws...
  # Get ranks for SDM only and for NEVA
  likes.temp<- data.frame(read_rds(paste(out.dir, res.files[i], sep = ""))[[1]])
  names(likes.temp)<- c("Likelihood", "Prior", "Posterior")
  likes.temp$Iteration<- rep(seq(from = 1, to = nrow(likes.temp), by = 1))
  likes.temp$SDM.Rank<- order(likes.temp$Prior)
  likes.temp$NEVA.Rank<- order(likes.temp$Posterior)
  
  # Need to make predictions from each of these iterations
  mods.temp<- data.frame(read_rds(paste(out.dir, res.files[i], sep = ""))[[2]])
  colnames(mods.temp)<- gam.coef
  
  # Get the SDM fits...
  sdm.best<- likes.temp[likes.temp$SDM.Rank<= 900, ]
  sdm.best<- sdm.best[order(sdm.best$SDM.Rank), ]
  sdm.best.params<- mods.temp[sdm.best$Iteration, ] %>%
    summarize_all(., funs(mean, sd), na.rm = T)
  
  # Get the NEVA predictions
  neva.best<- likes.temp[likes.temp$NEVA.Rank<= 900, ]
  neva.best<- neva.best[order(neva.best$NEVA.Rank), ]
  neva.best.params<- mods.temp[neva.best$Iteration, ] %>%
    summarize_all(., funs(mean, sd), na.rm = T)
  
  if(i == 1){
    sdm.result<- data.frame("COMNAME" = spp, "Season" = season, sdm.best.params)
    neva.result<- data.frame("COMNAME" = spp, "Season" = season, neva.best.params)
  } else {
    sdm.temp<- data.frame("COMNAME" = spp, "Season" = season, sdm.best.params)
    sdm.result<- bind_rows(sdm.result, sdm.temp)
    neva.temp<- data.frame("COMNAME" = spp, "Season" = season, neva.best.params)
    neva.result<- bind_rows(neva.result, neva.temp)
  }
}

# Make some plots...
sdm.plot<- sdm.result[,c(1:21)] %>%
  gather(., "Variable", "Mean", -COMNAME, -Season)
sdm.sd<- sdm.result[,c(1, 2, 22:40)] %>%
  gather(., "Variable", "SD", -COMNAME, -Season)
sdm.plot<- bind_cols(sdm.plot, sdm.sd[,-c(1:2)])
sdm.plot$Approach<- rep("SDM", nrow(sdm.plot))

neva.plot<- neva.result[,c(1:21)] %>%
  gather(., "Variable", "Mean", -COMNAME, -Season)
neva.sd<- neva.result[,c(1, 2, 22:40)] %>%
  gather(., "Variable", "SD", -COMNAME, -Season)
neva.plot<- bind_cols(neva.plot, neva.sd[,-c(1:2)])
neva.plot$Approach<- rep("SDM + NEVA", nrow(neva.plot))

params.plot<- bind_rows(sdm.plot, neva.plot)

sdm.plot.fall<- ggplot(na.omit(subset(params.plot, Season == "FALL" & Variable == "X.Intercept._mean")), aes(x = COMNAME, y = Mean, color = Approach)) + 
  scale_color_manual(name = "Approach", values = c("#33a02c", "#1f78b4")) +
  geom_point(position = position_dodge(0.8))+
  geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD), width=.2,
                position=position_dodge(0.8)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 

## That didn't really work...what about prediction curves? 
# Results
for(i in seq_along(res.files)){
  
  spp<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][1])
  season<- toupper(strsplit(gsub(".rds", "", gsub("mcmc_", "", res.files[i])), "_")[[1]][2])
  
  # Jon's data
  vuln<- paste(jon.df$SENSITIVITY[jon.df$COMNAME == spp], "_", jon.df$SENSITIVITY.CERTAINTY[jon.df$COMNAME == spp], sep = "")
  dir<- paste(jon.df$DIRECTIONAL.EFFECT[jon.df$COMNAME == spp], "_", jon.df$DIRECTIONAL.EFFECT.CERTAINTY[jon.df$COMNAME == spp], sep = "")
  plot.title<- paste("Sens: ", vuln, " and Dir: ", dir, sep = "")
  
  # Model fit -- presence and biomass
  gam.b.temp<- readRDS(paste(out.dir, gsub("mcmc_", "gamfitbio", res.files[i]), sep = ""))
  ilink<- family(gam.b.temp)$linkinv
  gam.coef<- names(coef(gam.b.temp))
  
  # Ranks of candidate draws...
  # Get ranks for SDM only and for NEVA
  likes.temp<- data.frame(read_rds(paste(out.dir, res.files[i], sep = ""))[[1]])
  names(likes.temp)<- c("Likelihood", "Prior", "Posterior")
  likes.temp$Iteration<- rep(seq(from = 1, to = nrow(likes.temp), by = 1))
  likes.temp$SDM.Rank<- order(likes.temp$Prior)
  likes.temp$NEVA.Rank<- order(likes.temp$Posterior)
  
  # Need to make predictions from each of these iterations
  mods.temp<- data.frame(read_rds(paste(out.dir, res.files[i], sep = ""))[[2]])
  colnames(mods.temp)<- gam.coef
  
  # Get the SDM fits...
  sdm.best<- likes.temp[likes.temp$SDM.Rank<= 900, ]
  sdm.best<- sdm.best[order(sdm.best$SDM.Rank), ]
  sdm.best.params<- mods.temp[sdm.best$Iteration, ] 
  
  # Get the NEVA predictions
  neva.best<- likes.temp[likes.temp$NEVA.Rank<= 900, ]
  neva.best<- neva.best[order(neva.best$NEVA.Rank), ]
  neva.best.params<- mods.temp[neva.best$Iteration, ] 
  
  # Plot the curves...
  # All possible prediction curves, the original smooth, and then the one we have selected...
  # Some prep
  pred.dat.use<- pred.dat[pred.dat$SEASON == season, ]
  rescaled.dat.use<- rescaled.dat[rescaled.dat$Season == season, ]
  
  # Predictor matrix
  pred.mat<- predict(gam.b.temp, newdata = pred.dat.use, type = "lpmatrix")
  
  # All predictions
  pred.all.sdm<- apply(sdm.best.params, 1, function(x) exp(ilink(pred.mat %*% x)))
  pred.all.neva<- apply(neva.best.params, 1, function(x) exp(ilink(pred.mat %*% x)))
  
  # Plotting
  # SST
  want<- 1:500
  sst.all.sdm<- data.frame(pred.all.sdm[want,])
  colnames(sst.all.sdm)<- paste("Mod.", seq(from = 1, to = ncol(sst.all.sdm)), sep = "")
  sst.all.sdm<- sst.all.sdm %>%
    gather(., "Model", "Pred")
  sst.all.sdm$Value<- rep(rescaled.dat.use$SST, length(unique(sst.all.sdm$Model)))
  sst.all.sdm$Parameter<- rep("SST", nrow(sst.all.sdm))
  
  sst.all.neva<- data.frame(pred.all.neva[want,])
  colnames(sst.all.neva)<- paste("Mod.", seq(from = 1, to = ncol(sst.all.neva)), sep = "")
  sst.all.neva<- sst.all.neva %>%
    gather(., "Model", "Pred")
  sst.all.neva$Value<- rep(rescaled.dat.use$SST, length(unique(sst.all.neva$Model)))
  sst.all.neva$Parameter<- rep("SST", nrow(sst.all.neva))
  
  # Mean, Max and Min at each SST value across all models...
  sst.summarized.sdm<- sst.all.sdm %>%
    group_by(., Value) %>%
    summarize_at(., "Pred", c("mean", "min", "max"), na.rm = T)
  
  sst.summarized.neva<- sst.all.neva %>%
    group_by(., Value) %>%
    summarize_at(., "Pred", c("mean", "min", "max"), na.rm = T)
  
  sst.dat<- data.frame("Parameter" = rep("SST", 1000), "Value" = rep(rescaled.dat.use$SST, 2), "Pred.Mean" = c(sst.summarized.sdm$mean, sst.summarized.neva$mean), "Pred.Min" = c(sst.summarized.sdm$min, sst.summarized.neva$min), "Pred.Max" = c(sst.summarized.sdm$max, sst.summarized.neva$max), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
  
  sst.out<- ggplot() +
    geom_line(data = sst.dat, aes(x = Value, y = log(Pred.Mean+1), group = Model, color = Model)) +
    scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
    xlab("SST") +
    theme_bw() 
  
  sst.out<- sst.out +
    geom_line(data = sst.dat, aes(x = Value, y = log(Pred.Min+1), group = Model, color = Model), lty = "dashed")
  
  sst.out<- sst.out +
    geom_line(data = sst.dat, aes(x = Value, y = log(Pred.Max+1), group = Model, color = Model), lty = "dashed") 
  
  # Depth
  want<- 501:1000
  depth.all.sdm<- data.frame(pred.all.sdm[want,])
  colnames(depth.all.sdm)<- paste("Mod.", seq(from = 1, to = ncol(depth.all.sdm)), sep = "")
  depth.all.sdm<- depth.all.sdm %>%
    gather(., "Model", "Pred")
  depth.all.sdm$Value<- rep(rescaled.dat.use$Depth, length(unique(depth.all.sdm$Model)))
  depth.all.sdm$Parameter<- rep("Depth", nrow(depth.all.sdm))
  
  depth.all.neva<- data.frame(pred.all.neva[want,])
  colnames(depth.all.neva)<- paste("Mod.", seq(from = 1, to = ncol(depth.all.neva)), sep = "")
  depth.all.neva<- depth.all.neva %>%
    gather(., "Model", "Pred")
  depth.all.neva$Value<- rep(rescaled.dat.use$Depth, length(unique(depth.all.neva$Model)))
  depth.all.neva$Parameter<- rep("Depth", nrow(depth.all.neva))
  
  # Mean, Max and Min at each depth value across all models...
  depth.summarized.sdm<- depth.all.sdm %>%
    group_by(., Value) %>%
    summarize_at(., "Pred", c("mean", "min", "max"), na.rm = T)
  
  depth.summarized.neva<- depth.all.neva %>%
    group_by(., Value) %>%
    summarize_at(., "Pred", c("mean", "min", "max"), na.rm = T)
  
  depth.dat<- data.frame("Parameter" = rep("Depth", 1000), "Value" = rep(rescaled.dat.use$Depth, 2), "Pred.Mean" = c(depth.summarized.sdm$mean, depth.summarized.neva$mean), "Pred.Min" = c(depth.summarized.sdm$min, depth.summarized.neva$min), "Pred.Max" = c(depth.summarized.sdm$max, depth.summarized.neva$max), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
  
  depth.out<- ggplot() +
    geom_line(data = depth.dat, aes(x = Value, y = log(Pred.Mean+1), group = Model, color = Model)) +
    scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
    xlab("Depth") +
    theme_bw() 
  
  depth.out<- depth.out +
    geom_line(data = depth.dat, aes(x = Value, y = log(Pred.Min+1), group = Model, color = Model), lty = "dashed")
  
  depth.out<- depth.out +
    geom_line(data = depth.dat, aes(x = Value, y = log(Pred.Max+1), group = Model, color = Model), lty = "dashed") 
  
  out<- plot_grid(sst.out + theme(legend.position="none"), depth.out + theme(legend.position="none"), nrow = 1, ncol = 2, align = "hv", scale = 1)
  legend<- get_legend(sst.out)
  out<- plot_grid(out, legend, rel_widths = c(3, 0.5))
  title<- ggdraw() + draw_label(plot.title, fontface='bold')
  out<- plot_grid(title, out, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  ggsave(paste(out.dir, tolower(spp), "_", tolower(season), "_PredCurvInt.jpg", sep = ""), out, width = 11, height = 8, dpi = 125, units = "in")
}

# Manuscript Figures — Study area map --------------------------------------
# Spatial projections
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19

nelme<- st_read(here::here("Data", "NELME_clipped.shp"))
st_crs(nelme)<- "+init=epsg:4326"
nelme.sp<- as(st_zm(nelme), "Spatial")

bstrat<- st_read(here::here("Data/BottomTrawlStrata", "BTS_Strata.shp")) %>%
  st_transform(., "+init=epsg:4326")
strata.ca<- c(1351, 1310, 1320, 1410, 1420, 1490, 1990, 1410, 1420, 1490, 5440, 5480, 5430) # Canada

bstrat<- bstrat %>%
  filter(., !STRATA %in% strata.ca & STRATA <= 3990) %>%
  filter(., !as.numeric(STRATUMA) %in% strata.ca & as.numeric(STRATUMA) <= 3990)

# GoM
gom<- st_read(here::here("Data", "GoM_sf.shp"))
st_crs(gom)<- "+init=epsg:4326"
gom.sp<- as(st_zm(gom), "Spatial")

# SNE/MAB
south<- st_read(here::here("Data", "SNEandMAB_sf.shp"))
st_crs(south)<- "+init=epsg:4326"
south<- as(st_zm(south), "Spatial")

# Spatial stuff -- gets us the states and shoreline

#Bounds
xlim.use<- c(-77, -65)
ylim.use<- c(35.05, 45.2)

states <- c("Maine", "New Hampshire", "Massachusetts", "Vermont", "New York", "Rhode Island", "Connecticut", "Delaware", "New Jersey", "Maryland", "Pennsylvania", "Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", "District of Columbia", "West Virgina")
provinces <- c("Ontario", "Québec", "Nova Scotia", "New Brunswick")

us<- raster::getData("GADM",country="USA",level=1)
us.states<- us[us$NAME_1 %in% states,]
us.states<- gSimplify(us.states, tol=0.01, topologyPreserve=TRUE)
us.states<- st_as_sf(us.states)
canada<- raster::getData("GADM",country="CAN",level=1)
ca.provinces<- canada[canada$NAME_1 %in% provinces,]
ca.provinces<- gSimplify(ca.provinces, tol=0.01, topologyPreserve=TRUE)
ca.provinces<- st_as_sf(ca.provinces)

us.states.f<- fortify(us.states, NAME_1)
ca.provinces.f<- fortify(ca.provinces, NAME_1)

# Alright, plot time
plot.out<- ggplot() + 
  geom_sf(data = us.states, fill = "white", lwd = 0.4, show.legend = FALSE) +
  geom_sf(data = ca.provinces.f, fill = "white", lwd = 0.4, show.legend = FALSE) +
  geom_sf(data = gom, aes(fill = "#377eb8"), color = NA, alpha = 0.75, show.legend = TRUE) +
  geom_sf(data = south, aes(fill = "#4daf4a"), color = NA, alpha = 0.75, show.legend = TRUE) +
  geom_sf(data = bstrat, fill = NA, color = "black", show.legend = FALSE) + 
  scale_fill_manual(name = "Region", values = c("#377eb8", "#ff7f00"), labels = c("Gulf of Maine", "Southern New England/Mid Atlantic Bight")) +
  xlim(xlim.use) +
  ylim(ylim.use) +
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"))
ggsave(paste(out.dir, "StudyArea.jpg", sep = ""), plot.out, width = 8, height = 11, dpi = 400)

us.large<- gSimplify(us, tol=0.1, topologyPreserve=TRUE)
us.large<- st_as_sf(us.large)

ca.large<- gSimplify(canada, tol=0.1, topologyPreserve=TRUE)
ca.large<- st_as_sf(ca.large)

plot.large<- ggplot() + 
  geom_sf(data = us.large, fill = "white", lwd = 0.7) +
  geom_sf(data = ca.large, fill = "white", lwd = 0.7) +
  xlim(c(-100, 0)) +
  ylim(c(25, 55)) +
  coord_sf(datum = NA) 
ggsave(paste(out.dir, "StudyAreaLarge.jpg", sep = ""), plot.large, width = 11, height = 8)

# Manuscript Figures — Combined approach figure ---------------------------
set.seed(1313) 
# Sensitivity/Vulnerability
plot <- ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1)) +
  ylab("Density") +
  xlab("Future Projected Relative Biomass\n Magnitude of Change") +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), labels = c("-2 SD", "-1 SD", "Mean", "+1 SD", "+2 SD")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14)) 

base.preds<- rnorm(101, mean = 0.1, sd = 1)
base.preds.probs<- round(quantile(base.preds, c(0, 0.2, 0.4, 0.6, 0.8, 1)), 2)

base.preds.labs<- data.frame("X" = as.numeric(base.preds.probs), "Y" = rep(0.5, length(base.preds.probs)), label = c("Baseline\n 0th pctl", "Baseline\n 20th pctl", "Baseline\n 40th pctl", "Baseline\n 60th pctl", "Baseline\n 80th pctl", "Baseline\n 100th pctl"))

plot2<- plot + 
  geom_vline(xintercept = base.preds.probs, lty = "dashed") +
  geom_label(data = base.preds.labs, aes(x = X, y = Y, label = label), size = 4)

exp.sens.out<- plot2 + stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), xlim = c(-3, base.preds.labs$X[1]), geom = "area", fill = '#d7301f', alpha = .75) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), xlim = c(base.preds.labs$X[6], 3), geom = "area", fill = '#d7301f', alpha = .75) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), xlim = c(base.preds.labs$X[1], base.preds.labs$X[2]), geom = "area", fill = '#fc8d59', alpha = .75) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), xlim = c(base.preds.labs$X[5], base.preds.labs$X[6]), geom = "area", fill = '#fc8d59', alpha = .75) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), xlim = c(base.preds.labs$X[2], base.preds.labs$X[3]), geom = "area", fill = '#fdcc8a', alpha = .75) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), xlim = c(base.preds.labs$X[4], base.preds.labs$X[5]), geom = "area", fill = '#fdcc8a', alpha = .75) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), xlim = c(base.preds.labs$X[3], base.preds.labs$X[4]), geom = "area", fill = 'white', alpha = .75) +
  geom_label(aes(x = -2.8, y = 0.45, label = "Very high"), size = 4, fill = '#d7301f', alpha = .6) +
  geom_label(aes(x = 3, y = 0.45, label = "Very high"), size = 4, fill = '#d7301f', alpha = .6) +
  geom_label(aes(x = -1.75, y = 0.45, label = "High"), size = 4, fill = '#fc8d59', alpha = .6) +
  geom_label(aes(x = 1.75, y = 0.45, label = "High"), size = 4, fill = "#fc8d59", alpha = .6) +
  geom_label(aes(x = -0.5, y = 0.45, label = "Moderate"), size = 4, fill = "#fdcc8a", alpha = .6) +
  geom_label(aes(x = 0.5, y = 0.45, label = "Moderate"), size = 4, fill = "#fdcc8a", alpha = .6) +
  geom_label(aes(x = 0, y = 0.45, label = "Low"), size = 4, fill = "white") 

# Dir Effect
plot <- ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1)) +
  ylab("Density") +
  xlab("Future Projected Relative Biomass\n") +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), labels = c("-2 SD", "-1 SD", "Mean", "+1 SD", "+2 SD")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14)) 

base.preds.probs<- round(quantile(base.preds, c(0.33, 0.67)), 2)

base.preds.labs<- data.frame("X" = as.numeric(base.preds.probs), "Y" = rep(0.5, length(base.preds.probs)), label = c("Baseline\n 33rd pctl", "Baseline\n 67th pctl"))

plot2<- plot + 
  geom_vline(xintercept = base.preds.probs, lty = "dashed") +
  geom_label(data = base.preds.labs, aes(x = X, y = Y, label = label), size = 4)

dir.out<- plot2 + stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), xlim = c(-3, base.preds.labs$X[1]), geom = "area", fill = '#0571b0', alpha = .75) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), xlim = c(base.preds.labs$X[1], base.preds.labs$X[2]), geom = "area", fill = 'white', alpha = .75) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), xlim = c(base.preds.labs$X[2], 3), geom = "area", fill = '#ca0020', alpha = .75) +
  geom_label(aes(x = -1.8, y = 0.45, label = "Negative"), size = 4, fill = '#0571b0', alpha = .6) +
  geom_label(aes(x = 0.05, y = 0.45, label = "Neutral"), size = 4, fill = 'white', alpha = .6) +
  geom_label(aes(x = 1.8, y = 0.45, label = "Positive"), size = 4, fill = '#ca0020', alpha = .6)

out<- plot_grid(exp.sens.out, dir.out, ncol = 1)
ggsave(paste(out.dir, "Distributions.jpg", sep = ""), out, width = 15, height = 15, dpi = 300)


# Manuscript Figures — Climate model projections --------------------------------
# What would we like to show -- sea surface temperature projected temperature changes relative to 1982-2011 and then temperature changes relative to baseline period 2011-2015
## Some spatial stuff for data visualiztion
# Spatial projections
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19

#Bounds
xlim.use<- c(-77, -65)
ylim.use<- c(35.05, 45.2)

states <- c("Maine", "New Hampshire", "Massachusetts", "Vermont", "New York", "Rhode Island", "Connecticut", "Delaware", "New Jersey", "Maryland", "Pennsylvania", "Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", "District of Columbia", "West Virgina")
provinces <- c("Ontario", "Québec", "Nova Scotia", "New Brunswick")

us <- raster::getData("GADM",country="USA",level=1)
us.states <- us[us$NAME_1 %in% states,]
us.states <- gSimplify(us.states, tol = 0.075, topologyPreserve = TRUE)
canada <- raster::getData("GADM",country="CAN",level=1)
ca.provinces <- canada[canada$NAME_1 %in% provinces,]
ca.provinces <- gSimplify(ca.provinces, tol = 0.075, topologyPreserve = TRUE)

us.states.f<- fortify(us.states, NAME_1)
ca.provinces.f<- fortify(ca.provinces, NAME_1)

### Applying anomalies -- mean
## Inspect the file using the ncdf4 libaray
## Get sst anomaly 
sst.anom.temp<- raster::stack(here::here("Data", "SST.CMIP5.1982-2099.anom.nc"), varname = "sstanom")
sst.anom<- raster::rotate(sst.anom.temp)

## Get oisst data
oisst.dat.temp<- raster::stack(here::here("Data", "EC_sst_1981_2015_OISST-V2-AVHRR_agg_combined.nc"))
oisst.dat<- raster::rotate(oisst.dat.temp)

# Need to get climatology from the OISST data -- set up OISST stack as time series
oisst.min<- gsub("X", "", min(names(oisst.dat)))
oisst.min.date<- as.Date(gsub("[.]", "-", oisst.min))
oisst.max<- gsub("X", "", max(names(oisst.dat)))
oisst.max.date<- as.Date(gsub("[.]", "-", oisst.max))

# Calculate monthly mean temperature -- this would be compared to the sstclim data (monthly climate ensemble)
oisst.dates<- seq.Date(from = oisst.min.date, to = oisst.max.date, by = "day")
oisst.dat<- setZ(oisst.dat, oisst.dates)

# Aggregate daily to monthly data
oisst.monthly <- zApply(oisst.dat, by = as.yearmon, mean)

# Lets check that
test.dat<- subset(oisst.dat, which(getZ(oisst.dat) >="1981-09-01" & getZ(oisst.dat) <= "1981-09-30"))
sept1981.mu<- calc(test.dat, mean)
plot(oisst.monthly[[1]]-sept1981.mu)

# Everything seems fine there, now need the monthly average for each month across baseline years (1982-2011)
dates<- getZ(oisst.monthly)
subset.vec<- which(dates > "Dec 1981" & dates < "Jan 2012", arr.ind = TRUE)
oisst.monthly.sub<- oisst.monthly[[subset.vec]]
oisst.monthly.sub<- setZ(oisst.monthly.sub, dates[subset.vec])

oisst.clim<- stack()
months<- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

for(i in seq_along(months)) {
  # Get all the monthly raster
  stack.temp<- subset(oisst.monthly.sub, which(grepl(months[i], names(oisst.monthly.sub))))
  month.clim<- calc(stack.temp, fun = mean)
  oisst.clim<- stack(oisst.clim, month.clim)
  names(oisst.clim)[i]<- months[i]
}

# Check that
test.dat<- subset(oisst.monthly.sub, which(grepl("Jan", names(oisst.monthly.sub))))
jan.mu<- calc(test.dat, mean)
summary(oisst.clim[[1]] - jan.mu)

# Looks good -- time to apply the anomalies to the oisst.clim
oisst.clim.coarse<- raster::resample(oisst.clim, sst.anom)
names(oisst.clim.coarse)<- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# Okay, now good to apply the anomalies from the climate models to climatology and get raw values
sst.model<- stack()

for(i in 1:nlayers(sst.anom)) {
  index.match<- which(gsub("X", "", names(oisst.clim.coarse)) == unlist(strsplit(names(sst.anom)[i], "[.]"))[2])
  rast.temp<- oisst.clim.coarse[[index.match]] + sst.anom[[i]]
  sst.model<- stack(sst.model, rast.temp)
  names(sst.model)[i]<- names(sst.anom)[i]
}

# One more step, need to fill in the missing coastline raster cell values.... Function below, i is corresponding to a three by three window
fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),0) )
  } else {
    return( round(x[i],0) )
  }
}  

# Now apply that function to each raster stack
for(i in 1:nlayers(sst.model)) {
  new.rast<- focal(sst.model[[i]], w = matrix(1, 3, 3), fun = fill.na, pad = TRUE, na.rm = FALSE)
  sst.model[[i]]<- new.rast
}

sst.model.proj<- projectRaster(sst.model, crs = proj.utm)
names(sst.model.proj)<- names(sst.anom)
sst.model<- projectRaster(sst.model.proj, crs = proj.wgs84)

# Okay....now can we get some season and regional averages? Would want to show NES LME, and then GoM vs. SNE changes in spring and fall...
nelme.sst<- raster::extract(sst.model, nelme, fun = mean)
south.sst<- raster::extract(sst.model, south, fun = mean)
gom.sst<- raster::extract(sst.model, gom, fun = mean)

gom.oisst<- raster::extract(oisst.dat, gom, fun = mean, na.rm = T)
south.oisst<- raster::extract(oisst.dat, south, fun = mean, na.rm = T)
reg.oisst.df<- data.frame("Date" = rep(gsub("X", "", dimnames(gom.oisst)[[2]]), 2), "Region" = c(rep("Gulf of Maine", length(dimnames(gom.oisst)[[2]])), rep("Southern New England", length(dimnames(gom.oisst)[[2]]))), "SST" = c(as.numeric(gom.oisst), as.numeric(south.oisst)))
reg.oisst.df$Date<- as.Date(gsub("[.]", "-", reg.oisst.df$Date))
reg.oisst.df$Month<- format(reg.oisst.df$Date, "%m")
reg.oisst.df$Year<- as.numeric(format(reg.oisst.df$Date, "%Y"))
reg.oisst.df.summ<- reg.oisst.df %>%
  group_by(., Region, Year, Month) %>%
  summarise("DailyTempAboveThresh" = sum(SST>20, na.rm = T),
            "DaysMonth" = n()) %>%
  mutate("Proportion" = DailyTempAboveThresh/DaysMonth)
reg.oisst.df.summ$Plot.Date<- as.Date(as.yearmon(paste(reg.oisst.df.summ$Year, reg.oisst.df.summ$Month, sep = "-")))

plot.dat<- reg.oisst.df.summ %>%
  dplyr::filter(., Plot.Date >= "1982-01-01" & Region == "Gulf of Maine")
ggplot(data = plot.dat, aes(x = Plot.Date, y = Proportion, group = Region)) +
  geom_line()


# Now --- want to get the average value by date and plot the time series
sst.mn<- data.frame("Date" = rep(gsub("X", "", dimnames(nelme.sst)[[2]]), 3), "Region" = c(rep("NES LME", length(dimnames(nelme.sst)[[2]])), rep("GoM", length(dimnames(nelme.sst)[[2]])), rep("Southern NES LME", length(dimnames(nelme.sst)[[2]]))), "SST" = c(as.numeric(nelme.sst), as.numeric(gom.sst), as.numeric(south.sst)))
sst.mn$Date<- as.Date(gsub("[.]", "-", sst.mn$Date))
sst.mn$Month<- as.numeric(format(sst.mn$Date, "%m"))
sst.mn$Year<- as.numeric(format(sst.mn$Date, "%Y"))
sst.mn$Season<- ifelse(sst.mn$Month >= 3 & sst.mn$Month <= 5, "SPRING", 
                       ifelse(sst.mn$Month >= 9 & sst.mn$Month <= 11, "FALL", NA))


## Percent of months that exceed 20 degrees
sst.mn.20deg<- sst.mn %>%
  mutate(., "Threshold" = ifelse(SST >= 20, "Yes", "No"))
sst.mn.20deg<- sst.mn.20deg %>%
  group_by(., Region, Year, Threshold) %>%
  tally() %>%
  filter(., Threshold == "Yes")

ggplot(data = sst.mn.20deg, aes(x = Year, y = n, color = Region)) +
  geom_line() 


## Inspect the file using the ncdf4 libaray
## Get sst anomaly 
sst.anom.temp<- raster::stack(here::here("Data", "SST.CMIP5.1982-2099.anom.nc"), varname = "sstpct05")
sst.anom<- raster::rotate(sst.anom.temp)

# Apply the new anomalies to the oisst.clim
oisst.clim.coarse<- raster::resample(oisst.clim, sst.anom)
names(oisst.clim.coarse)<- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# Okay, now good to apply the anomalies from the climate models to climatology and get raw values
sst.model<- stack()

for(i in 1:nlayers(sst.anom)) {
  index.match<- which(gsub("X", "", names(oisst.clim.coarse)) == unlist(strsplit(names(sst.anom)[i], "[.]"))[2])
  rast.temp<- oisst.clim.coarse[[index.match]] + sst.anom[[i]]
  sst.model<- stack(sst.model, rast.temp)
  names(sst.model)[i]<- names(sst.anom)[i]
}

# One more step, need to fill in the missing coastline raster cell values.... Function below, i is corresponding to a three by three window
fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),0) )
  } else {
    return( round(x[i],0) )
  }
}  

# Now apply that function to each raster stack
for(i in 1:nlayers(sst.model)) {
  new.rast<- focal(sst.model[[i]], w = matrix(1, 3, 3), fun = fill.na, pad = TRUE, na.rm = FALSE)
  sst.model[[i]]<- new.rast
}

sst.model.proj<- projectRaster(sst.model, crs = proj.utm)
names(sst.model.proj)<- names(sst.anom)
sst.model<- projectRaster(sst.model.proj, crs = proj.wgs84)

# Okay....now can we get some season and regional averages? Would want to show NES LME, and then GoM vs. SNE changes in spring and fall...
nelme.sst<- raster::extract(sst.model, nelme, fun = mean)
south.sst<- raster::extract(sst.model, south, fun = mean)
gom.sst<- raster::extract(sst.model, gom, fun = mean)

# Now --- want to get the average value by date and plot the time series
sst.pct05<- data.frame("Date" = rep(gsub("X", "", dimnames(nelme.sst)[[2]]), 3), "Region" = c(rep("NES LME", length(dimnames(nelme.sst)[[2]])), rep("GoM", length(dimnames(nelme.sst)[[2]])), rep("Southern NES LME", length(dimnames(nelme.sst)[[2]]))), "SST.05" = c(as.numeric(nelme.sst), as.numeric(gom.sst), as.numeric(south.sst)))
sst.pct05$Date<- as.Date(gsub("[.]", "-", sst.pct05$Date))
sst.pct05$Month<- as.numeric(format(sst.pct05$Date, "%m"))
sst.pct05$Year<- as.numeric(format(sst.pct05$Date, "%Y"))
sst.pct05$Season<- ifelse(sst.pct05$Month >= 3 & sst.pct05$Month <= 5, "SPRING", 
                          ifelse(sst.pct05$Month >= 9 & sst.pct05$Month <= 11, "FALL", NA))

## Inspect the file using the ncdf4 libaray
## Get sst anomaly 
sst.anom.temp<- raster::stack(here::here("Data", "SST.CMIP5.1982-2099.anom.nc"), varname = "sstpct95")
sst.anom<- raster::rotate(sst.anom.temp)

# Apply final anomalies to the oisst.clim
oisst.clim.coarse<- raster::resample(oisst.clim, sst.anom)
names(oisst.clim.coarse)<- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# Okay, now good to apply the anomalies from the climate models to climatology and get raw values
sst.model<- stack()

for(i in 1:nlayers(sst.anom)) {
  index.match<- which(gsub("X", "", names(oisst.clim.coarse)) == unlist(strsplit(names(sst.anom)[i], "[.]"))[2])
  rast.temp<- oisst.clim.coarse[[index.match]] + sst.anom[[i]]
  sst.model<- stack(sst.model, rast.temp)
  names(sst.model)[i]<- names(sst.anom)[i]
}

# One more step, need to fill in the missing coastline raster cell values.... Function below, i is corresponding to a three by three window
fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),0) )
  } else {
    return( round(x[i],0) )
  }
}  

# Now apply that function to each raster stack
for(i in 1:nlayers(sst.model)) {
  new.rast<- focal(sst.model[[i]], w = matrix(1, 3, 3), fun = fill.na, pad = TRUE, na.rm = FALSE)
  sst.model[[i]]<- new.rast
}

sst.model.proj<- projectRaster(sst.model, crs = proj.utm)
names(sst.model.proj)<- names(sst.anom)
sst.model<- projectRaster(sst.model.proj, crs = proj.wgs84)

# Okay....now can we get some season and regional averages? Would want to show NES LME, and then GoM vs. SNE changes in spring and fall...
nelme.sst<- raster::extract(sst.model, nelme, fun = mean)
south.sst<- raster::extract(sst.model, south, fun = mean)
gom.sst<- raster::extract(sst.model, gom, fun = mean)

# Now --- want to get the average value by date and plot the time series
sst.pct95<- data.frame("Date" = rep(gsub("X", "", dimnames(nelme.sst)[[2]]), 3), "Region" = c(rep("NES LME", length(dimnames(nelme.sst)[[2]])), rep("GoM", length(dimnames(nelme.sst)[[2]])), rep("Southern NES LME", length(dimnames(nelme.sst)[[2]]))), "SST.95" = c(as.numeric(nelme.sst), as.numeric(gom.sst), as.numeric(south.sst)))
sst.pct95$Date<- as.Date(gsub("[.]", "-", sst.pct95$Date))
sst.pct95$Month<- as.numeric(format(sst.pct95$Date, "%m"))
sst.pct95$Year<- as.numeric(format(sst.pct95$Date, "%Y"))
sst.pct95$Season<- ifelse(sst.pct95$Month >= 3 & sst.pct95$Month <= 5, "SPRING", 
                          ifelse(sst.pct95$Month >= 9 & sst.pct95$Month <= 11, "FALL", NA))
colnames(sst.pct95)[3]<- "SST.95"

# Merge em all together?
sst.all<- sst.mn %>%
  left_join(., sst.pct05) %>%
  left_join(., sst.pct95) %>%
  drop_na(Season)

# Fit linear model for trends
# linear model function
lm.fun <- function(df){
  df.mod<- df %>%
    group_by(Year) %>%
    summarize_at("SST", mean)
  mod.out<- lm(SST ~ Year, data = df.mod)
  return(mod.out)
}

sst.all.mod<- sst.all %>%
  group_by(Season, Region) %>%
  nest() %>%
  mutate(., "Linear.Model" = map(data, lm.fun),
         "Tidied" = map(Linear.Model, tidy),
         "Glanced" = map(Linear.Model, glance),
         "Augmented" = map(Linear.Model, augment))

# Fits
sst.all.mod %>%
  unnest(Tidied)

sst.all.mod %>%
  unnest(Glanced)

# Diffs to Kleisner
kleisner<- data.frame("Season" = c("FALL", "SPRING", "FALL", "SPRING", "FALL", "SPRING"), "Region" = c("NES LME", "NES LME", "GoM", "GoM", "Southern NES LME", "Southern NES LME"), "term" = rep("Year", 6), "estimate.k" = c(0.049, 0.048, 0.053, 0.052, 0.047, 0.047))
us<- sst.all.mod %>%
  unnest(Tidied)
us<- us[rep(c(FALSE, TRUE), 6),]
sst.comparison<- us %>%
  left_join(., kleisner) %>%
  mutate(., "SST.Trend.Diff" = estimate.k - estimate)

summary(sst.comparison$SST.Trend.Diff)
# Plot
sst.yrly<- sst.all %>%
  group_by(Season, Region, Year) %>%
  summarize_at(c("SST", "SST.05", "SST.95"), mean)

sst.base<- sst.yrly %>%
  filter(., Year >= 1982 & Year <= 2011) %>%
  group_by(Season, Region) %>%
  summarize_at(c("SST"), mean)
colnames(sst.base)[3]<- "SST.Base"

sst.plot<- sst.yrly %>%
  left_join(., sst.base)
sst.plot<- sst.plot %>%
  mutate(., "Anomaly.Mean" = SST - SST.Base,
         "Anomaly.5th" = SST.05 - SST.Base,
         "Anomaly.95th" = SST.95 - SST.Base)
sst.plot$Region<- factor(sst.plot$Region, levels = c("NES LME", "GoM", "Southern NES LME"), labels = c("Northeast Shelf Large Marine Ecosystem", "Gulf of Maine", "Southern New England/Mid Atlantic Bight"))
sst.plot$Season<- factor(sst.plot$Season, levels = c("FALL", "SPRING"), labels = c("Fall", "Spring"))

sst.plot<- sst.plot %>%
  dplyr::filter(., Year <= 2055)

clim.sst.plot<- ggplot() + 
  #geom_vline(xintercept = 2050, size = 2) +
  geom_line(data = sst.plot, aes(x = Year, y = Anomaly.Mean, color = Region)) +
  geom_ribbon(data = sst.plot, aes(x = Year, ymin = Anomaly.5th, ymax = Anomaly.95th, fill = Region), alpha=0.15) +
  scale_color_manual(name = "Region", values = c("black", '#377eb8', '#4daf4a')) +
  scale_fill_manual(name = "Region", values = c("black", '#377eb8', '#4daf4a')) +
  labs(x = "Year", y = "RCP 8.5 Models Ensemble Projected SST Anomaly\n Relative to 21982-2011 Climatology", fill = "Region") +
  facet_wrap(~Season) +
  theme(strip.background = element_blank())
ggsave(paste(out.dir, "ClimateModelSSTAnomaly.jpg", sep = ""), clim.sst.plot, width = 11, height = 8)

# Add observed temperature anomalies?
oisst.dat.temp<- raster::stack(here::here("Data", "EC_sst_1981_2015_OISST-V2-AVHRR_agg_combined.nc"))
oisst.dat<- raster::rotate(oisst.dat.temp)

# Need to get climatology from the OISST data -- set up OISST stack as time series
oisst.min<- gsub("X", "", min(names(oisst.dat)))
oisst.min.date<- as.Date(gsub("[.]", "-", oisst.min))
oisst.max<- gsub("X", "", max(names(oisst.dat)))
oisst.max.date<- as.Date(gsub("[.]", "-", oisst.max))

# Libraries
library(rts)

# Preliminaries
month.conv<- data.frame("Month.Chr" = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"), "Month.Num" = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"))

baseline.period = c("1982-01-01", "2011-12-31")
oisst.dates<- seq.Date(from = oisst.min.date, to = oisst.max.date, by = "day")
oisst.dat<- setZ(oisst.dat, oisst.dates)

## Baseline
## Need baseline daily average temp across years for every day
dates.unique<- unique(format(getZ(oisst.dat), "%m-%d"))
daily.means<- stack(lapply(seq(length(dates.unique)), function(x) calc(oisst.dat[[which(grepl(dates.unique[x], getZ(oisst.dat)) == TRUE)]], fun = mean)))
names(daily.means)<- dates.unique
daily.sd<- stack(lapply(seq(length(dates.unique)), function(x) calc(oisst.dat[[which(grepl(dates.unique[x], getZ(oisst.dat)) == TRUE)]], fun = sd)))
names(daily.sd)<- dates.unique

# Alright, now substract each day from the average to get the anomaly
daily.anoms<- stack(lapply(seq(1:nlayers(oisst.dat)), function(x) (oisst.dat[[x]] - daily.means[[match(format(getZ(oisst.dat)[x], "%m-%d"), gsub("[.]", "-", gsub("X", "", names(daily.means))))]])/daily.sd[[match(format(getZ(oisst.dat)[x], "%m-%d"), gsub("[.]", "-", gsub("X", "", names(daily.sd))))]]))

## Get the data ready for plotting
# Smoothing
daily.anoms<- setZ(daily.anoms, oisst.dates)
names(daily.anoms)<- oisst.dates
ts.wide.daily<- do.call("cbind", lapply(seq(1:nlayers(daily.anoms)), function(x) as.data.frame(daily.anoms[[x]], xy = TRUE)))
ts.df.daily<- ts.wide.daily %>%
  subset(., select=which(!duplicated(names(.)))) %>%
  gather(., Year, SST, -x, -y) 

ts.df.daily$Year<- gsub("[.]", "-", gsub("X", "", ts.df.daily$Year))
ts.df.dailymu<- ts.df.daily %>%
  group_by(Year) %>%
  summarize(., Mean.SST = mean(SST, na.rm = T)) %>%
  separate(., col = Year, into = c("Year", "Month", "Day"), sep = "-") %>%
  mutate(., Plot.Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>%
  data.frame

ts.df.dailymu$Season<- ifelse(ts.df.dailymu$Month %in% c("03", "04", "05"), "Spring",
                              ifelse(ts.df.dailymu$Month %in% c("09", "10", "11"), "Fall", NA))

ts.df.yearlymu<- ts.df.dailymu %>%
  group_by(Year, Season) %>%
  summarize(., 
            "Anomaly.Mean" = mean(Mean.SST, na.rm = T)) %>%
  drop_na(Season)
ts.df.yearlymu$Season<- factor(ts.df.yearlymu$Season, levels = c("Fall", "Spring"))
ts.df.yearlymu$Year<- as.numeric(ts.df.yearlymu$Year)

clim.sst.plot.out<- clim.sst.plot +
  geom_point(data = ts.df.yearlymu, aes(x = Year, y = Anomaly.Mean), pch = 21, color = "black", size = 1) +
  geom_line(data = ts.df.yearlymu, aes(x = Year, y = Anomaly.Mean), linetype = "dotted")
ggsave(paste(out.dir, "ClimateModelSSTAnomaly.jpg", sep = ""), clim.sst.plot.out, width = 11, height = 8)


# Manuscript Figures — Lobster used temperatures --------------------------
dat.lob<- dat %>%
  filter(., COMNAME == "AMERICAN LOBSTER")

dat.lob.k<- dat.lob %>%
  filter(., EST_YEAR >= 1982 & EST_YEAR <= 2013)
dat.lob.k$SEASON<- rep("Combined", nrow(dat.lob.k))

dat.lob.us<- dat.lob %>%
  filter(., EST_YEAR >= 1982 & EST_YEAR <= 2011)

lob.us<- ggplot(dat.lob.us) +
  geom_density(aes(x = SEASONALMU.OISST, group = SEASON, linetype = SEASON)) +
  geom_hline(yintercept = 0, color = "white", lwd = 1) +
  ylab("Density") + 
  xlab(expression("Seasonal sea surface temperature ("*{}^{o}*"C)")) +
  scale_linetype_manual(name = "Season", values = c("dashed", "dotted"))

lob.both<- lob.us +
  geom_density(data = dat.lob.k, aes(x = SEASONALMU.OISST, linetype = SEASON)) +
  geom_hline(yintercept = 0, color = "white", lwd = 1) +
  scale_linetype_manual(name = "Season", values = c("solid", "dotted", "dashed"))

sst.out<- lob.both

lob.us<- ggplot(dat.lob.us) +
  geom_density(aes(x = BOTTEMP, group = SEASON, linetype = SEASON)) +
  geom_hline(yintercept = 0, color = "white", lwd = 1) +
  ylab("Density") + 
  xlab(expression("In situ bottom temperature ("*{}^{o}*"C)")) +
  scale_linetype_manual(name = "Season", values = c("dashed", "dotted"))

lob.both<- lob.us +
  geom_density(data = dat.lob.k, aes(x = BOTTEMP, linetype = SEASON)) +
  geom_hline(yintercept = 0, color = "white", lwd = 1) +
  scale_linetype_manual(name = "Season", values = c("solid", "dotted", "dashed"))

bot.out<- lob.both

out<- plot_grid(sst.out, bot.out, labels = c("American lobster used\n sea surface temperatures", "American lobster used\n bottom temperatures"), nrow = 1)
ggsave(paste(out.dir, "LobsterTemperaturePlots.jpg", sep = ""), out, width = 11, height = 8, units = "in")


  
