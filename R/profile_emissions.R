source("./R/workspace.R")

#############################################
#### TEMPORAL PROFILING ANNUAL EMISSIONS ####
#############################################

## See ReadMe for more detail ##
## Also create a total emissions annual profile, for comparison with Sentinel-5p.   ##

## set the emissions year to model and the species
y_emis <- 2019 # the emissions year
## !! Need a way of choosing NFR/fuel combos to be modelled on specific year data, as opposed to means !! ##
y_spec_NFR <- NULL # the NFR codes that are to have year-specific profiles (otherwise `average` year = 0)
species <- "NOx"
classification <- "SNAP" # this can be anything, as long as the file exists

#v_num <<- "1B" # alphanumeric version ID, to carry through data/plots/gams outputs (max 6 char)

#### read and match the NAEI data to the temporal profile in the lookup table, along with SNAP sector ####
dt_naei_profs <- JoinNAEItoProfiles(year = y_emis, species, classification) 

#### create new sector-wide temporal profiles: ####
# (needs to build weighted data and perform GAM routine over all sectors/time - ~6 mins for 'SNAP' for one pollutant)
v_time <- c("yday","month","wday","hour","hourwday") # currently yday, month, wday, hourwday and hour

lapply(X=setNames(v_time, v_time), FUN = GAMProfileBySector, year = y_emis, species = species, classification = classification, emis = dt_naei_profs, yr_spec_NFR = y_spec_NFR)

#### create coefficients (centered on 1) as csv files for sectors: ####
v_time <- c("yday","month","wday","hour","hourwday") # currently yday, month, wday, hourwday and hour

lapply(X=setNames(v_time, v_time), FUN = sectorCoefficients, year = y_emis, species = species, classification = classification, emis = dt_naei_profs)

#### plot GAMs for any species, year & timesteps on the same plot ####
## the same is done for emissions.
GAMplots(v_years = c(2019), v_species = c("NOx"), classification = "SNAP")

####################################################################################################