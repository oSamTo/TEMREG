source("./R/workspace.R")

#############################################
#### TEMPORAL PROFILING ANNUAL EMISSIONS ####
#############################################

## Taking NAEI annual totals for a species and splitting into hourly fractions.     ##
## This is done with the MOY/DOW/HOD profiles. Match a profile to an NFR and split. ##
## Then sum up the NFRs into SNAPs and create a temporal profile fo the SNAP        ##
## Also create a total emissions annual profile, for comparison with Sentinel-5p.   ##

## set the year to model and the species
y_emis <- 2019 # the emissions year
## !! Need a way of choosing NFR/fuel combos to be modelled on specific year data, as opposed to means !! ##
y_spec_NFR <- NULL # the NFR codes that are to have year-specific profiles (otherwise `average` year = 0)
species <- "NOx"
classification <- "SNAP"

v_num <<- "1B" # alphanumeric version ID, to carry through data/plots/gams outputs (max 6 char)

#### read and match the NAEI data to the temporal profile in the lookup table, along with SNAP sector ####
dt_naei_profs <- JoinNAEItoProfiles(year = y_emis, species = species) 
  
#### create new sector-wide temporal profiles: ####
## weighted emissions from all NFR/profile_ID combinations, per sector. recreate temporal profiles. 
## also outputs .csvs of coefficients and plot. 
v_time <- c("yday","month","wday","hour","hourwday")

GAMofGAMs(year = y_emis, species = species, classification = classification, emis = dt_naei_profs)


l_DUKEMs_profiles <- lapply(X=setNames(v_time, v_time), FUN = TempProfileBySector, year = y_emis, species = species, classification = classification, emis = dt_naei_profs, yr_spec_NFR = y_spec_NFR)

## using the new sector-level profiles, visualize emissions (also outputs some metadata)
PlotEmissionsOverTime(year = y_emis, species = species, classification = classification, sec_profs = l_DUKEMs_profiles, emis = dt_naei_profs)


####################################################################################################
