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

#### read and match the NAEI data to the temporal profile in the lookup table, along with SNAP sector ####
dt_naei_profs <- JoinNAEItoProfiles(year = y_emis, species = species) 
  
#### create sector-wide temporally profiled emissions: weighted emissions from all NFR/profile combinations, per sector ####
## choose either integers for SNAP, or capital letters for GNFR
v_sectors <- list(1,2,3,4,5,6,7,8,9,10,11)

l_DUKEMs_profiles <- lapply(X=setNames(v_sectors, v_sectors), FUN = EmissionsProfileBySector, year = y_emis, species = species, classification = classification, emis = dt_naei_profs, yr_spec_NFR=NULL, hod_by_dow=F, hour_emis=T)

saveRDS(l_DUKEMs_profiles, paste0("./doc/profiled_emissions_sector/",species,"_",classification,"_",y_emis,"_profemis.rds"))

#### Using temporally distributed emissions (made above); ####
## format the data to standard output tables: yday, month, wday, hour and hour by wday 
FormatToSectorCoeffs(year = y_emis, species = species, classification = classification)

## create another fitted profile per sector (GAM) based on the weighted emissions
FitGAMsToSectorTotals()


####################################################################################################



