source("./R/dukem.R")

#############################################
#### TEMPORAL PROFILING ANNUAL EMISSIONS ####
#############################################

## See ReadMe for more detail ##
## Also create a total emissions annual profile, for comparison with Sentinel-5p.   ##

##########################################################################################
#### Create Profile_ID level GAMs from the formatted raw data, with accompanying csv: ####
GAMbyProfile(timescale = "hourwday", exclude = c("STAT_DOM_GEN","STAT_PIC_GEN"))

# plot the Profile_ID GAMs (plots Profiles that are passed in v_profiles, this is in parallel):
GAMparallelPlot(v_profiles = unique(NFR_to_Profile[,Profile_ID]))
### The above only needs doing from new raw data to generate Profile GAMs for sectors ####
##########################################################################################

##########################################################################################
#### set the emissions year to model and the species ####
y_emis <- NA # the emissions year: SET TO NA for generic year representation
species <- "NOX" # species to profiles: SET TO NA for generic pollutant representation
classification <- "SNAP" # this can be anything, as long as the lookup file "NFR_to_xxxxx.csv" exists
y_spec_NFR <- NULL # the NFR codes that are to have year-specific profiles (otherwise `average` year = 0)

v_time <- c("yday","month","wday","hour","hourwday") # currently yday, month, wday, hourwday and hour

#### Create new sector-wide temporal profiles: ####
# GAMs by sector in loop and saved to list object - currently this is the method being used. 
lapply(v_time, GAMBYSectorLOOP, year = y_emis, species = species, classification = classification, yr_spec_NFR = NULL)

# plot the new sector-level GAMs:
GAMsectorPlot(year = y_emis, species = species, classification = classification)
