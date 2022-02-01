source("./R/dukem.R")

#############################################
#### TEMPORAL PROFILING ANNUAL EMISSIONS ####
#############################################

## See ReadMe for more detail ##
## Also create a total emissions annual profile, for comparison with Sentinel-5p.   ##

##########################################################################################
#### Create Profile_ID level GAMs from the formatted raw data, with accompanying csv: ####
GAMbyProfile(timescale = "hourwday", exclude = c("STAT_DOM_GEN","STAT_PIC_GEN"))

# plot the Profile_ID GAMs (does all of them):
GAMprofilePlot()
### The above only needs doing from new raw data to generate Profile GAMs for sectors ####
##########################################################################################

##########################################################################################
#### set the emissions year to model and the species ####
y_emis <- NA # the emissions year: SET TO NA for generic year representation
species <- NA # species to profiles: SET TO NA for generic pollutant representation
classification <- "SNAP" # this can be anything, as long as the lookup file "NFR_to_xxxxx.csv" exists
y_spec_NFR <- NULL # the NFR codes that are to have year-specific profiles (otherwise `average` year = 0)

v_time <- c("yday","month","wday","hour","hourwday") # currently yday, month, wday, hourwday and hour

#### Create new sector-wide temporal profiles: ####
# (needs to build weighted data and perform GAM routine over all sectors/time 

##!!## Need to decide on one GAM per timescale, with by= sector, or list of sector GAMs per timescale
##!!## The former is too computational for a laptop. 
        # option 1: do on modeling PC and move back to this Git. 
        # option 2: move everything to JASMIN
        # option 3: set up a dataLabs via the DUKEMs project. 

# GAMs by sector in loop and saved to list object - this is the only thing that will work on a laptop. 
lapply(v_time, GAMBYSectorLOOP,year = y_emis, species = species, classification = classification, yr_spec_NFR = NULL)
#GAMBYSectorLOOP(year = y_emis, species = species, timescale = "hour", classification = classification, yr_spec_NFR = NULL)


# plot the new sector-level GAMs:
GAMsectorPlot(year = y_emis, species = species, classification = classification)


###
lapply(v_time, GAMProfileBySector,year = y_emis, species = species, classification = classification, yr_spec_NFR = NULL)

