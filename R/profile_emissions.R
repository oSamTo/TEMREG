source("./R/re-structure.R")

#############################################
#### TEMPORAL PROFILING ANNUAL EMISSIONS ####
#############################################

## See ReadMe for more detail ##
## Also create a total emissions annual profile, for comparison with Sentinel-5p.   ##

#### Create Profile_ID level GAMs from the formatted raw data: ####
GAMbyProfile(timescale = "hourwday", exclude = c("STAT_DOM_GEN","STAT_PIC_GEN"))

# plot the Profile_ID GAMs (does all of them):
GAMprofilePlot()

## set the emissions year to model and the species
y_emis <- c(2019) # the emissions year: SET TO NA for generic year representation
species <- "NOx" # species to profiles: SET TO NA for generic pollutant representation
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


lapply(v_time, GAMProfileBySector,year = y_emis, species = species, classification = classification, yr_spec_NFR = NULL)

GAMBySector(year = y_emis, species = species, timescale="yday", classification = classification, yr_spec_NFR = NULL)

# plot the new sector-level GAMs:
GAMsectorPlot(year = y_emis, species = species, timescale = "hour", classification = classification)





#### plot GAMs for any species, year & timesteps on the same plot ####
## the same is done for emissions.
GAMplots(year = c(2019), v_species = c("NOx","CO2","CH4","N2O","NH3","SOx"), classification = "SNAP")

####################################################################################################