library(here) # construct file paths relative to project root
here::i_am("./_targets.R")

##########################################
#### Processing the UK inventory data ####
##########################################

#####################
#### Load package dependencies 

library(targets)

#"%!in%" <- Negate("%in%")

tar_option_set(
  packages = c("terra", "data.table", "mgcv", "metagam", "readxl", "stringr", "ggplot2", "vctrs", "cowplot")
)

# PARAMETERS FILE
# set year to run - either "genYr" (smooths over last 5 years) or a specific year (e.g. 2020)
# i_a <- as.numeric(commandArgs(trailingOnly = TRUE)[1]) # array number
# i_a <- 2
# yr <- 2019
# pollutant <- c("co","nox")[i_a] # either "genPoll" (smooths over all 7 EMEP pollutants) or specific pollutant
# vector to be run in dynamic targets (?) or over cluster (?)

# the above loads these packages at a global level for all targets. You can also choose to load them separately 

#####################
#### Load the functions to be used in targets

# Functions used by the targets.
source(here::here("R", "TEMREG_functions.R"))
#source(here::here("R", "xxxxxxxxxxxx.R")) # maybe further analyses with Sentinel etc

#####################
#### Define pipeline

# List of target objects.
list(
    
  ####################
  ## Look-up tables ##
  tar_target(lookup_sectors,  "data/lookups/TEMREG_sectors.csv",   format = "file"),
  tar_target(data_sources,    "data/lookups/profile_id_names.csv", format = "file"),
  #tar_target(lookup_CRF,  "data/lookups/CRF_to_SNAP.csv", format = "file"),
  tar_target(lookup_PID,  "./../../inventory_processor/data/lookups/pollutants.xlsx", format = "file"),
      
  tar_target(dt_NFR_PRO_AGG,  sector_to_ProID(lookup_sectors, data_sources)),
  #tar_target(dt_PRO_AGG,  ProID_to_Agg(lookup_sectors, agg_sectors = c("SNAP","GNFR"))),
  #tar_target(dt_CRF,  fread(lookup_CRF)),
  #tar_target(dt_ISO,  fread(lookup_ISO)),
  tar_target(dt_PID_NAEI,  as.data.table(read_excel(lookup_PID, sheet = "NAEI_pollutants"))),
  
  #tar_target(dt_SNAPGNFR,   as.data.table(read_excel(lookup_sectors, sheet = "SNAPtoGNFR"))),
  #tar_target(dt_GNFRSNAP,   as.data.table(read_excel(lookup_sectors, sheet = "GNFRtoSNAP"))),
  #tar_target(dt_SNAPnames,  as.data.table(read_excel(lookup_sectors, sheet = "SNAPnames"))),
    
  #################################
  ## vectors for dynamic targets ##
  # create a target of aggregation types (required for dynamic targets later)
  # https://books.ropensci.org/targets/dynamic.html
  tar_target(v_yr, c(2018:2021, "genYr")),
  tar_target(v_ProfileIDs, vector_sources(data_sources, level = "Profile_ID")), # this allows data to be read in again if the file(s) change
  tar_target(v_sources, vector_sources(data_sources, level = "data_source")), # this allows data to be read in again if the file(s) change
  tar_target(v_pollutant,  c("nox","nh3","co", "pm25", "sox", "pmco", "nmvoc", "genPoll") ),
  tar_target(v_aggregations, c("SNAP", "GNFR")),
    
  ###########################################
  ## file targets for formatted input data ##
  # making separate targets for separate time scales. 
  # ensures re-read if they are recreated at any point (using v_sources above)
  # also ensures that an update to a month file (e.g.) wont re-run GAM generation for all time scales (down the line). 
  tar_target(v_fname_sourceFiles_hour,     paste0("../TEMREG_data/Data_formatted/hour/",v_sources,"_input_hour.csv"),
                                                  format = "file"),
  tar_target(v_fname_sourceFiles_hourwday, paste0("../TEMREG_data/Data_formatted/hourwday/",v_sources,"_input_hourwday.csv"),
                                                  format = "file"),
  tar_target(v_fname_sourceFiles_wday,     paste0("../TEMREG_data/Data_formatted/wday/",v_sources,"_input_wday.csv"),
                                                  format = "file"),
  tar_target(v_fname_sourceFiles_month,    paste0("../TEMREG_data/Data_formatted/month/",v_sources,"_input_month.csv"),
                                                  format = "file"),
  tar_target(v_fname_sourceFiles_yday,     paste0("../TEMREG_data/Data_formatted/yday/",v_sources,"_input_yday.csv"),
                                                  format = "file"),  
  											  
  ##################################
  ## create Profile ID level GAMs ##
  # once again keep the time scales as separate targets, using sourceFile targets to read all data
  # this means that if a source file is not updated (for that time scale), the GAMs are not regenerated
  # These are activity data GAMs, at the Profile_ID level, and therefore do NOT require a pollutant
  # returns a vector of ProfileID GAM filenames (length = n. Profile IDs), which can be further aggregated into SNAP/GNFR etc 
  tar_target(v_fname_ProIDgams_hour,     GAM_by_ProfileID(v_IDs = v_ProfileIDs, v_fname = v_fname_sourceFiles_hour,
                                                          y = v_yr, time_scale = "hour"), pattern = map(v_yr), format = "file"),
  tar_target(v_fname_ProIDgams_hourwday, GAM_by_ProfileID(v_IDs = v_ProfileIDs, v_fname = v_fname_sourceFiles_hourwday, 
                                                          y = v_yr, time_scale = "hourwday"), pattern = map(v_yr), format = "file"),
  tar_target(v_fname_ProIDgams_wday,     GAM_by_ProfileID(v_IDs = v_ProfileIDs, v_fname = v_fname_sourceFiles_wday,
                                                          y = v_yr, time_scale = "wday"), pattern = map(v_yr), format = "file"),
  tar_target(v_fname_ProIDgams_month,    GAM_by_ProfileID(v_IDs = v_ProfileIDs, v_fname = v_fname_sourceFiles_month,
                                                          y = v_yr, time_scale = "month"), pattern = map(v_yr), format = "file"),
  tar_target(v_fname_ProIDgams_yday,     GAM_by_ProfileID(v_IDs = v_ProfileIDs, v_fname = v_fname_sourceFiles_yday,
                                                          y = v_yr, time_scale = "yday"), pattern = map(v_yr), format = "file"),
  # separate target for year-specific traffic gams from HE data
  tar_target(v_fname_ProIDgams_HE,       list.files("../TEMREG_data/Data/RoadTransport/HE_API/HE_ProID_GAMs", pattern = ".rds$",
                                                    full.names=T), format = "file"), 
    
  ## plot Profile ID level GAMs ##
  # plot all the Profile ID gams - use all Profile ID gams and group by Profile ID into one plot. 
  tar_target(v_fname_ProIDgams_plots, GAM_by_Profile_plot(v_IDs = v_ProfileIDs, y = v_yr, v_fname = vec_c(v_fname_ProIDgams_hour,
																							             v_fname_ProIDgams_hourwday,
																							             v_fname_ProIDgams_wday,
																							             v_fname_ProIDgams_month,
																							             v_fname_ProIDgams_yday)),
                                                          pattern = map(v_yr), format = "file"),
  #################################
  ## Pollutant-based sector GAMs ##  
  # Weighted emissions contribution, via NFR --> AggSector, applied to Profile_IDs
  tar_target(l_agg_GAMs_hour,     GAM_by_Agg_sector(v_fname = v_fname_ProIDgams_hour, y = v_yr, 
                                                    pollutant = v_pollutant, time_scale = "hour", classification = v_aggregations, 
												    dt_lookup = dt_NFR_PRO_AGG, dt_PID = dt_PID_NAEI, inv_year = 2022), 
												    pattern = cross(v_aggregations, v_pollutant, v_yr), format = "file"),
  tar_target(l_agg_GAMs_hourwday, GAM_by_Agg_sector(v_fname = v_fname_ProIDgams_hourwday, y = v_yr,
                                                    pollutant = v_pollutant, time_scale = "hourwday", classification = v_aggregations, 
											     	dt_lookup = dt_NFR_PRO_AGG, dt_PID = dt_PID_NAEI, inv_year = 2022), 
												    pattern = cross(v_aggregations, v_pollutant, v_yr), format = "file"),
  tar_target(l_agg_GAMs_wday,     GAM_by_Agg_sector(v_fname = v_fname_ProIDgams_wday, y = v_yr,
                                                    pollutant = v_pollutant, time_scale = "wday", classification = v_aggregations, 
												    dt_lookup = dt_NFR_PRO_AGG, dt_PID = dt_PID_NAEI, inv_year = 2022), 
												    pattern = cross(v_aggregations, v_pollutant, v_yr), format = "file"),
  tar_target(l_agg_GAMs_month,    GAM_by_Agg_sector(v_fname = v_fname_ProIDgams_month, v_fname_HE = v_fname_ProIDgams_HE, y = v_yr,
                                                    pollutant = v_pollutant, time_scale = "month", classification = v_aggregations, 
												    dt_lookup = dt_NFR_PRO_AGG, dt_PID = dt_PID_NAEI, inv_year = 2022), 
												    pattern = cross(v_aggregations, v_pollutant, v_yr), format = "file"),
  tar_target(l_agg_GAMs_yday,     GAM_by_Agg_sector(v_fname = v_fname_ProIDgams_yday, v_fname_HE = v_fname_ProIDgams_HE, y = v_yr,
                                                    pollutant = v_pollutant, time_scale = "yday", classification = v_aggregations, 
											    	dt_lookup = dt_NFR_PRO_AGG, dt_PID = dt_PID_NAEI, inv_year = 2022), 
											    	pattern = cross(v_aggregations, v_pollutant, v_yr), format = "file"),												
												
  ## plot Agged GAMs ##
  # plot all the aggregated sector gams into composite plots - one plot for all sectors in an aggregation scheme. 
  tar_target(v_fname_AggGAM_plots, GAM_agg_plot(y = v_yr, pollutant = v_pollutant,
                                                classification = v_aggregations, dt_lookup = dt_NFR_PRO_AGG,
												dt_PID = dt_PID_NAEI, inv_year = 2022,
                                                v_fname = vec_c(l_agg_GAMs_hour,
												                l_agg_GAMs_hourwday,
												                l_agg_GAMs_wday,
																l_agg_GAMs_month,
																l_agg_GAMs_yday)),
												pattern = cross(v_aggregations, v_pollutant, v_yr), format = "file"),
   
  #######################
  ## Model input files ##
  # production of model ready data files; target for each model  
  tar_target(v_emep_pollutant,  c("nox","nh3","co", "pm25", "sox", "pmco", "voc") ),
  
  tar_target(v_fname_EMEPinput_hourwday, EMEP4UK_profiles(y = v_yr, pollutant = "genPoll", v_iso = 27, 
                                                          classification = "SNAP", time_scale = "hourwday",
														  v_fname = vec_c(l_agg_GAMs_hour,
												                l_agg_GAMs_hourwday,
												                l_agg_GAMs_wday,
																l_agg_GAMs_month,
																l_agg_GAMs_yday)),
														  pattern = map(v_yr), format = "file"),
  tar_target(v_fname_EMEPinput_wday,     EMEP4UK_profiles(y = v_yr, pollutant = v_emep_pollutant, v_iso = 27, 
                                                          classification = "SNAP", time_scale = "wday",
														  v_fname = vec_c(l_agg_GAMs_hour,
												                l_agg_GAMs_hourwday,
												                l_agg_GAMs_wday,
																l_agg_GAMs_month,
																l_agg_GAMs_yday)),
														  pattern = cross(v_emep_pollutant, v_yr), format = "file"),
  tar_target(v_fname_EMEPinput_month,    EMEP4UK_profiles(y = v_yr, pollutant = v_emep_pollutant, v_iso = 27, 
                                                          classification = "SNAP", time_scale = "month",
														  v_fname = vec_c(l_agg_GAMs_hour,
												                l_agg_GAMs_hourwday,
												                l_agg_GAMs_wday,
																l_agg_GAMs_month,
																l_agg_GAMs_yday)),
														  pattern = cross(v_emep_pollutant, v_yr), format = "file"),
  
  # plot the model input
  tar_target(v_fname_EMEPinput_hourwday_plot, EMEP4UK_profiles_plot(y = v_yr, classification = "SNAP", 
                                                                    time_scale = "hourwday",
																	v_fname = v_fname_EMEPinput_hourwday),
														            pattern = map(v_yr) ),
  tar_target(v_fname_EMEPinput_wday_plot,     EMEP4UK_profiles_plot(y = v_yr, classification = "SNAP", 
                                                                    time_scale = "wday",
																	v_fname = v_fname_EMEPinput_wday),
														            pattern = map(v_yr) ),
  tar_target(v_fname_EMEPinput_month_plot,    EMEP4UK_profiles_plot(y = v_yr, classification = "SNAP", 
                                                                    time_scale = "month",
																	v_fname = v_fname_EMEPinput_month),
														            pattern = map(v_yr) )
  
  
  
)