library(targets)
#####################
#### Inspect the pipeline
#####################

#tar_manifest(fields = all_of("command"))
#tar_glimpse()
#tar_visnetwork()
#tar_outdated() # what is out of date?

# To establish a different store and script per project, write a top-level _targets.yaml 
# configuration to specify these paths explicitly. You can do this from R with tar_config_set().

#tar_config_set(script = "target_flows/targets_NAEI.R", store = "target_stores/store_NAEI", project = "project_NAEI")

yr <- "genYr"
# serial

#Sys.setenv(TAR_PROJECT = "project_NAEI")

yr_text <- ifelse(yr == "genYr", "non-specific year", yr)

print(paste0(Sys.time(),": Running TEMREG for ", yr_text))

tar_make(envir = "envir1")
tar_meta(fields = warnings, complete_only = T)
tar_meta(fields = error, complete_only = TRUE)

# parallelise for pollutants
#system.time(tar_make_future(workers = 6L))

