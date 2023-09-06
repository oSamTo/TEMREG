
library(targets)
#library(options)

define_option(
  "thing",
  #default = "shout",
  #desc = "Print output in uppercase ('shout') or lowercase ('whisper')",
  option_name = "thing",
  envvar_name = "THING"
)

define_option(
  "array",
  #default = "shout",
  #desc = "Print output in uppercase ('shout') or lowercase ('whisper')",
  option_name = "array",
  envvar_name = "ARRAY"
)

# To establish a different store and script per project, write a top-level _targets.yaml 
# configuration to specify these paths explicitly. You can do this from R with tar_config_set().
i_a <- as.numeric(commandArgs(trailingOnly = TRUE)[1]) # array number

v_thing <- c("sam","jim","tom")[i_a]

print(v_thing)

v_store   <- paste0("test/target_stores/store_",v_thing)
#v_project <- paste0("project_",v_thing)

tar_config_set(script = "test/targets_stuff.R", store = v_store, project = "project_sam")
#tar_config_set(script = "target_flows/targets_EMEP.R", store = "target_stores/store_EMEP", project = "project_EMEP")
#tar_config_set(script = "target_flows/targets_MapEire.R", store = "target_stores/store_MapEire", project = "project_MapEire")

# serial
Sys.setenv(TAR_PROJECT = "project_sam", THING = v_thing, ARRAY = i_a)

tar_make()
#tar_make(envir = "envir1")


tar_meta(fields = warnings, complete_only = T)


# parallel
#system.time(tar_make_future(workers = 6L))

#source("test/stuff.R")


