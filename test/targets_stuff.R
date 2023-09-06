library(targets)

tar_option_set(
  packages = c("data.table")
)

#i_a <- as.numeric(commandArgs(trailingOnly = TRUE)[1]) # array number
#v_thing <- c("sam","jim","tom")[i_a]

source("test/stuff.R")


#####################
#### Define pipeline

# List of target objects.
list(
    
  ####################
  ## Look-up tables ##
  tar_target(dir_targ,  create_DIR(name = v_thing),   format = "file"),
  tar_target(data_targ, make_output(name = v_thing, loc = dir_targ, n_sector = 11, i_a = i_a), format = "file")
  
  )