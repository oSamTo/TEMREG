## initialise the workspace for profiling annual emissions into sector profiles ##

## packages and workspace
packs <- c("metagam","ungeviz","mgcv","stringr","grid","plyr","ggplot2","data.table","stats","lubridate","units","sf", "cowplot","foreach","doParallel","terra")
lapply(packs, require, character.only = TRUE)

"%!in%" <- Negate("%in%")

source("./R/plotting.R")
#spatial
#BNG <<- suppressWarnings(CRS("+init=epsg:27700"))
#LL <<- suppressWarnings(CRS("+init=epsg:4326"))
r_uk_BNG <<- rast(xmin=-230000, xmax = 750000, ymin = -50000, ymax = 1300000, res=1000, crs="EPSG:27700", vals=NA)

# lookup NFR to Profile ID - only one file for this to link to sector lookups
NFR_to_Profile <<- fread("./data/lookup/NFR19_to_GAMsectors.csv")
dt_pollutants <<- fread("./data/lookup/pollutants.csv")


##################################################################################################
#### function to create single Profile ID GAMs that build up into the sectoral GAMs           ####
#### Option to run across timescales and to EXCLUDE Profile_IDs in a vector                   ####
#### Can only run domestic profile GAMs on a modelling PC                                     ####

GAMbyProfile <- function(timescale, exclude = NULL){
  
  print(paste0(Sys.time(), ": Processing ",timescale,"..."))
  
  # read data
  dt <- fread(paste0("./data/input_formatted/",timescale,"_formatted_input.csv"))
  
  v_IDs <- unique(dt[,Profile_ID])
  
  # loop through profile IDs and create a GAM for each one. 
  for(i in v_IDs){
    
    print(paste0(Sys.time(), ":         ",i))
    
    # subset the formatted data
    dt_GAM_data <- dt[Profile_ID == i]
    
    if(timescale != "hourwday") setnames(dt_GAM_data, paste0(timescale), "time")
    
    # re-scale the data around the mean
    if(timescale != "hourwday"){
      dt_GAM_data[, N := N/mean(N,na.rm=T), by=.(Profile_ID)]
    }else{
      dt_GAM_data[, N := N/mean(N,na.rm=T), by=.(Profile_ID, wday)]
    }
    
    dt_GAM_data[,Profile_ID := as.character(Profile_ID)]
    if(timescale == "hourwday") dt_GAM_data[,wday := factor(wday, levels=1:7)]
    
    if(length(unique(dt_GAM_data$N))==1){
      dt_GAM_data[,N := sample(seq(0.99,1.01,0.01), nrow(dt_GAM_data), replace=T)]
    }
    
    # make GAM - some specifics per timescale type.
    
    if(timescale == "hour"){
      gam_sect <- gam(N ~ s(time, bs="cc"), data = dt_GAM_data, method="REML")
    }else if(timescale=="hourwday"){
      gam_sect <- gam(N ~ s(hour, bs="cc", by=wday) + wday, data=dt_GAM_data, method="REML")
    }else if(timescale=="wday"){
      gam_sect <- gam(N ~ s(time, bs="cr", k=6), data = dt_GAM_data, method="REML")
    }else if(timescale=="month"){
      if(i=="RAIL_GEN"){
        gam_sect <- gam(N ~ s(time, bs="cr",k=3), data = dt_GAM_data, method="REML")
      }else{
        gam_sect <- gam(N ~ s(time, bs="cr"), data = dt_GAM_data, method="REML")
      }
      
    }else{
      if(i=="RAIL_GEN"){
        gam_sect <- gam(N ~ s(time, bs="cc",k=3), data = dt_GAM_data, method="REML")
      }else{
        gam_sect <- gam(N ~ s(time, bs="cc"), data = dt_GAM_data, method="REML")
      }
      
    }
    
    #### strip data and save ####
    gam_sr <- strip_rawdata(gam_sect)
    saveRDS(gam_sr, paste0("N:/dump/ukem/",timescale,"/",timescale,"_",i,"_GAM.rds"))
    
    #l_gams[[paste0(timescale,"_",i)]] <- gam_sr
    
    
  } # end of profile ID loop
  
  #saveRDS(l_gams, paste0("N:/dump/ukem/",timescale,"_GAMs_ProfileID.rds"))
  
} # END of function


##################################################################################################
#### function to 1. return the NAEI emissions data, formatted, for the given year & Species.  ####
####             2. match the annual NAEI data to the temporal profile codes and GAm files.   ####

JoinNAEItoProfiles <- function(v_year = NA, species = NA, classification){
  
  ########################################################
  
  # year           = *numeric* year to process. Will determine calendar structure plus year specific profiles.
  #if(!is.numeric(v_year)) stop ("Year is not numeric")
  # species        = *character* name of air pollutant or GHG or metal etc. Needs to conform to a list of options.
  if(species %!in% c("NOx","SOx","CH4","CO2","N2O","NH3", "CO", "NMVOC", "PM25", "PM10", "PMCO","HCL") ) stop ("Species must be in: 
                                            AP:    BaP, CO, NH3, NMVOC, NOx, SO2
                                            PM:    PM25, PM10
                                            GHG:   CH4, CO2, N2O
                                            Metal: Cd, Cu, Hg, Ni, Pb, Zn")
  
  #########################################################
  
  ## lookup tables ##
  # The classification name must have an 'NFR_to_xxxx.csv' file ready for the matching to work. Can be custom. 
  dt_NFR_to_sect <- fread(paste0("./data/lookup/NFR19_to_",classification,".csv")) # NFRs, Profile_IDs & sector groupings (e.g. SNAP)
  dt_NFR_to_sect[NFR19 == "4.00E+01", NFR19 := "4E1"] # how to fix this?? options(scipen=999) not working
  dt_NFR_to_sect[NFR19 == "4.00E+02", NFR19 := "4E2"] # how to fix this?? options(scipen=999) not working
  
  # latest available NAEI data files
  naei_files <- list.files("./data/NAEI_total", pattern="-2019.csv", full.names = T)
  naei_files <- naei_files[grep(paste(dt_pollutants[, file_name], collapse = "|"), naei_files)]
  l_naei <- lapply(naei_files, fread, header=T)
  
  # some formatting (remove blank rows and superfluous columns)
  colskeep <- c("Gas","NFR/CRF Group","Source","Activity",2010:2019)
  l_naei <- lapply(l_naei, function(x) x[Source != ""])
  l_naei <- lapply(l_naei, function(x) x[ ,..colskeep])
  dt_naei <- rbindlist(l_naei, use.names=T)
  suppressMessages(dt_naei[, Gas := plyr::mapvalues(Gas, c(dt_pollutants[,naei_longname]), c(dt_pollutants[,file_name]))])
  
  # restructure data
  dt_naei <- melt(dt_naei, id.vars = c("Gas","NFR/CRF Group","Source","Activity"), variable.name = "year", value.name = "emission")
  setnames(dt_naei, c("NFR/CRF Group"), c("NFR19"))
  suppressWarnings(dt_naei[,emission := as.numeric(emission)])
  dt_naei <- dt_naei[!is.na(emission)]
  dt_naei[Gas == "pmco" & emission <0, emission := 0]
  #if(species == "CO2") dt_naei[, emission := emission/12*44]
  
  # subset according to choice
  if(!(is.na(v_year))) dt_naei <- dt_naei[year %in% c(v_year)]
  if(!(is.na(species))) dt_naei <- dt_naei[Gas %in% dt_pollutants[upper_name == species, file_name]]
  
  # join the NAEI emissions to the classification lookup and to the GAM name lookup
  dt_joined <- dt_NFR_to_sect[dt_naei, on = c("NFR19","Source","Activity")][NFR_to_Profile , on = c("NFR19","Source","Activity")]
  dt_joined <- dt_joined[!is.na(emission)]
  setnames(dt_joined, classification, "sector")
  
  # remove avi.cruise in SNAP classes using anti-join
  if(classification=="SNAP"){
    avi_NFRs <- dt_joined[sector == "avi.cruise"]
    dt_naei <- dt_naei[!avi_NFRs, on=.(NFR19, Source, Activity)]
    dt_joined <- dt_joined[sector != "avi.cruise"]
  } 
  
  # reformat and sum by NFR19/Profile.
  dt_joined <- dt_joined[,c("Gas","NFR19","sector","Profile_ID","year","emission")]
  
  # Calculate the Profile_ID fraction of sector, and take a mean over all pollutants and years
  dt_joined_agg <- dt_joined[, .(emission = sum(emission, na.rm=T)), by=.(Gas, sector, Profile_ID, year)]
  dt_joined_agg[, emis_sec_frac := emission/sum(emission), by=.(Gas, sector, year)]
  dt_joined_agg <- dt_joined_agg[!(is.nan(emis_sec_frac))]
  
  # take the mean contribution by each Profile to each sector. Due to variability in years and pollutants, re-adjust to 1. 
  dt_mean_frac <- dt_joined_agg[, .(emis_sec_frac = mean(emis_sec_frac)), by = .(sector, Profile_ID)]
  dt_mean_frac[, emis_sec_frac := emis_sec_frac/sum(emis_sec_frac), by=sector]
  
  if(classification=="SNAP") dt_mean_frac[, sector := as.numeric(sector)]
  
  setkeyv(dt_mean_frac, c("sector","Profile_ID"))
  
  return(dt_mean_frac)
  
}


#########################################################################################################
#### function to create the sector GAMs as directly above, but in a loop and saving 11 GAM object as list

GAMBYSectorLOOP <- function(year,species, timescale, classification, yr_spec_NFR = NULL){
  
  print(paste0(Sys.time(),": Generating ",classification," sector GAMs for ",timescale,"..."))
  
  ## multiplier to create coefficients centred on 1
  c_mult <- ifelse(timescale == "yday", 365, ifelse(timescale == "month", 12, ifelse(timescale == "wday", 7, ifelse(timescale == "hour", 24, 24))))
  
  ## fetch NAEI data and match to sectors etc. 
  emis <- JoinNAEItoProfiles(v_year = year, species = species, classification)
  
  # change this to read from the emissions/sector lookup - call the profileNAEI function. 
  v_sectors <- emis[,unique(sector)]
  
  l_dt_csvs <- list()
  l_GAM_sectors <- list()
  
  for(s in v_sectors){
    
    print(paste0(Sys.time(),":                ",classification," ",s))
    
    # for all the relevant GAM files:
    ## use the same lapply theory but over profile GAMs to build a table, which you then GAM again
    v_IDs <- emis[sector==s,unique(Profile_ID)]
    
    l_ID_data <- list()
    
    # loop through sectors set above and create a GAM from weighted sub-sector GAMs
    for(i in v_IDs){
      
      i_gam <- readRDS(paste0("./data/GAM_profileID/",timescale,"/",timescale,"_",i,"_GAM.rds"))
      
      # create a table to fit to the model - will need all the Profile_IDs in the GAM here, even if not present in emissions data
      # hourwday needs a separate one making
      
      dt_fit <- blankResponse(v_sectors = i, sectorColName = "Profile_ID", timeColName = "time", timescale)
      
      # fit the GAM to the blank data many times to get a sampling space. 
      # subset the Profile_IDs actually present in the emissions data (possibly different due to pollutant)
      gam_sample <- sampleGAM(s = s, gam = i_gam, dt = dt_fit, n = 50, timescale)
      
      # subset the weighting by Profile_ID and add (for sector GAM)
      prof_weight <- emis[sector==s & Profile_ID == i, emis_sec_frac]
      gam_sample[, w := prof_weight]
      
      # add to list for whole sector
      l_ID_data[[paste0(s,"_",i)]] <- gam_sample
      
    } # ID loop
    
    # rbind it into one table. 
    dt_GAM_data <- rbindlist(l_ID_data, use.names = T)
    # setnames back to standard posit names. 
    if(timescale != "hourwday") setnames(dt_GAM_data, "time", paste0(timescale))
    
    #if(timescale == "hourwday") dt_GAM_data[,wday := factor(wday, levels=1:7)] # think this should be set in blankResponse
    
    ## Make a GAM for the whole sector, weighted by emissions contribution of NFR codes grouped by Profile_ID
    # separate if clauses for little changes in knot selection, bs term etc
    if(timescale == "hour"){
      gam_sect <- gam(N ~ s(hour, bs="cc"), data = dt_GAM_data, weights = w, method="REML")
    }else if(timescale=="hourwday"){
      gam_sect <- gam(N ~ s(hour, bs="cc", by=wday) + wday, data=dt_GAM_data, weights = w, method="REML")
    }else if(timescale=="wday"){
      gam_sect <- gam(N ~ s(wday, bs="cr", k=6), data = dt_GAM_data, weights = w, method="REML")
    }else if(timescale=="month"){
      gam_sect <- gam(N ~ s(month, bs="cr"), data = dt_GAM_data, weights = w, method="REML")
    }else{
      gam_sect <- gam(N ~ s(yday, bs="cc"), data = dt_GAM_data, weights = w, method="REML")
    }
    
    #### strip data and put into list ####
    gam_sr <- strip_rawdata(gam_sect)
    l_GAM_sectors[[s]] <- gam_sr
    
    ## go on to make coefficients table and add to another list, writing at the end. 
    dt_fit <- blankResponse(v_sectors = s, sectorColName = "sector", timeColName = timescale, timescale = timescale)
    
    # fit sector model to empty data
    dt_coeffs <- SectorCoeffs(gam = gam_sr, dt = dt_fit, timescale = timescale)
    
    l_dt_csvs[[paste0(s,"_csv")]] <- dt_coeffs
    
  } # end of sector loop
  
  dt_csvs <- rbindlist(l_dt_csvs, use.names = T)
  
  ## set filenames and write data
  
  if(((year %in% 2010:2019) & (species %in% dt_pollutants[,upper_name]))){
    filename <- paste0("GAM_",timescale,"_",classification,"_",species,"_",year,"_LIST")
  }else if(is.na(year) & (species %in% dt_pollutants[,upper_name])){
    filename <- paste0("GAM_",timescale,"_",classification,"_",species,"_allYr_LIST")
  }else if((year %in% 2010:2019) & is.na(species)){
    filename <- paste0("GAM_",timescale,"_",classification,"_allGas_",year,"_LIST")
  }else{
    filename <- paste0("GAM_",timescale,"_",classification,"_allGas_allYr_LIST")
  }
  
  if(is.na(species)){
    dir.create(file.path(paste0("./output/GAM_sector/",classification,"/allGas")), showWarnings = F)
    dir.create(file.path(paste0("./output/coeffs_sector/",classification,"/allGas")), showWarnings = F)
    
    saveRDS(l_GAM_sectors, paste0("./output/GAM_sector/",classification,"/allGas","/",filename,".rds"))
    fwrite(dt_csvs, paste0("./output/coeffs_sector/",classification,"/allGas","/",filename,".csv"))
  }else{
    dir.create(file.path(paste0("./output/GAM_sector/",classification,"/",species)), showWarnings = F)
    dir.create(file.path(paste0("./output/coeffs_sector/",classification,"/",species)), showWarnings = F)
    
    saveRDS(l_GAM_sectors, paste0("./output/GAM_sector/",classification,"/",species,"/",filename,".rds"))
    fwrite(dt_csvs, paste0("./output/coeffs_sector/",classification,"/",species,"/",filename,".csv"))
  }
  
  #print(paste0(Sys.time(),": DONE."))
  
  
} # end of function


#########################################################################################################
#### function to make a blank fit table, to apply sectors to

blankResponse <- function(v_sectors, sectorColName, timeColName, timescale){
  
  c_mult <- ifelse(timescale == "yday", 365, ifelse(timescale == "month", 12, ifelse(timescale == "wday", 7, ifelse(timescale == "hour", 24, 24))))
  
  if(timescale == "hour"){
    dt <- data.table( name = rep(v_sectors, each=24), time = rep(0:23,length(v_sectors)))
  }else if(timescale=="hourwday"){
    dt <- data.table( name = rep(v_sectors, each=24*7), wday = rep(rep(1:7,each = 24), length(v_sectors)), hour = rep(rep(0:23,7), length(v_sectors)))
  }else{
    dt <- data.table( name = rep(v_sectors, each=c_mult), time = rep(1:c_mult,length(v_sectors)))
  }
  
  dt[, name := factor( name ) ]
  
  if(timescale!="hourwday") setnames(dt, "time", paste0(timeColName))
  if(timescale=="hourwday") dt[,wday := factor(wday, levels = 1:7)]
  
  setnames(dt, "name", sectorColName )
  
  return(dt)
  
}

#########################################################################################################
#### function to sample a GAMs model n. amount of times to sample the uncertainty - for sector GAM   ####

sampleGAM <- function(s, gam, dt, n, timescale){
  
  # s is the sector, for completion
  # gam is the gam to be sampled repeatedly
  # dt provides the response variable to be predicted (i.e. time)
  # n is number of samples to be taken 
  
  i_gam_sampled <- suppressWarnings(as.data.table(ungeviz::sample_outcomes(gam, dt, times = n)))
  
  # centre the samples around 1, as the subsequent GAM is massively effected by different scales.
  if(timescale != "hourwday"){
    i_gam_sampled[, N := N/mean(N, na.rm=T), by=.(Profile_ID, .draw)]
    i_gam_sampled[,c(".draw") := NULL]
  }else{
    i_gam_sampled[, N := N/mean(N, na.rm=T), by=.(Profile_ID, wday, .draw)]
    i_gam_sampled[,c(".draw") := NULL]
  }
  
  # add sector & weights
  i_gam_sampled[, sector := s]
  
  return(i_gam_sampled)
  
}


#########################################################################################################
#### function to fit a sector GAM to the blank table

SectorCoeffs <- function(gam, dt, timescale){
  
  # gam is the gam to build a response table from 
  # dt is the blank response table to fit
  # timescale is the timescale being worked on
  
  # fit sector model to empty data
  fits = predict(gam, newdata=dt, type='response', se=T)
  predicts = as.data.table(data.frame(dt, fits) %>% mutate(lower = fit - 1.96*se.fit, upper = fit + 1.96*se.fit))
  
  if(timescale=="hourwday"){
    predicts[, coeff   := fit/mean(fit, na.rm=T), by=.(sector,wday) ]
    predicts[, c("coeff_l","coeff_u") := list((lower/fit) * coeff , (upper/fit) * coeff), by=.(sector,wday)]
  }else{
    predicts[, coeff   := fit/mean(fit, na.rm=T), by=.(sector)]
    predicts[, c("coeff_l","coeff_u") := list((lower/fit) * coeff , (upper/fit) * coeff), by=.(sector)]
  }
  
  #if(timescale != "hourwday") setnames(predicts, "time", timescale)
  predicts[,c("fit","se.fit","lower","upper") := NULL]
  
  if(timescale=="hourwday"){
    setkeyv(predicts, c("sector","hour","wday"))
  }else{
    setkeyv(predicts, c("sector",paste(timescale)))
  }
  
  return(predicts)
  
} # end of function

#########################################################################################################
#### function to create EMEP4UK inputs - either original simple .dat files or fully 3D netcdfs

EMEP4UKprofiles <- function(year, species = c("NOx","SOx","CH4","CO2","N2O","NH3", "CO", "NMVOC", "PM25", "PMCO", "PM10","HCL"), classification, filetype = c("data file","data cube"), hourly, comp_plot){
  
  filetype <- match.arg(filetype)
  species <- match.arg(species, several.ok = T)
  
  if(is.na(year)) year <- "AllYr" ; print("Using generic year representation")
  if(is.numeric(year) & year < 2010) stop ("If year is defined, it needs to be 2010 or later")
  
  #########################################################
  
  ### large if clause for the filetype;
  
  if(filetype == "data file"){
    
    ## current understanding: 
    # Month and wday files: for each pollutant & SNAP code & country
    # hourly files: by day for each country & SNAP, pollutant generic
    
    ### Hourly ###
    if(hourly == T){
      # read in the hourwday data  - this is pollutant generic
      dt_hourwday <- fread(paste0("./output/coeffs_sector/",classification,"/allGas/GAM_hourwday_",classification,"_allGas_",year,"_LIST.csv"))[,c("coeff_l","coeff_u") := NULL][,coeff := round(coeff, 3)]
      
      # reshape and write out - it'll get overwritten for every species processed but will be identical
      dt_hw <- dcast(dt_hourwday, wday+sector ~ hour, value.var = "coeff")
      
      header <- readLines("./output/model_inputs/EMEP4UK/pre_DUKEMS/HourlyFacs.INERIS",n=3)
      
      f_hourly <- "./output/model_inputs/EMEP4UK/DUKEMS/HourlyFacs_DUKEMS.INERIS"
      file.create(f_hourly)
      cat(header, file = f_hourly, sep="\n")
      fwrite(dt_hw, f_hourly, append = T)
      
    }else{
      
    }
    
    ### Daily and monthly ###
    # read in the wday and month profiles into a list of lists (for each species)
    l_profiles <- lapply(species, collectWdayMonthProfiles, year = year, classification = classification)
    names(l_profiles) <- species
    
    # write the profiles - this includes reading the old ones and inserting UK data
    lapply(l_profiles, writeWdayMonthProfiles, year = year, classification = classification)
    
    if(comp_plot == T){
      
      lapply(species, plotWdayMonthProfiles, year = year, classification = classification)
      
    }else{
      
    }
    
    
    
    
    
    
  }else{
    
    # nothing here yet 
    print(".nc data cube format not available yet.")
    
  } # if else for file type
  
  
} # end of function


#########################################################################################################
#### functions to 1. return list of csv or .rds versions of the sector-level temporal profiles (only wday/month)
####              2. read in original EMEP structure and write with updated UK values

collectWdayMonthProfiles <- function(year, species, classification){
  
  ## read in wday and month data  - this is pollutant specific
  dt_wday <- fread(paste0("./output/coeffs_sector/",classification,"/",species,"/GAM_wday_",classification,"_",species,"_",year,"_LIST.csv"))[,c("coeff_l","coeff_u") := NULL]#[,species := species]
  dt_month <- fread(paste0("./output/coeffs_sector/",classification,"/",species,"/GAM_month_",classification,"_",species,"_",year,"_LIST.csv"))[,c("coeff_l","coeff_u") := NULL]#[,species := species]
  
  l <- list(dt_wday, dt_month)
  l <- lapply(l, function(x) x[,coeff := round(coeff, 3)])
  
  names(l) <- paste0(species,"_",c("wday","month"))
  
  return(l)
}


writeWdayMonthProfiles <- function(year, classification, profiles){
  
  #timestep <- match.arg(timestep)
  #file_time <- ifelse(timestep == "wday", "DailyFac", "MonthlyFac")
  
  spec_func <- str_split(names(profiles),"_")[[1]][1]
  spec_emep <- tolower(spec_func)
  if(spec_func == "NMVOC") spec_emep <- "voc"
  
  ## bring in existing EMEP4UK temporal data
  dt_wday_exist  <- fread(paste0("./output/model_inputs/EMEP4UK/pre_DUKEMS/DailyFac.",spec_emep))
  dt_month_exist <- fread(paste0("./output/model_inputs/EMEP4UK/pre_DUKEMS/MonthlyFac.",spec_emep))
  
  ## WDAY
  ## subset the data according to the species and the time
  dt <- profiles[[paste0(spec_func,"_wday")]]
  dtw <- dcast(dt, sector ~ wday, value.var = "coeff")
  dtw <- data.table(country = 27, dtw)
  names(dtw) <- names(dt_wday_exist)
  
  ## MONTH
  ## subset the data according to the species and the time
  dt <- profiles[[paste0(spec_func,"_month")]]
  dtm <- dcast(dt, sector ~ month, value.var = "coeff")
  dtm <- data.table(country = 27, dtm)
  names(dtm) <- names(dt_month_exist)
  
  
  ## replace UK numbers - remove from original file, append new and order
  ## WDAY
  dt_wday_exist <- dt_wday_exist[V1 != 27]
  dt_wday_write <- rbindlist(list(dt_wday_exist, dtw), use.names = T)
  dt_wday_write <- dt_wday_write[order(V2, V1)]
  
  ## MONTH
  dt_month_exist <- dt_month_exist[V1 != 27]
  dt_month_write <- rbindlist(list(dt_month_exist, dtm), use.names = T)
  dt_month_write <- dt_month_write[order(V2, V1)]
  
  
  # write out data
  fwrite(dt_wday_write, paste0("./output/model_inputs/EMEP4UK/DUKEMS/DailyFac_DUKEMS.",spec_emep),sep = " ", col.names = F) 
  fwrite(dt_month_write, paste0("./output/model_inputs/EMEP4UK/DUKEMS/MonthlyFac_DUKEMS.",spec_emep),sep = " ", col.names = F) 
  
}


plotWdayMonthProfiles <- function(species, year, classification){
  
  spec_emep <- tolower(species)
  if(species == "NMVOC") spec_emep <- "voc"
 
  
  ##########
  ## HOUR ##
  dt_hour_exist  <- suppressWarnings(fread(paste0("./output/model_inputs/EMEP4UK/pre_DUKEMS/HourlyFacs.INERIS")))
  names(dt_hour_exist) <- c("wday","SNAP",1:24)
  dt_hour_exist_m <- melt(dt_hour_exist, id.vars = c("wday","SNAP"), variable.name = "hour", value.name = "coeff")
  dt_hour_exist_m[, source := "EMPE4UK"]
  
  
  dt_hour_duke <- suppressWarnings(fread(paste0("./output/model_inputs/EMEP4UK/DUKEMS/HourlyFacs_DUKEMS.INERIS")))
  names(dt_hour_duke) <- c("wday","SNAP",1:24)
  dt_hour_duke_m <- melt(dt_hour_duke, id.vars = c("wday","SNAP"), variable.name = "hour", value.name = "coeff")
  dt_hour_duke_m[, source := "DUKEMS"]
  
  dt_hour <- rbindlist(list(dt_hour_exist_m, dt_hour_duke_m), use.names = T)
  
  g1 <- ggplot(data = dt_hour[SNAP < 11], aes(x = hour, y = coeff, group = source, colour = source))+
    geom_line()+
    facet_grid(wday~SNAP)+
    theme_bw()
  
  ##########
  ## WDAY ##
  dt_wday_exist  <- suppressWarnings(fread(paste0("./output/model_inputs/EMEP4UK/pre_DUKEMS/DailyFac.",spec_emep)))[V1==27]
  names(dt_wday_exist) <- c("country","SNAP",1:7)
  dt_wday_exist_m <- melt(dt_wday_exist, id.vars = c("country","SNAP"), variable.name = "wday", value.name = "coeff")
  dt_wday_exist_m[, source := "EMPE4UK"]
  
  
  dt_wday_duke <- suppressWarnings(fread(paste0("./output/model_inputs/EMEP4UK/DUKEMS/DailyFac_DUKEMS.",spec_emep)))[V1==27]
  names(dt_wday_duke) <- c("country","SNAP",1:7)
  dt_wday_duke_m <- melt(dt_wday_duke, id.vars = c("country","SNAP"), variable.name = "wday", value.name = "coeff")
  dt_wday_duke_m[, source := "DUKEMS"]
  
  dt_wday <- rbindlist(list(dt_wday_exist_m, dt_wday_duke_m), use.names = T)
  
  g2 <- ggplot(data = dt_wday[SNAP < 11], aes(x = wday, y = coeff, group = source, colour = source))+
    geom_line()+
    facet_wrap(~SNAP, nrow=2)+
    theme_bw()
  
  ###########
  ## MONTH ##
  dt_month_exist  <- suppressWarnings(fread(paste0("./output/model_inputs/EMEP4UK/pre_DUKEMS/MonthlyFac.",spec_emep)))[V1==27]
  names(dt_month_exist) <- c("country","SNAP",1:12)
  dt_month_exist_m <- melt(dt_month_exist, id.vars = c("country","SNAP"), variable.name = "month", value.name = "coeff")
  dt_month_exist_m[, source := "EMPE4UK"]
  
  
  dt_month_duke <- suppressWarnings(fread(paste0("./output/model_inputs/EMEP4UK/DUKEMS/MonthlyFac_DUKEMS.",spec_emep)))[V1==27]
  names(dt_month_duke) <- c("country","SNAP",1:12)
  dt_month_duke_m <- melt(dt_month_duke, id.vars = c("country","SNAP"), variable.name = "month", value.name = "coeff")
  dt_month_duke_m[, source := "DUKEMS"]
  
  dt_month <- rbindlist(list(dt_month_exist_m, dt_month_duke_m), use.names = T)
  
  g3 <- ggplot(data = dt_month[SNAP < 11], aes(x = month, y = coeff, group = source, colour = source))+
    geom_line()+
    facet_wrap(~SNAP, nrow=2)+
    theme_bw()
  
  
  l_p <- list(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3))
  
  # write out data
  c1 <- plot_grid(l_p[[1]], l_p[[2]], l_p[[3]], ncol=1, nrow=3,rel_widths = 1, rel_heights = c(1.2,0.8,0.8), align="h", axis="l")
  
  save_plot(paste0("./output/model_inputs/EMEP4UK/comparison_plots/",species,"_EMPE4UK_DUKEMS_profiles.png"), c1, base_height = 15, base_width = 9)
  
  
}


#######################################################################################################
## NOT IN USE ##
#########################################################################################################
#### function to sample the GAM model space for every relevant profile_ID in a sector
#### also assigns the weighting, which is the Profile_ID emissions share of parent sector
#### returns a list of data tables, one for each sector passed

#sampledProfileGAMs <- function(v_sectors, emis, timescale){
  
#  c_mult <- ifelse(timescale == "yday", 365, ifelse(timescale == "month", 12, ifelse(timescale == "wday", 7, ifelse(timescale == "hour", 24, 24))))
  
  # get unqiue sector names from the classification system
  #setnames(emis, classification, "sector")
#  v_IDs <- emis[sector==v_sectors,unique(Profile_ID)]
  
#  l_ID_data <- list()
  
  # loop through sectors set above and create a GAM from weighted sub-sector GAMs
#  for(i in v_IDs){
    
    # subset the weighting by Profile_ID
#    prof_weight <- emis[sector==v_sectors & Profile_ID == i, emis_sec_frac]
    
#    i_gam <- readRDS(paste0("./data/GAM_profileID/",timescale,"/",timescale,"_",i,"_GAM.rds"))
    
    # create a table to fit to the model - will need all the Profile_IDs in the GAM here, even if not present in emissions data
    # hourwday needs a separate one making
#    if(timescale == "hourwday"){
#      dt_fit <- data.table(Profile_ID = rep(i, each=24*7), wday = rep(1:7,each=24), hour = rep(1:24,7) )
#    }else{
#      dt_fit <- data.table(Profile_ID = i, time = 1:c_mult)
#    }
    
#    dt_fit[, Profile_ID := factor(Profile_ID)] 
#    if(timescale == "hour") dt_fit[, time := time-1]
#    if(timescale == "hourwday") dt_fit[, hour := hour-1]
#    if(timescale == "hourwday") dt_fit[, wday := factor(wday, levels=1:7)]
    
    # fit the GAM to the blank data many times to get a sampling space. 
    # subset the Profile_IDs actually present in the emissions data (possibly different due to pollutant)
#    dt_g_sampled <- suppressWarnings(as.data.table(ungeviz::sample_outcomes(i_gam, dt_fit, times = 50)))
    
    # centre the samples around 1, as the subsequent GAM is massively effected by different scales.
#    if(timescale != "hourwday"){
#      dt_g_sampled[, N := N/mean(N, na.rm=T), by=.(Profile_ID, .draw)]
#      dt_g_sampled[,c(".draw") := NULL]
#    }else{
#      dt_g_sampled[, N := N/mean(N, na.rm=T), by=.(Profile_ID, wday, .draw)]
#      dt_g_sampled[,c(".draw") := NULL]
#    }
    
    # add sector & weights
#    dt_g_sampled[, sector := v_sectors]
#    dt_g_sampled[, w := prof_weight]
    
    # add to list for whole sector
#    l_ID_data[[paste0(v_sectors,"_",i)]] <- dt_g_sampled
    
#  } # ID loop
  
#  dt_ID_data <- rbindlist(l_ID_data, use.names = T)
#  return(dt_ID_data)
  
  
#} #  end of function

#######################################################################################################
#### Function to create weighted sector level GAMs via classification system, using NFR/Profiles.  ####
#### Result: weighted temporal profiles for nominated Sectors, based on emissions per profile      ####
####       : flat csv outputs of coefficients                                                      #### 

## NOT IN USE - by= argument is far too computational to work, especially for yday and hourwday. 
#GAMBySector <- function(year, species, timescale, classification, yr_spec_NFR = NULL){
  
  ####################################################
  
  # year           = *numeric* year to process OR NA for generic.
  # species        = *character* name of air pollutant or GHG or metal etc OR NA for generic
#  if(species %!in% c("NOx","SOx","CH4","CO2","N2O","NH3") & !is.na(species) ) stop ("POllutant must be in: 
#                                            AP:    CO, NH3, NMVOC, NOx, SO2
#                                            PM:    PM2.5, PM10
#                                            GHG:   CH4, CO2, N2O
#                                            OR pollutant must be NA to allow for generic pollutant type")
  
  # timescale       = *character vector* timesteps to recreate profiles for
#  if(timescale %!in% c("yday","month","wday","hour","hourwday")) stop("Not currently a valid timestep to run profiles for")
  
  # yr_spec_NFR    = *character* optional vector of NFR codes to be profiled by year specific data
  
  ####################################################
  
#  print(paste0(Sys.time(),": Generating ",classification," sector GAMs for ",timescale,"..."))
  
  ## multiplier to create coefficients centred on 1
#  c_mult <- ifelse(timescale == "yday", 365, ifelse(timescale == "month", 12, ifelse(timescale == "wday", 7, ifelse(timescale == "hour", 24, 24))))
  
  ## fetch NAEI data and match to sectors etc. 
#  emis <- JoinNAEItoProfiles(v_year = year, species = species, classification)
  
  ## read in the correct GAM library (by timescale)
  #l_gam_time <- readRDS(paste0("./data/GAM_profileID/",timescale,"_GAMs_ProfileID.rds"))
  
  # change this to read from the emissions/sector lookup - call the profileNAEI function. 
#  v_sectors <- emis[,unique(sector)]
  
#  print(paste0(Sys.time(),":                Collate & sample data..."))
  
#  l_sector_data <- lapply(v_sectors, sampledProfileGAMs, emis = emis, timescale = timescale )
  
  # combine all the data together - this is the new sector GAM data
#  dt_GAM_data <- rbindlist(l_sector_data, use.names = T)
#  dt_GAM_data[, w := as.numeric(as.character(w))]
#  dt_GAM_data[, sector := factor(sector)]
  #if(timescale == "hourwday") dt_GAM_data[,wday := factor(wday, levels=1:7)]
  
#  print(paste0(Sys.time(),":                Run GAM model..."))
  
  ## Make a GAM for the whole sector, weighted by emissions contribution of NFR codes grouped by Profile_ID
  # separate if clauses for little changes in knot selection, bs term etc
#  if(timescale == "hour"){
#    gam_sect <- gam(N ~ s(time, bs="cc", by = sector) + sector, data = dt_GAM_data, weights = w, method="REML")
#  }else if(timescale=="hourwday"){
#    gam_sect <- gam(N ~ s(hour, bs="cc", by = sector) + s(wday, bs="cr", k=6, by = sector) + sector, data=dt_GAM_data, weights = w, method="REML")
#  }else if(timescale=="wday"){
#    gam_sect <- gam(N ~ s(time, bs="cr", k=6, by = sector) + sector, data = dt_GAM_data, weights = w, method="REML")
#  }else if(timescale=="month"){
#    gam_sect <- gam(N ~ s(time, bs="cr", by = sector) + sector, data = dt_GAM_data, weights = w, method="REML")
#  }else{
#    gam_sect <- gam(N ~ s(time, bs="cc", by = sector) + sector, data = dt_GAM_data, weights = w, method="REML")
#  }
#  print(Sys.time())
  
  #### strip data and save ####
#  gam_sr <- strip_rawdata(gam_sect)
  
#  rm(gam_sect)
#  gc()
  
  ## make a coefficients table as well:
  ## create a blank table for GAM estimation
#  print(paste0(Sys.time(),":                Coefficients csv and write data..."))
  
#  dt_gc <- blankResponse(v_sectors = v_sectors, timescale = timescale)
#  dt_gc[, sector := factor(sector)]
  
  # fit sector model to empty data
#  dt_coeffs <- SectorCoeffs(gam_sr = gam_sr, dt_gc = dt_gc, timescale = timescale)
#  if(timescale != "hourwday") setnames(dt_coeffs, "time", timescale)
#  dt_coeffs[,c("fit","se.fit","lower","upper") := NULL]
  
#  if(timescale=="hourwday"){
#    setkeyv(dt_coeffs, c("sector","hour","wday"))
#  }else{
#    setkeyv(dt_coeffs, c("sector",paste(timescale)))
#  }
  
  ## set filenames and write data
  
#  if(((year %in% 2010:2019) & (species %in% dt_pollutants[,upper_name]))){
#    filename <- paste0("GAM_",timescale,"_",classification,"_",species,"_",year)
#  }else if(is.na(year) & (species %in% dt_pollutants[,upper_name])){
#    filename <- paste0("GAM_",timescale,"_",classification,"_",species,"_allYr")
#  }else if((year %in% 2010:2019) & is.na(species)){
#    filename <- paste0("GAM_",timescale,"_",classification,"_allGas_",year)
#  }else{
#    filename <- paste0("GAM_",timescale,"_",classification,"_allGas_allYr")
#  }
  
#  if(is.na(species)){
#    dir.create(file.path(paste0("./output/GAM_sector/",classification,"/allGas")), showWarnings = F)
#    dir.create(file.path(paste0("./output/coeffs_sector/",classification,"/allGas")), showWarnings = F)
    
#    saveRDS(gam_sr, paste0("./output/GAM_sector/",classification,"/allGas","/",filename,".rds"))
#    fwrite(dt_coeffs, paste0("./output/coeffs_sector/",classification,"/allGas","/",filename,".csv"))
#  }else{
#    dir.create(file.path(paste0("./output/GAM_sector/",classification,"/",species)), showWarnings = F)
#    dir.create(file.path(paste0("./output/coeffs_sector/",classification,"/",species)), showWarnings = F)
    
#    saveRDS(gam_sr, paste0("./output/GAM_sector/",classification,"/",species,"/",filename,".rds"))
#    fwrite(dt_coeffs, paste0("./output/coeffs_sector/",classification,"/",species,"/",filename,".csv"))
#  }
  
#  print(paste0(Sys.time(),": DONE."))
  
#} # end of function

