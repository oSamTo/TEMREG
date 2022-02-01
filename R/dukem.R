## initialise the workspace for profiling annual emissions into sector profiles ##

## packages and workspace
packs <- c("metagam","ungeviz","mgcv","sp","raster","stringr","gdalUtils","rgeos","rgdal","grid","plyr","ggplot2","data.table","stats","sf","lubridate","units", "cowplot")
lapply(packs, require, character.only = TRUE)

"%!in%" <- Negate("%in%")

#spatial
BNG <<- suppressWarnings(CRS("+init=epsg:27700"))
LL <<- suppressWarnings(CRS("+init=epsg:4326"))
r_uk_BNG <<- raster(xmn=-230000, xmx = 750000, ymn = -50000, ymx = 1300000, res=1000, crs=BNG, vals=NA)

# lookup NFR to Profile ID - only one file for this to link to sector lookups
NFR_to_Profile <<- fread("./data/lookup/NFR_to_ProfileID.csv")
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
#### function to plot the profile_ID GAMs for comparison and to inspect the shapes of the     ####
#### GAMs going into the sector level processing (built on activity data)                     ####

GAMprofilePlot <- function(){
  
  print(paste0(Sys.time(),": Plotting all Profile_IDs in composite..."))
  
  v_IDs <- unique(NFR_to_Profile[,Profile_ID])
  
  # big loop to produce lots of plots. Dont need to do this but it's helpful.
  
  for(i in v_IDs){
    
    print(paste0(Sys.time(),":                   ", i))
    
    l_plot <- list()
    
    for(timescale in c("hour","hourwday","wday","yday")){
      
      #l_plot_data <- list()
      
      gam_sr <- readRDS(paste0("./data/GAM_profileID/",timescale,"/",timescale,"_",i,"_GAM.rds"))
        
      dt_gc <- blankFitSector(v_sectors = i, timescale = timescale)
      dt_gc[, sector := factor(sector)]
      if(timescale=="hourwday") dt_gc[,wday := factor(wday, levels = 1:7)]
      
      # fit sector model to empty data
      dt_coeffs <- SectorCoeffs(gam_sr = gam_sr, dt_gc = dt_gc, timescale = timescale)
      #if(timescale != "hourwday") setnames(dt_coeffs, "time", timescale)
      dt_coeffs[,c("fit","se.fit","lower","upper") := NULL]
      
      if(timescale=="hourwday"){
        setkeyv(dt_coeffs, c("sector","hour","wday"))
      }else{
        setkeyv(dt_coeffs, c("sector","time"))
      }
      
      
      g1 <- ggplot()+
        {if(timescale == "hourwday")geom_line(data=dt_coeffs, aes(x=hour,y=coeff, group=wday, colour=wday))}+
        {if(timescale != "hourwday")geom_line(data=dt_coeffs, aes(x=time,y=coeff), colour="red")}+
        {if(timescale == "hourwday")geom_ribbon(data=dt_coeffs, aes(x=hour, ymin=coeff_l, ymax=coeff_u, group=wday, fill=wday), alpha=0.3)}+
        {if(timescale != "hourwday")geom_ribbon(data=dt_coeffs, aes(x=time, ymin=coeff_l, ymax=coeff_u), alpha=0.3)}+
        labs(x=timescale,y="Coefficient")+
        #{if(timescale == "hourwday")facet_wrap(~wday, nrow=2)}+
        theme_bw()+
        theme(axis.text.x = element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.title = element_text(size=16),
              legend.title=element_blank())
      
      l_plot[[timescale]] <- ggplotGrob(g1)
      
    } # time loop
    
    # stitch plots together
    
    d1 <- plot_grid(l_plot[["hourwday"]], l_plot[["hour"]], l_plot[["wday"]], l_plot[["yday"]], ncol=1, nrow=4,rel_heights = c(1.5,1,1,1), rel_widths = c(1,1,1,1), align="h", axis="l")
    
    save_plot(paste0("./data/GAM_profileID/plots/GAM_",i,"_plot.png"), d1, base_height = 14, base_width = 14)
    
  } # end of profile ID loop
  
  print(paste0(Sys.time(),": DONE."))
  
} # end of function

##################################################################################################
#### function to 1. return the NAEI emissions data, formatted, for the given year & Species.  ####
####             2. match the annual NAEI data to the temporal profile codes and GAm files.   ####

JoinNAEItoProfiles <- function(v_year = NA, species = NA, classification){
  
  ########################################################
  
  # year           = *numeric* year to process. Will determine calendar structure plus year specific profiles.
  if(!is.numeric(v_year)) stop ("Year is not numeric")
  # species        = *character* name of air pollutant or GHG or metal etc. Needs to conform to a list of options.
  if(species %!in% c("NOx","SOx","CH4","CO2","N2O","NH3") ) stop ("Species must be in: 
                                            AP:    BaP, CO, NH3, NMVOC, NOx, SO2
                                            PM:    PM2.5, PM10
                                            GHG:   CH4, CO2, N2O
                                            Metal: Cd, Cu, Hg, Ni, Pb, Zn")
  
  #########################################################
  
  ## lookup tables ##
  # The classification name must have an 'NFR_to_xxxx.csv' file ready for the matching to work. Can be custom. 
  dt_NFR_to_sect <- fread(paste0("./data/lookup/NFR_to_",classification,".csv")) # NFRs, Profile_IDs & sector groupings (e.g. SNAP)
  dt_NFR_to_sect[NFR19 == "4.00E+01", NFR19 := "4E1"]
  dt_NFR_to_sect[NFR19 == "4.00E+02", NFR19 := "4E2"]
  
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
  if(species == "CO2") dt_naei[, emission := emission/12*44]
  
  # subset according to choice
  if(!(is.na(year))) dt_naei <- dt_naei[year %in% c(v_year)]
  if(!(is.na(species))) dt_naei <- dt_naei[Gas %in% dt_pollutants[upper_name %in% species, file_name]]
  
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
  
  #if(identical(dt_naei[,sum(emission, na.rm=T)], dt_joined[,sum(emission, na.rm=T)])){
  #  NULL
  #}else{
  #  print("Emissions total has changed since joining to profile table - CHECK")
  #  break
  #}
 
  # reformat and sum by NFR19/Profile.
  dt_joined <- dt_joined[,c("Gas","NFR19","sector","Profile_ID","year","emission")]
  
  # Calculate the Profile_ID fraction of sector, and take a mean over all pollutants and years
  dt_joined_agg <- dt_joined[, .(emission = sum(emission, na.rm=T)), by=.(Gas, sector, Profile_ID, year)]
  dt_joined_agg[, emis_sec_frac := emission/sum(emission), by=.(Gas, sector, year)]
  dt_joined_agg <- dt_joined_agg[!(is.nan(emis_sec_frac))]
  
  # take the mean contribution by each Profile to each sector. Due to variability in years and pollutants, re-adjust to 1. 
  dt_mean_frac <- dt_joined_agg[, .(emis_sec_frac = mean(emis_sec_frac)), by = .(sector, Profile_ID)]
  dt_mean_frac[, emis_sec_frac := emis_sec_frac/sum(emis_sec_frac), by=sector]
  
  setkeyv(dt_mean_frac, c("sector","Profile_ID"))
  
  return(dt_mean_frac)
  
}

#########################################################################################################
#### function to sample the GAM model space for every relevant profile_ID in a sector
#### also assigns the weighting, which is the Profile_ID emissions share of parent sector
#### returns a list of data tables, one for each sector passed

sampledProfileGAMs <- function(v_sectors, emis, timescale){
  
  c_mult <- ifelse(timescale == "yday", 365, ifelse(timescale == "month", 12, ifelse(timescale == "wday", 7, ifelse(timescale == "hour", 24, 24))))
  
    # get unqiue sector names from the classification system
    #setnames(emis, classification, "sector")
    v_IDs <- emis[sector==v_sectors,unique(Profile_ID)]
    
    l_ID_data <- list()
    
    # loop through sectors set above and create a GAM from weighted sub-sector GAMs
    for(i in v_IDs){
      
      # subset the weighting by Profile_ID
      prof_weight <- emis[sector==v_sectors & Profile_ID == i, emis_sec_frac]
      
      i_gam <- readRDS(paste0("./data/GAM_profileID/",timescale,"/",timescale,"_",i,"_GAM.rds"))
      
      # create a table to fit to the model - will need all the Profile_IDs in the GAM here, even if not present in emissions data
      # hourwday needs a separate one making
      if(timescale == "hourwday"){
        dt_fit <- data.table(Profile_ID = rep(i, each=24*7), wday = rep(1:7,each=24), hour = rep(1:24,7) )
      }else{
        dt_fit <- data.table(Profile_ID = i, time = 1:c_mult)
      }
      
      dt_fit[, Profile_ID := factor(Profile_ID)] 
      if(timescale == "hour") dt_fit[, time := time-1]
      if(timescale == "hourwday") dt_fit[, hour := hour-1]
      
      # fit the GAM to the blank data many times to get a sampling space. 
      # subset the Profile_IDs actually present in the emissions data (possibly different due to pollutant)
      dt_g_sampled <- suppressWarnings(as.data.table(ungeviz::sample_outcomes(i_gam, dt_fit, times = 50)))
      
      # centre the samples around 1, as the subsequent GAM is massively effected by different scales.
      if(timescale != "hourwday"){
        dt_g_sampled[, N := N/mean(N, na.rm=T), by=.(Profile_ID, .draw)]
        dt_g_sampled[,c(".draw") := NULL]
      }else{
        dt_g_sampled[, N := N/mean(N, na.rm=T), by=.(Profile_ID, wday, .draw)]
        dt_g_sampled[,c(".draw") := NULL]
      }
      
      # add sector & weights
      dt_g_sampled[, sector := v_sectors]
      dt_g_sampled[, w := prof_weight]
      
      # add to list for whole sector
      l_ID_data[[paste0(v_sectors,"_",i)]] <- dt_g_sampled
      
    } # ID loop
    
    dt_ID_data <- rbindlist(l_ID_data, use.names = T)
    return(dt_ID_data)
    
 
} #  end of function


#########################################################################################################
#### function to make a blank fit table, to apply a sectoral level GAM to

blankFitSector <- function(v_sectors, timescale){
  
  c_mult <- ifelse(timescale == "yday", 365, ifelse(timescale == "month", 12, ifelse(timescale == "wday", 7, ifelse(timescale == "hour", 24, 24))))
  
  if(timescale == "hour"){
    dt <- data.table(sector = rep(v_sectors, each=24), time = rep(0:23,length(v_sectors)))
  }else if(timescale=="hourwday"){
    dt <- data.table(sector = rep(v_sectors, each=24*7), wday = rep(rep(1:7,each = 24), length(v_sectors)), hour = rep(rep(0:23,7), length(v_sectors)))
  }else{
    dt <- data.table(sector = rep(v_sectors, each=c_mult), time = rep(1:c_mult,length(v_sectors)))
  }
  
  dt[, sector := factor(sector)]
  
  return(dt)
  
}

#########################################################################################################
#### function to fit a sector GAM to the blank table

SectorCoeffs <- function(gam_sr, dt_gc, timescale){
  
  # fit sector model to empty data
  fits = predict(gam_sr, newdata=dt_gc, type='response', se=T)
  predicts = as.data.table(data.frame(dt_gc, fits) %>% mutate(lower = fit - 1.96*se.fit, upper = fit + 1.96*se.fit))
  
  if(timescale=="hourwday"){
    predicts[, coeff   := fit/mean(fit, na.rm=T), by=.(sector,wday) ]
    predicts[, c("coeff_l","coeff_u") := list((lower/fit) * coeff , (upper/fit) * coeff), by=.(sector,wday)]
  }else{
    predicts[, coeff   := fit/mean(fit, na.rm=T), by=.(sector)]
    predicts[, c("coeff_l","coeff_u") := list((lower/fit) * coeff , (upper/fit) * coeff), by=.(sector)]
  }
  
  return(predicts)
  
} # end of function



#######################################################################################################
#### Function to create weighted sector level GAMs via classification system, using NFR/Profiles.  ####
#### Result: weighted temporal profiles for nominated Sectors, based on emissions per profile      ####
####       : flat csv outputs of coefficients                                                      #### 

GAMBySector <- function(year, species, timescale, classification, yr_spec_NFR = NULL){
  
  ####################################################
  
  # year           = *numeric* year to process OR NA for generic.
  # species        = *character* name of air pollutant or GHG or metal etc OR NA for generic
  if(species %!in% c("NOx","SOx","CH4","CO2","N2O","NH3") & !is.na(species) ) stop ("POllutant must be in: 
                                            AP:    CO, NH3, NMVOC, NOx, SO2
                                            PM:    PM2.5, PM10
                                            GHG:   CH4, CO2, N2O
                                            OR pollutant must be NA to allow for generic pollutant type")
  
  # timescale       = *character vector* timesteps to recreate profiles for
  if(timescale %!in% c("yday","month","wday","hour","hourwday")) stop("Not currently a valid timestep to run profiles for")
  
  # yr_spec_NFR    = *character* optional vector of NFR codes to be profiled by year specific data
  
  ####################################################
  
  print(paste0(Sys.time(),": Generating ",classification," sector GAMs for ",timescale,"..."))
  
  ## multiplier to create coefficients centred on 1
  c_mult <- ifelse(timescale == "yday", 365, ifelse(timescale == "month", 12, ifelse(timescale == "wday", 7, ifelse(timescale == "hour", 24, 24))))
  
  ## fetch NAEI data and match to sectors etc. 
  emis <- JoinNAEItoProfiles(v_year = year, species = species, classification)
  
  ## read in the correct GAM library (by timescale)
  #l_gam_time <- readRDS(paste0("./data/GAM_profileID/",timescale,"_GAMs_ProfileID.rds"))
  
  # change this to read from the emissions/sector lookup - call the profileNAEI function. 
  v_sectors <- emis[,unique(sector)]
  
  print(paste0(Sys.time(),":                Collate & sample data..."))
  
  l_sector_data <- lapply(v_sectors, sampledProfileGAMs, emis = emis, timescale = timescale )
  
  # combine all the data together - this is the new sector GAM data
  dt_GAM_data <- rbindlist(l_sector_data, use.names = T)
  dt_GAM_data[, w := as.numeric(as.character(w))]
  dt_GAM_data[, sector := factor(sector)]
  
  print(paste0(Sys.time(),":                Run GAM model..."))
  
  ## Make a GAM for the whole sector, weighted by emissions contribution of NFR codes grouped by Profile_ID
  # separate if clauses for little changes in knot selection, bs term etc
  if(timescale == "hour"){
    gam_sect <- gam(N ~ s(time, bs="cc", by = sector) + sector, data = dt_GAM_data, weights = w, method="REML")
  }else if(timescale=="hourwday"){
    gam_sect <- gam(N ~ s(hour, bs="cc", by = sector) + s(wday, bs="cr", k=6, by = sector) + sector, data=dt_GAM_data, weights = w, method="REML")
  }else if(timescale=="wday"){
    gam_sect <- gam(N ~ s(time, bs="cr", k=6, by = sector) + sector, data = dt_GAM_data, weights = w, method="REML")
  }else if(timescale=="month"){
    gam_sect <- gam(N ~ s(time, bs="cr", by = sector) + sector, data = dt_GAM_data, weights = w, method="REML")
  }else{
    gam_sect <- gam(N ~ s(time, bs="cc", by = sector) + sector, data = dt_GAM_data, weights = w, method="REML")
  }
  print(Sys.time())
  
  #### strip data and save ####
  gam_sr <- strip_rawdata(gam_sect)
  
  rm(gam_sect)
  gc()
  
  #######################################################################################
  #######################################################################################
  
  ### BACKUP CODE FOR LIST OFF 11 SECTORS
  # loop through sectors set above and create a GAM from weighted sub-sector GAMs
  for(s in v_sectors){
    
    print(paste0(Sys.time(),":                ",classification," ",s))
    
    # for all the relevant GAM files:
    # read the GAM and 
    
    l_GAM_data <- list()
    
    ## use the same lapply theory but over profile GAMs to build a table, which you then GAM again
    
    l_sector_data <- lapply(emis[sector==s,Profile_ID], sampledProfileGAMs, emis = emis, timescale = timescale )
    
    # rbind it
    dt_GAM_data <- rbindlist(l_sector_data, use.names = T)
    dt_GAM_data[, w := as.numeric(as.character(w))]
    
    ## Make a GAM for the whole sector, weighted by emissions contribution of NFR codes grouped by Profile_ID
    # separate if clauses for little changes in knot selection, bs term etc
    if(timestep == "hour"){
      gam_sect <- gam(N ~ s(time, bs="cc"), data = dt_GAM_data, weights = w, method="REML")
    }else if(timestep=="hourwday"){
      gam_sect <- gam(N ~ s(hour, bs="cc") + s(wday, bs="cr", k=6), data=dt_GAM_data, weights = w, method="REML")
    }else if(timestep=="wday"){
      gam_sect <- gam(N ~ s(time, bs="cr", k=6), data = dt_GAM_data, weights = w, method="REML")
    }else if(timestep=="month"){
      gam_sect <- gam(N ~ s(time, bs="cr"), data = dt_GAM_data, weights = w, method="REML")
    }else{
      gam_sect <- gam(N ~ s(time, bs="cc"), data = dt_GAM_data, weights = w, method="REML")
    }
    
    #### strip data and put into list ####
    gam_sr <- strip_rawdata(gam_sect)
    l_GAM_sectors[[s]] <- gam_sr
    
    ## go on to make coefficients table and add to another list, writing at the end. 
  
  } # end of sector loop
  
  
  # svae the list as rda
  saveRDS(gam_sr, paste0("./output/GAM_sector/",timestep,"/",species,"_",classification,"_",s,"_GAM_",timestep,".rds"))
  
  #######################################################################################
  #######################################################################################
  
  ## make a coefficients table as well:
  ## create a blank table for GAM estimation
  print(paste0(Sys.time(),":                Coefficients csv and write data..."))
  
  dt_gc <- blankFitSector(v_sectors = v_sectors, timescale = timescale)
  dt_gc[, sector := factor(sector)]
  
  # fit sector model to empty data
  dt_coeffs <- SectorCoeffs(gam_sr = gam_sr, dt_gc = dt_gc, timescale = timescale)
  if(timescale != "hourwday") setnames(dt_coeffs, "time", timescale)
  dt_coeffs[,c("fit","se.fit","lower","upper") := NULL]
  
  if(timescale=="hourwday"){
    setkeyv(dt_coeffs, c("sector","hour","wday"))
  }else{
    setkeyv(dt_coeffs, c("sector",paste(timescale)))
  }
  
  ## set filenames and write data
  
  if(((year %in% 2010:2019) & (species %in% dt_pollutants[,upper_name]))){
    filename <- paste0("GAM_",timescale,"_",classification,"_",species,"_",year)
  }else if(isFALSE(year) & (species %in% dt_pollutants[,upper_name])){
    filename <- paste0("GAM_",timescale,"_",classification,"_",species,"_allYr")
  }else if((year %in% 2010:2019) & isFALSE(species)){
    filename <- paste0("GAM_",timescale,"_",classification,"_allGas_",year)
  }else{
    filename <- paste0("GAM_",timescale,"_",classification,"_allGas_allYr")
  }
  
  dir.create(file.path(paste0("./output/GAM_sector/",classification)), showWarnings = FALSE)
  dir.create(file.path(paste0("./output/coeffs_sector/",classification)), showWarnings = FALSE)
  
  saveRDS(gam_sr, paste0("./output/GAM_sector/",classification,"/",filename,".rds"))
  fwrite(dt_coeffs, paste0("./output/coeffs_sector/",classification,"/",filename,".csv"))
  
  print(paste0(Sys.time(),": DONE."))
  
} # end of function



#########################################################################################################
#### function to plot GAM/coefficient data per sector, per year(?) per pollutant (?)
## which is better comparison? all pollutants for one year on a plot?
## or same pollutant across many years on one plot?

GAMsectorPlot <- function(year, species, classification){
  
  #### Needs data to exist to work. ####
  
  ####################################################
  
  # year            = *numeric* years to process. NA is for generic. 
  # species         = *character* name of air pollutant or GHG or metal etc. NA is for generic. 
  # classification  = *character* a sector classification system (e.g SNAP) for which GAMs have already been produced
  
  ####################################################
  
  print(paste0(Sys.time(),": Plotting..."))
  
  l_plot <- list()
  
  for(ts in c("hour","wday","yday")){
    
    l_plot_data <- list()
    
    # for(y in v_years){
    
    for(s in species){
      
      if(((year %in% 2010:2019) & (species %in% dt_pollutants[,upper_name]))){
        filename <- paste0("GAM_",timescale,"_",classification,"_",species,"_",year)
      }else if(isFALSE(year) & (species %in% dt_pollutants[,upper_name])){
        filename <- paste0("GAM_",timescale,"_",classification,"_",species,"_allYr")
      }else if((year %in% 2010:2019) & isFALSE(species)){
        filename <- paste0("GAM_",timescale,"_",classification,"_allGas_",year)
      }else{
        filename <- paste0("GAM_",timescale,"_",classification,"_allGas_allYr")
      }
      
      # read in the data
      dt <- fread(paste0("./output/coeffs_sector/",classification,"/",filename,".csv"))
      dt[,time := ts]
      dt[,Year := year]
      dt[, Gas:= s]
      setnames(dt, ts, "time_value")
      
      l_plot_data[[paste0(s,"_",year,"_",ts)]] <- dt
      
    } # pollutant loop
    
    #  } # year loop
    
    dt_plot_data <- rbindlist(l_plot_data, use.names = T)
    
    g1 <- ggplot()+
      geom_line(data=dt_plot_data, aes(x=time_value, y=coeff, group = Gas, colour=Gas))+
      geom_ribbon(data=dt_plot_data, aes(x=time_value, ymin=coeff_l, ymax=coeff_u, group = Gas, fill=Gas), alpha=0.4)+
      labs(x=ts,y="Coefficient")+
      geom_hline(yintercept = 1, linetype="dashed")+
      facet_wrap(~sector, nrow=1)+
      theme_bw()+
      theme(axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=12),
            axis.title = element_text(size=16),
            legend.title=element_blank())
    
    l_plot[[ts]] <- ggplotGrob(g1)
    
  } # time loop
  
  # stitch plots together
  
  d1 <- plot_grid(l_plot[["hour"]], l_plot[["wday"]], l_plot[["yday"]], ncol=1, nrow=3,rel_heights = c(1,1,1), rel_widths = c(1,1,1), align="h", axis="l")
  
  dir.create(file.path(paste0("./output/plots/",classification)), showWarnings = FALSE)
  
  save_plot(paste0("./output/plots/",classification,"/NOx_coeff_plots_emis",year,".png"), d1, base_height = 14, base_width = 22)
  
  print(paste0(Sys.time(),": DONE."))
  
}
