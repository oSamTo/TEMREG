## initialise the workspace for profiling annual emissions into sector profiles ##

"%!in%" <- Negate("%in%")


###############################################################################
#### series of functions to create temporally profiled emissions in the UK ####
###############################################################################

#################
## LOOKUPS etc ##
#################

## functions to set the various lookup tables needed
sector_to_ProID <- function(lookup_sectors, data_sources){
  
  # read
  dt_lu <- fread(lookup_sectors)
  dt_id <- fread(data_sources)
  
  # subset
  #dt_lu <- dt_lu[, c("NFR19","Source","Activity","Profile_ID")]
  
  # attach data source names
  dt <- dt_id[dt_lu, on = c("Profile_ID")]
    
  # remove duplicates
  dt <- dt[!duplicated(dt)]
  
  return(dt)
  
}

#ProID_to_Agg <- function(lookup_sectors, agg_sectors = c("SNAP","GNFR")){
#  
#  # read
#  dt <- fread(lookup_sectors)
#  
#  # subset
#  col_keep <- c("NFR19","Profile_ID", agg_sectors) 
#  dt <- dt[, ..col_keep]
#  
#  # remove duplicates
#  dt <- dt[!duplicated(dt)]
#  
#  return(dt)
#  
#}

#########################
#### VECTOR CREATION ####
#########################

## function to set the high level data sources as a vector
vector_sources <- function(data_sources, 
                           level = c("Profile_ID","data_source")){
  
  level <- match.arg(level)
  # read
  dt <- fread(data_sources)
  
  # extract unique data sources and return
  if(level == "Profile_ID"){
  
    v_sources <- sort(unique(dt[,Profile_ID]))
  
  }else if(level == "data_source"){
  
    v_sources <- sort(unique(dt[,data_source]))
  
  }
  
  
  return(v_sources)
  
}

##############################
### PROFILE ID GAM CREATION ##
##############################

## function to create single Profile ID GAMs that build up into the sectoral GAMs
## runs on fname, which is one target of a collection; collection (n = ~13) is all sectors in a time_scale (n = 5)
## !! Independent of Pollutant ; Activity driven GAMs !!
## Warning: domestic profile GAMs need high-memory, consider sampling       

GAM_by_ProfileID <- function(v_IDs, v_fname, 
                             time_scale = c("hour","hourwday","wday","month","yday"), y){
  
  time_scale <- match.arg(time_scale)
     
  # read data
  l_dt <- lapply(v_fname, fread)
  dt <- rbindlist(l_dt, use.names = T)
  
  # subset the data for year, if a year specific run is required. 
  # only some data has year specific info, keep all if it doesn't.
  if(y != "genYr"){
    
	# establish Profile_IDs that have the year info
	dt_have <- dt[year == y] ; v_pro_have <- dt_have[, unique(Profile_ID)]
	
	# then attach the other profiles that didn't have any year specific stuff. 
	dt <- rbindlist(list(dt[Profile_ID %!in% v_pro_have], dt_have), use.names = T)
  
  }
   
  # make a gam for each unique Profile ID
  v_gam_fname <- sapply(v_IDs, runGAM_ID, dt = dt, y = y, time_scale )
  
  # return
  return(v_gam_fname)  
  
} 

## function to run the GAM for the parent function above
runGAM_ID <- function(ID, dt, y, time_scale){

  # subset the formatted data
  dt_GAM_data <- dt[Profile_ID == ID]
  
  # standardise around the mean
  if(time_scale != "hourwday") dt_GAM_data[, N := N/mean(N, na.rm=T), by = year]
  if(time_scale == "hourwday") dt_GAM_data[, N := N/mean(N, na.rm=T), by = .(wday, year)]
  
  # particular settings for hourwday time scale
  if(time_scale != "hourwday") setnames(dt_GAM_data, paste0(time_scale), "time")
  if(time_scale == "hourwday") dt_GAM_data[,wday := factor(wday, levels=1:7)]
      
  # make GAM - some specifics per time_scale type.
  if(time_scale == "hourwday"){
     gam_ID <- gam(N ~ s(hour, bs = "cc", by = wday) + wday, 
	               data = dt_GAM_data, method = "REML") 
  }else if(time_scale == "wday"){
     gam_ID <- gam(N ~ s(time, bs = "cr", k = 6), 
	               data = dt_GAM_data, method="REML")
  }else{
    
     # until rail data is improved, it requires its own gam with k=4 (quarterly data) 
     if(ID == "RAIL_GEN"){
	   gam_ID <- gam(N ~ s(time, bs = "cc", k=4), 
	                 data = dt_GAM_data, method="REML")
	 }else{
	   gam_ID <- gam(N ~ s(time, bs = "cc"), 
	                 data = dt_GAM_data, method="REML")
	 }	 
	 
  }
  
  # strip data and save
  gam_strip <- strip_rawdata(gam_ID)
  
  # set filename, write and return fname
  fname <- paste0("data/GAM_ProfileID/",time_scale,"/gam_",ID,"_",y,"_",time_scale,".rds")
  
  saveRDS(gam_strip, fname)
  
  return(fname)
  
  
}

## function to plot the Profile ID level GAM
GAM_by_Profile_plot <- function(v_IDs, y, v_fname){
    
  # make a composite plot for every ID, with all time_scales on (not month)
  # need to do this in serial, one after the other, to avoid X11 errors 
  # (you still get X11 errors on save_plot, but it seems to work ok)
  
  l <- list()
  
  for(ID in v_IDs){
        
	## HOUR
	fname <- paste0("data/GAM_ProfileID/hour/gam_",ID,"_",y,"_hour.rds")
	gam_strip <- readRDS(fname)
	
	# set data up for plot
	dt <- data.table(Profile_ID = ID, hour = 0:23, alpha = 0)	  
	dt[, "alpha" := mgcv::predict.gam(gam_strip, newdata = list(time = hour)), by = seq_len(nrow(dt))]
		
	# plot
    gg_hour <- runPLOT_ID(dt, x_label = "hour")
	
	## HOURWDAY
	fname <- paste0("data/GAM_ProfileID/hourwday/gam_",ID,"_",y,"_hourwday.rds")
	gam_strip <- readRDS(fname)
	
	# set data up for plot
	dt <- data.table(Profile_ID = ID, wday = rep(1:7, each = 24), hour = rep(0:23, 7), alpha = 0)
	dt[, wday := factor(wday, levels = 1:7)]
	dt[, "alpha" := mgcv::predict.gam(gam_strip, newdata = list(wday = wday, hour = hour)), by = seq_len(nrow(dt))]
	
	# plot
    gg_hourwday <- runPLOT_ID(dt, x_label = "hourwday")
		
	## WDAY
	fname <- paste0("data/GAM_ProfileID/wday/gam_",ID,"_",y,"_wday.rds")
	gam_strip <- readRDS(fname)
	
	# set data up for plot
	dt <- data.table(Profile_ID = ID, wday = 1:7, alpha = 0)	  
	dt[, "alpha" := mgcv::predict.gam(gam_strip, newdata = list(time = wday)), by = seq_len(nrow(dt))]
		
	# plot
    gg_wday <- runPLOT_ID(dt, x_label = "wday")
		
	## YDAY
	fname <- paste0("data/GAM_ProfileID/yday/gam_",ID,"_",y,"_yday.rds")
	gam_strip <- readRDS(fname)
	
	# set data up for plot
	dt <- data.table(Profile_ID = ID, yday = 1:365, alpha = 0)	  
	dt[, "alpha" := mgcv::predict.gam(gam_strip, newdata = list(time = yday)), by = seq_len(nrow(dt))]
		
	# plot
    gg_yday <- runPLOT_ID(dt, x_label = "yday")
	
	d1 <- plot_grid(gg_hour, gg_hourwday, gg_wday, gg_yday, ncol=1, nrow=4,rel_heights = c(1.5,1,1,1), rel_widths = c(1,1,1,1), align="h", axis="l")
    
	fname <- paste0("./data/GAM_ProfileID/a_plots/plot_",ID,"_",y,".png")
	
    save_plot(fname, d1, base_height = 14, base_width = 14)
	
	l[[ID]] <- fname
	
		    
   } # profile loop
        
	v_fname <- unlist(l, use.names=F)
	
	return(v_fname)
	
} 
  
## function to plot Profile ID GAMs for parent function above
runPLOT_ID <- function(dt, x_label){
  
  if(x_label != "hourwday"){
  
     g1 <- ggplot()+
           geom_line(data = dt, aes(x = get(x_label), y = alpha), colour="red")+
		   labs(x = x_label, y = "Coefficient")+
		   theme_bw()+
		     theme(axis.text.x = element_text(size=12),
             axis.text.y = element_text(size=12),
             axis.title = element_text(size=16),
             legend.title=element_blank())	
   
  }else{
  
     g1 <- ggplot()+
           geom_line(data = dt, aes(x = hour,y = alpha, group = wday, colour=wday))+
		   labs(x = "hour", y = "Coefficient")+
		   theme_bw()+
		     theme(axis.text.x = element_text(size=12),
             axis.text.y = element_text(size=12),
             axis.title = element_text(size=16),
             legend.title=element_blank())	
            
  }
   
  return(g1)  
 
}

##############################
### AGG-SECTOR GAM CREATION ##
##############################

## function to create the sector GAMs
## uses Profile IDs and NFR contributions to sector totals to weight those GAMs
GAM_by_Agg_sector <- function(v_fname, v_fname_HE, y, pollutant, 
                              time_scale = c("hour","hourwday","wday","month","yday"), 
							  classification = c("SNAP","GNFR"), dt_lookup, dt_PID,
							  inv_year){
  
 # read in inventory_processor info
  
 # sector classification - SNAP, GNFR, other (?)
 # year - this is "genYr" or e.g. 2020
 # time_scale - could have all and do lapply?
    
  classification <- match.arg(classification)
  timescale <- match.arg(time_scale)
    
  # fetch NAEI data and match to sectors etc. (write for posterity)
  dt_emis <- join_NAEI_to_Profiles(y, pollutant, classification, dt_lookup, dt_PID, inv_year)
  dt_emis[, year := y]
  fwrite(dt_emis, paste0("output/profile_shares/Pro_ID_percent_",pollutant,"_",y,".csv"))
  
  # change this to read from the emissions/sector lookup - call the profileNAEI function. 
  v_sectors <- dt_lookup[,unique(get(classification))]
  if(classification == "SNAP") v_sectors <- v_sectors[ !v_sectors == 99]
  
      
  l_gam <- lapply(X = v_sectors, FUN = create_sector_gam, y = y,
                       dt_emis = dt_emis, classification = classification, time_scale = time_scale,
					   v_fname = v_fname, v_fname_HE = v_fname_HE, n_samples = 100)
  
  names(l_gam) <- as.character(v_sectors)
  
  ## create tables of central estimates for csv outputs
  l_csv <- list()
  
  for(s in v_sectors){
       
    # go on to make coefficients table and add to another list, writing at the end. 
	dt_central <- blank_response(time_scale, n_samples = 1, n_steps = "full")
	  
    # fit sector model to empty data
	if(time_scale != "hourwday"){
    
	dt_central[, N := mgcv::predict.gam(l_gam[[as.character(s)]], newdata = list(time = time))]
	dt_central[, N := N/mean(N)] # tiny fluctuations from rounded sum. Adjust back.
    	  
    } else {
  
    dt_central[, N := mgcv::predict.gam(l_gam[[as.character(s)]], newdata = list(wday = wday, hour = hour))]
    dt_central[, N := N/mean(N), by = "wday"]
	
    }
	
  l_csv[[paste0(s,"_csv")]] <- dt_central
    
  } # end of sector loop for csvs 
  
  l_results <- list(l_gam, l_csv)
  
  
  # write out the sector GAMs and the central estimate data tables - pass to write function
  l_fname <- lapply(X = v_sectors, write_sector_gam, l_results, y, pollutant, classification, time_scale)
   
  v_fname <- unlist(l_fname)
  
  return(v_fname)  
  
} # end of function

## function to 1. return the NAEI emissions data, formatted, for the given year & Species. 
##             2. match the annual NAEI data to the temporal profile codes and GAM files. 

join_NAEI_to_Profiles <- function(y, pollutant, classification = c("SNAP","GNFR"), dt_lookup, dt_PID, inv_year){
  
  ########################################################
  
  # year           = *numeric* year to process. Will determine calendar structure plus year specific profiles.
  #if(!is.numeric(v_year)) stop ("Year is not numeric")
  # species        = *character* name of air pollutant or GHG or metal etc. Needs to conform to a list of options.
  #if(pollutant %!in% c("ch4","co2","n2o","bap","bz","hcl","nox","so2","nh3", "co", "voc","cd","cu","pb","hg","ni","zn", "pm0_1","pm1","pm2_5","pm10", NA) ) stop ("Species must be in: 
  #                                          AP:    bap, bz, co, hcl, nh3, nox, so2, voc
  #                                          PM:    pm0_1, pm1, pm2_5, pm10
  #                                          GHG:   ch4, co2, n2o
  #                                          Metal: cd, cu, hg, ni, pb, zn
  #                                          (Or NA)")
  
  classification <- match.arg(classification)
  
  #########################################################
  
  ## adjustment to lookup table ##  
  dt_lookup[NFR19 == "4.00E+01", NFR19 := "4E1"] # how to fix this?? options(scipen=999) not working
  dt_lookup[NFR19 == "4.00E+02", NFR19 := "4E2"] # how to fix this?? options(scipen=999) not working
  dt_lookup[NFR19 == "40", NFR19 := "4E1"] # how to fix this?? options(scipen=999) not working
  dt_lookup[NFR19 == "400", NFR19 := "4E2"] # how to fix this?? options(scipen=999) not working
  
  # read data from inventory processor, subset accordingly
  dt_naei <- fread(paste0("../../inventory_processor/data/NAEI/inv",inv_year,"/totals/NAEI_AllPoll_TOTALS_inv",inv_year,"_emis_1970-",inv_year-2,"_NFR_t.csv"))
  
  # keep all pollutants if making an 'generic' pollutant 
  if(pollutant == "genPoll") dt_naei <- dt_naei[Pollutant %in% c("nox","nh3","co", "pm25", "sox", "pmco", "nmvoc")]
  if(pollutant != "genPoll") dt_naei <- dt_naei[Pollutant == pollutant]
  
  # keep specific year, or last 5 years if year is 'genYr'
  if(y == "genYr") dt_naei <- dt_naei[Year %in% (max(Year) - 5):max(Year)]
  if(y != "genYr") dt_naei <- dt_naei[Year == y]
    
  # at the moment, road re-suspension (z_11C )is set to MISC_ANY and is dominating PM and making it flat(ish)
  # when doing pollutant generic, it effects that also due to the mean contribution to SN07 being high
  # it comes from all sources of traffic, so by removing it, we are maintaining it within the other sectors?
  dt_naei <- dt_naei[!(NFR19 == "z_11C" & Source == "Road transport - resuspension")]
  
  # join the NAEI emissions to the classification lookup and to the GAM name lookup
  dt_joined <- dt_lookup[dt_naei, on = c("NFR19","Source","Activity")]
  dt_joined <- dt_joined[!is.na(emis_t)]
    
  # reformat and sum by NFR19/Profile.
  cols_keep <- c("Pollutant","NFR19",classification,"Profile_ID","Year","emis_t")
  dt_joined <- dt_joined[,..cols_keep]
  setnames(dt_joined, classification, "sector")
      
  # Calculate the Profile_ID fraction of sector, and take a mean over all pollutants and years
  dt_joined_agg <- dt_joined[, .(emis_t = sum(emis_t, na.rm=T)), by=.(Pollutant, sector, Profile_ID, Year)]
  dt_joined_agg[, emis_sec_frac := emis_t/sum(emis_t), by= .(Pollutant, sector, Year)]
  
  dt_joined_agg <- dt_joined_agg[!(is.nan(emis_sec_frac))]
  dt_joined_agg <- dt_joined_agg[!(is.na(emis_sec_frac))]
  
  # take the mean contribution by each Profile to each sector. Due to variability in years and pollutants, re-adjust to 1. 
  dt_mean_frac <- dt_joined_agg[, .(emis_sec_frac = mean(emis_sec_frac)), by = .(sector, Profile_ID)]
  dt_mean_frac[, emis_sec_frac := emis_sec_frac/sum(emis_sec_frac), by = sector]
  
  if(classification=="SNAP") dt_mean_frac[, sector := as.numeric(sector)]
  
  setnames(dt_mean_frac, "sector", classification)
    
  return(dt_mean_frac)
  
}

## function to return a sector-level GAM, using weightings of sector contribution per Profile_ID.
create_sector_gam <- function(i_sector, y = y, dt_emis, classification, time_scale, v_fname, v_fname_HE, n_samples){

  # for all the relevant GAM files:
  ## use the same lapply theory but over profile GAMs to build a table,
  ## which you then GAM again
  v_ids <- dt_emis[get(classification) == i_sector, unique(Profile_ID)]
  
  dt_blank <- blank_response(time_scale, n_samples, n_steps = "half")

  if (length(v_ids) == 0) {
    
	## create a 'flat' GAM if there are no representative sub-sector GAMs
    dt_blank[, N := sample(seq(1 * ((100 - 2) / 100), 1 * ((100 + 2) / 100), 0.005), 1), 
         by = seq_len(nrow(dt_blank))]
  
    # make GAM - some specifics per time_scale type.
    if(time_scale == "hourwday"){
       m_gam_i <- gam(N ~ s(hour, bs = "cc", by = wday) + wday, 
	                  data = dt_blank, method = "REML") 
    }else if(time_scale == "wday"){
       m_gam_i <- gam(N ~ s(time, bs = "cr", k = 6), 
	                  data = dt_blank, method="REML")
    }else{
       m_gam_i <- gam(N ~ s(time, bs = "cc"), 
	                  data = dt_blank, method="REML") 	 	 
    }
	
  } else {
  
    # normal dist. sample across all Profile_ID gams and add on the error term as well
    l_dt <- lapply(X = v_ids, FUN =  sample_profile_gam, i_sector = i_sector, y = y, v_fname = v_fname,
	                                 dt = dt_blank, classification = classification, dt_emis = dt_emis, time_scale = time_scale )
		
    # rbind it into one table.
    dt_gam_data <- rbindlist(l_dt, use.names = TRUE)

    # take the weighted mean for values of every sample draw, at each t
    if(time_scale != "hourwday") dt_gam_data <- dt_gam_data[, .(N = weighted.mean(N, w)), by = .(time, draw)]
	if(time_scale == "hourwday") dt_gam_data <- dt_gam_data[, .(N = weighted.mean(N, w)), by = .(wday, hour, draw)]

    ## Make a GAM for the whole sector
    # separate if clauses for little changes in knot selection, bs term, time_scale etc
	if(time_scale == "hourwday"){
      m_gam_i <- gam(N ~ s(hour, bs = "cc", by = wday) + wday, 
	                  data = dt_gam_data, method = "REML") 
    }else if(time_scale == "wday"){
      m_gam_i <- gam(N ~ s(time, bs = "cr", k = 6), 
	                  data = dt_gam_data, method="REML")
    }else{
      m_gam_i <- gam(N ~ s(time, bs = "cc"), 
	                  data = dt_gam_data, method="REML") 	 	 
    }
	
  }

  return(m_gam_i)

}                    

blank_response <- function(time_scale, n_samples, n_steps = c("half","full")){
  
  if (time_scale == "yday") {
    if(n_steps == "half"){
	    time_steps <- seq(1, 365, by = 2)
	}else{
	    time_steps <- 1:365
	}
  } else if (time_scale == "month") {
    time_steps <- 1:12
  } else if (time_scale == "wday") {
    time_steps <- 1:7
  } else if (time_scale == "hour") {
    time_steps <- 0:23
  } else {
    time_steps <- 0:23
  }
  
  # make blank dt  
  if(time_scale != "hourwday"){
    dt <- data.table( time = rep(time_steps, each = n_samples))
  }else if(time_scale=="hourwday"){
    dt <- data.table( wday = rep(rep(1:7, each = 24), each = n_samples), hour = rep(rep(0:23,7), each = n_samples))
  }
  
  if(time_scale == "hourwday") dt[,wday := factor(wday, levels=1:7)]
    
  return(dt)  
}

sample_profile_gam <- function(ID, i_sector, y, v_fname, dt, classification, dt_emis, time_scale){

  # collect GAM
  if(ID %in% c("ROAD_CAT_URB", "ROAD_CAT_MOT", "ROAD_CAT_RUR", 
                "ROAD_LGV_URB", "ROAD_LGV_MOT", "ROAD_LGV_RUR", 
	            "ROAD_HGV_URB", "ROAD_HGV_MOT", "ROAD_HGV_RUR") & y != "genYr" & time_scale %in% c("month", "yday")){
				
				m_gam_ID <- readRDS(paste0("../TEMREG_data/Data/RoadTransport/HE_API/HE_ProID_GAMs/gam_",ID,"_",y,"_",time_scale,".rds"))
  }else{
				m_gam_ID <- readRDS(v_fname[grep(paste0("_",ID,"_",y,"_",time_scale), v_fname)])
  }
  
  dt_gam <- copy(dt)
  
  # sample normally across prediction interval
  if(time_scale != "hourwday"){
    
	dt_gam[, N := mgcv::predict.gam(m_gam_ID, newdata = list(time = time))
            + rnorm(nrow(dt_gam), 0, sqrt(mgcv::predict.gam(m_gam_ID,
            newdata = list(time = time),
            se.fit = TRUE)$se.fit^2 + m_gam_ID$sig2)) ]
    
	# add in weights, draw number
    dt_gam[, w := dt_emis[get(classification) == i_sector & Profile_ID == ID, emis_sec_frac]]
    dt_gam[, draw := 1:.N, by = time]
	  
  } else {
  
    dt_gam[, N := mgcv::predict.gam(m_gam_ID, newdata = list(wday = wday, hour = hour))
            + rnorm(nrow(dt_gam), 0, sqrt(mgcv::predict.gam(m_gam_ID,
            newdata = list(wday = wday, hour = hour),
            se.fit = TRUE)$se.fit^2 + m_gam_ID$sig2)) ]
    
	# add in weights, draw number
    dt_gam[, w := dt_emis[get(classification) == i_sector & Profile_ID == ID, emis_sec_frac]]
    dt_gam[, draw := 1:.N, by = .(wday, hour)] 
  
  }
  
  return(dt_gam)      
	 
}

## function to write out all the aggregated sector GAM information
write_sector_gam <- function(i_sector, l_results, y, pollutant, classification, time_scale){

  # if SNAP, pad the SNAP
  if(classification == "SNAP"){
    fname_sector <- str_pad(i_sector, 2, "0", side = "left")
  }else{
    fname_sector <- i_sector
  }

  ## the first part of the list, is a list of sector GAMs
  i_gam <- l_results[[1]][[as.character(i_sector)]]
  
  # set up folder
  folname <- paste0("output/GAM_sector/",classification,"/",pollutant,"/",y)
  suppressWarnings(dir.create(folname, recursive = T))
  
  # filename
  fname_gam <- paste0(folname,"/",classification,"_",fname_sector,"_",time_scale,"_",pollutant,"_",y,".rds")    
  
  # save as RDS
  saveRDS(i_gam, fname_gam)
    
  ## the second part of the list, is a list of tables of central estimates from GAMs
  dt_ce <- l_results[[2]][[paste0(as.character(i_sector),"_csv")]]
  
  # set up folder
  folname <- paste0("output/coeff_sector/",classification,"/",pollutant,"/",y)
  suppressWarnings(dir.create(folname, recursive = T))
  
  # filename
  fname_csv <- paste0(folname,"/",classification,"_",fname_sector,"_",time_scale,"_",pollutant,"_",y,".csv")    
  
  #save as RDS
  fwrite(dt_ce, fname_csv)
  
  return(fname_gam)

}


## function to plot the aggregated sector GAM
GAM_agg_plot <- function(y, pollutant, classification = c("SNAP","GNFR"), dt_lookup, dt_PID, inv_year, v_fname){
      
  # at present, this does NOT use the v_fname passed as an argument, it simply 
  # reads from drive, but it is there to prompt the function to run if needed.
	
  # make a composite plot for the classification (all sectors), with all time_scales on (not month)
  # need to do this in serial, one after the other, to avoid X11 errors 
  # (you still get X11 errors on save_plot, but it seems to work ok)
    
  # make a sector list based on the classification, and the pollutant (do NAEI matching again)
  dt_emis <- join_NAEI_to_Profiles(y, pollutant, classification, dt_lookup, dt_PID, inv_year)
    
  # change this to read from the emissions/sector lookup - call the profileNAEI function. 
  #v_sectors <- dt_emis[,unique(get(classification))]
  v_sectors <- dt_lookup[,unique(get(classification))]
  if(classification == "SNAP") v_sectors <- v_sectors[ !v_sectors == 99]
  if(classification == "SNAP") v_sectors <- str_pad(v_sectors, 2, 0, side = "left")
    	
  # lists for plots
  l_hour     <- list()
  l_hourwday <- list()
  l_wday     <- list()
  l_yday     <- list()
  
  for(v in v_sectors){
        
	## HOUR
	fname <- paste0("output/GAM_sector/",classification,"/",pollutant,"/",y,"/",classification,"_",v,"_hour_",pollutant,"_",y,".rds")
	gam_strip <- readRDS(fname)
	
	# set data up for plot
	dt <- data.table(Class = classification, Sector = v, hour = 0:23, alpha = 0)	  
	dt[, "alpha" := mgcv::predict.gam(gam_strip, newdata = list(time = hour)), by = seq_len(nrow(dt))]
    dt[, alpha := alpha/mean(alpha)]
	
	l_hour[[v]] <- dt
		
	## HOURWDAY
	fname <- paste0("output/GAM_sector/",classification,"/",pollutant,"/",y,"/",classification,"_",v,"_hourwday_",pollutant,"_",y,".rds")
	gam_strip <- readRDS(fname)
	
	# set data up for plot
	dt <- data.table(Class = classification, Sector = v, wday = rep(1:7, each = 24), hour = rep(0:23, 7), alpha = 0)
	dt[, wday := factor(wday, levels = 1:7)]
	dt[, "alpha" := mgcv::predict.gam(gam_strip, newdata = list(wday = wday, hour = hour)), by = seq_len(nrow(dt))]
	dt[, alpha := alpha/mean(alpha)]
	
	l_hourwday[[v]] <- dt
	
	## WDAY
	fname <- paste0("output/GAM_sector/",classification,"/",pollutant,"/",y,"/",classification,"_",v,"_wday_",pollutant,"_",y,".rds")
	gam_strip <- readRDS(fname)
	
	# set data up for plot
	dt <- data.table(Class = classification, Sector = v, wday = 1:7, alpha = 0)	  
	dt[, "alpha" := mgcv::predict.gam(gam_strip, newdata = list(time = wday)), by = seq_len(nrow(dt))]
	dt[, alpha := alpha/mean(alpha)]
	
	l_wday[[v]] <- dt
		
	## YDAY
	fname <- paste0("output/GAM_sector/",classification,"/",pollutant,"/",y,"/",classification,"_",v,"_yday_",pollutant,"_",y,".rds")
	gam_strip <- readRDS(fname)
	
	# set data up for plot
	dt <- data.table(Class = classification, Sector = v, yday = 1:365, alpha = 0)	  
	dt[, "alpha" := mgcv::predict.gam(gam_strip, newdata = list(time = yday)), by = seq_len(nrow(dt))]
	dt[, alpha := alpha/mean(alpha)]
		
	l_yday[[v]] <- dt	
		    
  } # profile loop
   
  # combine to one table and pass to plot function. 
  dt_hour     <- rbindlist(l_hour,     use.names = T)
  dt_hourwday <- rbindlist(l_hourwday, use.names = T)
  dt_wday     <- rbindlist(l_wday,     use.names = T)
  dt_yday     <- rbindlist(l_yday,     use.names = T)
   
  gg_hour     <- runPLOT_agg(dt = dt_hour,     x_label = "hour")
  gg_hourwday <- runPLOT_agg(dt = dt_hourwday, x_label = "hourwday")
  gg_wday     <- runPLOT_agg(dt = dt_wday,     x_label = "wday")
  gg_yday     <- runPLOT_agg(dt = dt_yday,     x_label = "yday")
   
  # make one composite plot	
  d1 <- plot_grid(gg_hour, gg_hourwday, gg_wday, gg_yday, ncol=1, nrow=4,rel_heights = c(1,1,1,1), rel_widths = c(1,1,1,1), align="h", axis="l")
    
  # set folder & filename
  folname <- paste0("output/plots/",classification,"/",pollutant)
  suppressWarnings(dir.create(folname, recursive = T))
    
  fname <- paste0(folname,"/GAM_",classification,"_",pollutant,"_",y,"_plots.png")
  
  # write and return
  save_plot(fname, d1, base_height = 20, base_width = 14)	
	  
  return(fname)
	
} 
  
## function to plot Profile ID GAMs for parent function above
runPLOT_agg <- function(dt, x_label){
  
  if(x_label != "hourwday"){
  
     g1 <- ggplot()+
           geom_line(data = dt, aes(x = get(x_label), y = alpha), colour = "red")+
		   labs(x = x_label, y = "Coefficient")+
		   facet_wrap(~Sector)+
		   theme_bw()+
		     theme(axis.text.x = element_text(size=12),
             axis.text.y = element_text(size=12),
             axis.title = element_text(size=16),
             legend.title=element_blank())	
   
  }else{
  
     g1 <- ggplot()+
           geom_line(data = dt, aes(x = hour, y = alpha, group = wday, colour = wday))+
		   labs(x = "hour", y = "Coefficient")+
		   facet_wrap(~Sector)+
		   theme_bw()+
		     theme(axis.text.x = element_text(size=12),
             axis.text.y = element_text(size=12),
             axis.title = element_text(size=16),
             legend.title=element_blank())	
            
  }
   
  return(g1)  
 
}

#################################
### MODEL-READY INPUT CREATION ##
#################################

## function to create the temporal input files for models
# only changing for countries nominated, at the moment always UK
# current understanding: 
   # Month and wday files: one file for each pollutant, SNAP code & country (however modern EMEP has monthly NCDFs)
   # hourly files: one file. Each country & SNAP across each hour; pollutant generic

EMEP4UK_profiles <- function(y, pollutant, v_iso = c(27), classification = "SNAP",
                             v_fname, time_scale = c("hourwday","wday","month")){
  
  # at present, this does NOT use the v_fname passed as an argument, it simply 
  # reads from drive, but it is there to prompt the function to run if needed.
  
  # need an emmep pollutanr variable, just for VOCs
  if(pollutant == "voc"){
    temreg_pollutant <- "nmvoc"
  }else{
    temreg_pollutant <- pollutant
  }
  
  classification <- match.arg(classification)
  #v_iso <- match.arg(v_iso)
  time_scale <- match.arg(time_scale)
  
  # set folder name
  folname <- paste0("output/model_inputs/EMEP4UK/",y)
  suppressWarnings(dir.create(folname, recursive = T))
  
  if(time_scale == "hourwday"){
  
  ## HOURLY ##
  # one file. hourly coeffs by wday & SNAP; pollutant generic, country generic
  dt_model_old <- fread("output/model_inputs/EMEP4UK/pre_TEMREG/HourlyFacs.INERIS", skip = 3)
  names(dt_model_old) <- c("wday","SNAP",1:24)
 
  # bring in the pollutant generic hourwday file
  l <- list()
  for(s in 1:11){
  
    dt <- fread(paste0("output/coeff_sector/SNAP/genPoll/",y,"/SNAP_",str_pad(s, 2, "0", side = "left"),"_hourwday_genPoll_",y,".csv"))
    dt[,SNAP := str_pad(s, 2, "0", side = "left")]
	dt[, N := format(round(N, 2), nsmall = 2) ]
	
	dtw <- dcast(dt, wday + SNAP ~ hour, value.var = "N")
	l[[s]] <- dtw
  
  }
  
  dt_hourwday <- rbindlist(l)
  keycol <- c("wday","SNAP")
  setorderv(dt_hourwday, keycol)
    
  # make a header, append 
  header_1 <- paste0("# Hourly factors generated by TEMREG for 'generic' pollutant, for ", y)
  header_2 <- "# day snap f1 f2 ... f24  (note f1 from 00:00 to 00:59)"
  header_3 <- "#"
  dt_header <- rbind(header_1, header_2, header_3)
  
  # set filenames, write the header and write/append the data
  fname <- paste0("HourlyFacs.INERIS")
  fwrite(dt_header, paste0(folname,"/",fname), col.names = F)
  fwrite(dt_hourwday, paste0(folname,"/",fname), append = T, sep = " ")
      
  }else if(time_scale == "wday"){
  
  ## DAILY ##
  # one file per pollutant. Daily coeffs, per SNAP and per country
  dt_model_old <- fread(paste0("output/model_inputs/EMEP4UK/pre_TEMREG/DailyFac.", pollutant))
  names(dt_model_old) <- c("iso","SNAP",1:7)
  
  for(i in v_iso){ # this is wrong, it will keep bringing in the old data and undoing the changes made in the previous loop
    
    # bring in the pollutant specific wday file
    l <- list()
    for(s in 1:11){  
      
	  dt <- fread(paste0("output/coeff_sector/SNAP/",temreg_pollutant,"/",y,"/SNAP_",
	                     str_pad(s, 2, "0", side = "left"),"_wday_",temreg_pollutant,"_",y,".csv"))
      #dt[,SNAP := str_pad(s, 2, "0", side = "left")]
	  dt[,SNAP := s]
	  dt[, iso := i]
	  
	  l[[s]] <- dt
	  rm(dt)
	  gc()
    }
  
    dt_wday_m <- rbindlist(l)
      
  	# add the (melted) original data to the new data (minus the iso of interest)
	# we are metling because making the new data wide, and rbindlist() after, messes up all the formatting
	dt_model_old_m <- melt(dt_model_old[iso != i], id.vars = c("iso","SNAP"), variable.name = "time", value.name = "N")
	dt_wday_m_bind <- rbindlist(list(dt_wday_m, dt_model_old_m), use.names = T)
    dt_wday_m_bind[, N := format(round(N, 3), nsmall = 3) ]
	
	dt_wday <- dcast(dt_wday_m_bind, iso + SNAP ~ time, value.var = "N")
	keycol <- c("iso","SNAP")
    setorderv(dt_wday, keycol)
  
  }
  
  # set filenames, write the header and write/append the data
  fname <- paste0("DailyFac.", pollutant)
  fwrite(dt_wday, paste0(folname,"/",fname), col.names = F, sep = " ")
    
  }else if(time_scale == "month"){
  
  ## MONTHLY ##
  # one file per pollutant. monthly coeffs, per SNAP and per country
  # pre-TEMREG monthly coeffs only have SNAPs 1:4 & 10, replacing with full set
  dt_model_old <- fread(paste0("output/model_inputs/EMEP4UK/pre_TEMREG/MonthlyFacs.", pollutant))
  names(dt_model_old) <- c("iso","SNAP",1:12)
    
  for(i in v_iso){ # this is wrong, it will keep bringing in the old data and undoing the changes made in the previous loop
    
    # bring in the pollutant generic hourwday file
    l <- list()
    for(s in 1:11){
  
      dt <- fread(paste0("output/coeff_sector/SNAP/",temreg_pollutant,"/",y,"/SNAP_",
	                      str_pad(s, 2, "0", side = "left"),"_month_",temreg_pollutant,"_",y,".csv"))
      #dt[,SNAP := str_pad(s, 2, "0", side = "left")]
	  dt[,SNAP := s]
	  dt[, iso := i]
	  
	  l[[s]] <- dt
	  rm(dt)
	  gc()
    }
  
    dt_month_m <- rbindlist(l)
      
  	# add the (melted) original data to the new data (minus the iso of interest)
	# we are metling because making the new data wide, and rbindlist() after, messes up all the formatting
	#str(melt(dt_model_old[iso != i], id.vars = c("iso","SNAP"), variable.name = "time", value.name = "N"))
	dt_model_old_m <- melt(dt_model_old[iso != i], id.vars = c("iso","SNAP"), variable.name = "time", value.name = "N")
	dt_month_m_bind <- rbindlist(list(dt_month_m, dt_model_old_m), use.names = T)
    dt_month_m_bind[, N := format(round(N, 3), nsmall = 3) ]
	#dt_month_m[, N := as.numeric(N)]
	dt_month_m_bind[, SNAP := str_pad(SNAP, 2, "0", side = "left")]
	#dt_month_m[, SNAP := as.numeric(SNAP)]
	
	dt_month <- dcast(dt_month_m_bind, iso + SNAP ~ time, value.var = "N")
	keycol <- c("iso","SNAP")
    setorderv(dt_month, keycol)
  
  }
  
  # set filenames, write the header and write/append the data
  fname <- paste0("MonthlyFacs.", pollutant)
  fwrite(dt_month, paste0(folname,"/",fname), col.names = F, sep = " ")
    
  }
    
  return(paste0(folname,"/",fname)) 
    
} # end of function

## function to plot old vs new profiles
# UK plots only
EMEP4UK_profiles_plot <- function(y, classification = "SNAP",
                                  time_scale = c("hourwday","wday","month"), v_fname){
  
   if(time_scale == "hourwday"){
	  
    dt_old <- fread("output/model_inputs/EMEP4UK/pre_TEMREG/HourlyFacs.INERIS", skip = 3)
    names(dt_old) <- c("wday","SNAP",1:24)
	dt_old_m <- melt(dt_old, id.vars = c("wday", "SNAP"), variable.name = "time", value.name = "alpha")
		
    dt_new <- fread(paste0("output/model_inputs/EMEP4UK/",y,"/HourlyFacs.INERIS"), skip = 3)
    names(dt_new) <- c("wday","SNAP",1:24)
	dt_new_m <- melt(dt_new, id.vars = c("wday", "SNAP"), variable.name = "time", value.name = "alpha")
	dt_new_m[, alpha := as.numeric(alpha)]
	
	# combine
	dt_old_m[, dataset := "pre_TEMREG"]
	dt_new_m[, dataset := paste0("TEMREG_",y)]
	dt <- rbindlist(list(dt_old_m, dt_new_m), use.names = T)
		
	# plot
	g1 <- ggplot()+
          geom_line(data = dt, aes(x = time, y = alpha, group = dataset, colour = dataset))+
		  scale_color_manual(values=c("black","red"))+
	      labs(x = "t", y = "Coefficient")+
		  facet_grid(SNAP~wday)+
		  theme_bw()+
		    theme(axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=12),
            axis.title = element_text(size=16),
            legend.title = element_blank(),
		    legend.text = element_text(size=16))	
    
    # name and save
    fname <- paste0("output/model_inputs/EMEP4UK/",y,"/HourlyFacs_plot.png")
    ggsave(fname, g1, width = 20, height = 16)
		
   }else if(time_scale == "wday"){
		
    fnames_old <- list.files("output/model_inputs/EMEP4UK/pre_TEMREG", pattern = "^DailyFac", full.names = T)	
	l_old <- lapply(fnames_old, fread)
	l_old <- lapply(l_old, function(x) setNames(x, c("iso","SNAP",1:7)) )
	l_old <- lapply(1:length(l_old), function(x) l_old[[x]][, pollutant := str_split(fnames_old[x],"\\.")[[1]][2]])
	dt_old <- rbindlist(l_old, use.names=T)[iso == 27]
	dt_old[, iso := NULL]
	dt_old_m <- melt(dt_old, id.vars = c("pollutant", "SNAP"), variable.name = "time", value.name = "alpha")
	
	fnames_new <- list.files(paste0("output/model_inputs/EMEP4UK/",y), pattern = "^DailyFac\\.", full.names = T)	
	l_new <- lapply(fnames_new, fread)
	l_new <- lapply(l_new, function(x) setNames(x, c("iso","SNAP",1:7)) )
	l_new <- lapply(1:length(l_new), function(x) l_new[[x]][, pollutant := str_split(fnames_new[x],"\\.")[[1]][2]])
	dt_new <- rbindlist(l_new, use.names=T)[iso == 27]
	dt_new[, iso := NULL]
	dt_new_m <- melt(dt_new, id.vars = c("pollutant", "SNAP"), variable.name = "time", value.name = "alpha")
	dt_new_m[, alpha := as.numeric(alpha)]
	
	# combine
	dt_old_m[, dataset := "pre_TEMREG"]
	dt_new_m[, dataset := paste0("TEMREG_",y)]
	dt <- rbindlist(list(dt_old_m, dt_new_m), use.names = T)
		
	# plot
	g2 <- ggplot()+
          geom_line(data = dt, aes(x = time, y = alpha, group = dataset, colour = dataset))+
		  scale_color_manual(values=c("black","red"))+
	      labs(x = "t", y = "Coefficient")+
		  facet_grid(SNAP~pollutant)+
		  theme_bw()+
		    theme(axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=12),
            axis.title = element_text(size=16),
            legend.title = element_blank(),
		    legend.text = element_text(size=16))	
    
    # name and save
    fname <- paste0("output/model_inputs/EMEP4UK/",y,"/DailyFacs_plot.png")
    ggsave(fname, g2, width = 20, height = 16)
	  
   }else if(time_scale == "month"){
	
	fnames_old <- list.files("output/model_inputs/EMEP4UK/pre_TEMREG", pattern = "^MonthlyFacs.", full.names = T)	
	l_old <- lapply(fnames_old, fread)
	l_old <- lapply(l_old, function(x) setNames(x, c("iso","SNAP",1:12)) )
	l_old <- lapply(1:length(l_old), function(x) l_old[[x]][, pollutant := str_split(fnames_old[x],"\\.")[[1]][2]])
	dt_old <- rbindlist(l_old, use.names=T)[iso == 27]
	dt_old[, iso := NULL]
	dt_old_m <- melt(dt_old, id.vars = c("pollutant", "SNAP"), variable.name = "time", value.name = "alpha")
	
	fnames_new <- list.files(paste0("output/model_inputs/EMEP4UK/",y), pattern = "^MonthlyFacs\\.", full.names = T)	
	l_new <- lapply(fnames_new, fread)
	l_new <- lapply(l_new, function(x) setNames(x, c("iso","SNAP",1:12)) )
	l_new <- lapply(1:length(l_new), function(x) l_new[[x]][, pollutant := str_split(fnames_new[x],"\\.")[[1]][2]])
	dt_new <- rbindlist(l_new, use.names=T)[iso == 27]
	dt_new[, iso := NULL]
	dt_new_m <- melt(dt_new, id.vars = c("pollutant", "SNAP"), variable.name = "time", value.name = "alpha")
	dt_new_m[, alpha := as.numeric(alpha)]
	
	# combine
	dt_old_m[, dataset := "pre_TEMREG"]
	dt_new_m[, dataset := paste0("TEMREG_",y)]
	dt <- rbindlist(list(dt_old_m, dt_new_m), use.names = T)
		
	# plot
	g3 <- ggplot()+
          geom_line(data = dt, aes(x = time, y = alpha, group = dataset, colour = dataset))+
		  scale_color_manual(values=c("black","red"))+
	      labs(x = "t", y = "Coefficient")+
		  facet_grid(SNAP~pollutant)+
		  theme_bw()+
		    theme(axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=12),
            axis.title = element_text(size=16),
            legend.title = element_blank(),
		    legend.text = element_text(size=16))	
    
    # name and save
    fname <- paste0("output/model_inputs/EMEP4UK/",y,"/MonthlyFacs_plot.png")
    ggsave(fname, g3, width = 20, height = 16)
	
   }
  
    
 
 }







hide <- function(){

  require(targets)
  require(terra)
  require(stringr)
  require(data.table)
  require(readxl)
  require(vctrs)
  require(mgcv)
  require(metagam)
  require(ggplot2)
  require(cowplot)
  
  source("R/TEMREG_functions.R")  
    
  
  #dt_NFR <- as.data.table(read_excel("data/lookups/sector_lookups.xlsx", sheet = "NAEI_lookup"))
  #dt_PID <- as.data.table(read_excel("data/lookups/pollutants.xlsx", sheet = "NAEI_pollutants"))
  #dt_SIC <- fread("data/lookups/points_sectors_to_SNAP.csv")
  v_yr <- c(2021, "genYr")
  y <- v_yr[1]
  pollutant <- c("nox")
  time_scale <- "yday"
  classification <- "SNAP"
  
  inv_year <- 2022

  v_sources <- tar_read_raw("v_sources")
  v_time_scales <- tar_read_raw("v_time_scales")
  v_ProfileIDs <- tar_read_raw("v_ProfileIDs")
  v_IDs <- v_ProfileIDs
  
  lookup_sectors <- "data/lookups/TEMREG_sectors.csv"
  data_sources <- "data/lookups/profile_id_names.csv"
  dt_NFR_PRO_AGG <- sector_to_ProID(lookup_sectors, data_sources)
  dt_lookup <- copy(dt_NFR_PRO_AGG)
  
  lookup_PID <- "./../../inventory_processor/data/lookups/pollutants.xlsx"
  dt_PID_NAEI <- as.data.table(read_excel(lookup_PID, sheet = "NAEI_pollutants"))
  dt_PID <- copy(dt_PID_NAEI)
  
  v_iso <- c(27)
  i <- 27
  
  v_fname <- tar_read_raw("v_fname_ProIDgams_yday")
  
  v_fname <- vec_c(tar_read_raw("v_fname_ProIDgams_yday"), tar_read_raw("v_fname_ProIDgams_HE"))
  
  
  
  v_fname <- vec_c(tar_read_raw("l_agg_GAMs_hour"), tar_read_raw("l_agg_GAMs_hourwday"),tar_read_raw("l_agg_GAMs_wday"),
        tar_read_raw("l_agg_GAMs_month"), tar_read_raw("l_agg_GAMs_yday"))
  
  
  tar_read_raw("v_fname_AggGAM_plots")
  
  
  l <- tar_read_raw("l_agg_GAMs_hour")
  
  length(l)
  lapply(l, length)
  
  
  v_fname <- vec_c(tar_read_raw("v_fname_ProIDgams_hour"), tar_read_raw("v_fname_ProIDgams_hourwday"),tar_read_raw("v_fname_ProIDgams_wday"),
        tar_read_raw("v_fname_ProIDgams_month"), tar_read_raw("v_fname_ProIDgams_yday"))


  
  fname <- v_fnames[12]
  
  dt_PID_NAEI <- tar_read_raw("dt_PID_NAEI")
  
  tar_read_raw("v_fname_ProIDgams_hour")
    
  ##############

}



