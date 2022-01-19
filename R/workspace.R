## initialise the workspace for splitting annual emissions into hourly profiles ##

## packages and workspace
packs <- c("sp","raster","stringr","gdalUtils","rgeos","rgdal","grid","plyr","ggplot2","ggrepel","data.table","stats","readr","ggplot2","sf","lubridate","units", "cowplot")
lapply(packs, require, character.only = TRUE)

"%!in%" <- Negate("%in%")

#spatial
BNG <<- suppressWarnings(CRS("+init=epsg:27700"))
LL <<- suppressWarnings(CRS("+init=epsg:4326"))
r_uk_BNG <<- raster(xmn=-230000, xmx = 750000, ymn = -50000, ymx = 1300000, res=1000, crs=BNG, vals=NA)

NAEI_pt_names <<- as.data.table(readxl::read_excel(paste0("//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/NAEI_data_and_SNAPS/NAEI_data/point/raw_data/NAEIPointsSources_2019.xlsx"), sheet="Data"))
NAEI_pt_names[,c("PollutantID","Datatype") := NULL]

## lookup table ##
dt_sect_to_prof <<- fread("./Data/Sectors.csv")

## profile tables ##
dt_prof_hour <<- melt(fread("./Data/hour_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year"), variable.name = "hour", value.name = "hour_coeff", variable.factor = F, value.factor = F ) %>% .[, hour := as.numeric(hour)]

dt_prof_hourwday <<- melt(fread("./Data/hourwday_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year","wday"), variable.name = "hour", value.name = "hourwday_coeff", variable.factor = F, value.factor = F ) %>% .[, hour := as.numeric(hour)]

dt_prof_wday <<- melt(fread("./Data/wday_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year"), variable.name = "wday", value.name = "wday_coeff", variable.factor = F, value.factor = F ) %>% .[, wday := as.numeric(wday)]

dt_prof_month <<- melt(fread("./Data/month_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year"), variable.name = "month", value.name = "month_coeff", variable.factor = F, value.factor = F ) %>% .[, month := as.numeric(month)]

dt_prof_yday <<- melt(fread("./Data/yday_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year"), variable.name = "yday", value.name = "yday_coeff", variable.factor = F, value.factor = F ) %>% .[, yday := as.numeric(yday)]

####################################################################################################
#### function to 1. return the NAEI emissions data, formatted, for the given year & Species.    ####
####             2. match the annual NAEI data to the NFR codes and the temporal profile codes. ####

JoinNAEItoProfiles <- function(year, species){
  
  ########################################################
  
  # year           = *numeric* year to process. Will determine calendar structure plus year specific profiles.
  if(!is.numeric(year)) stop ("Year is not numeric")
  # species        = *character* name of air pollutant or GHG or metal etc. Needs to conform to a list of options.
  if(species %!in% c("NOx","SOx","CH4","CO2","N2O") ) stop ("Species must be in: 
                                            AP:    BaP, CO, NH3, NMVOC, NOx, SO2
                                            PM:    PM2.5, PM10
                                            GHG:   CH4, CO2, N2O
                                            Metal: Cd, Cu, Hg, Ni, Pb, Zn")
  
  #########################################################
  
  
  colskeep <- c("NFR/CRF Group","Source","Activity",year)
  
  # read NAEI emissions data - currently on latest year = 2019. way to automate this?
  naei_files <- list.files("./data/NAEI_info", pattern="-2019.csv", full.names = T)
  dt_naei <- fread(naei_files[grep(tolower(species), naei_files)], header=T)
  
  # subset data (remove blank rows and superfluous columns)
  dt_naei <- dt_naei[Source != ""]
  dt_naei <- dt_naei[ ,..colskeep]
  
  # some formatting
  setnames(dt_naei, c("NFR/CRF Group",paste0(year)), c("NFR19","emission"))
  suppressWarnings(dt_naei[,emission := as.numeric(emission)])
  dt_naei <- dt_naei[!is.na(emission)]
  
  # set units
  if(species == "BaP"){
    units(dt_naei$emission) <- "Kg /yr"
  }else{
    units(dt_naei$emission) <- "Gg /yr"
  }
  
  # join the NAEI emissions to the temporal profile & classification data
  dt_joined <- dt_naei[dt_sect_to_prof, on = c("NFR19","Source","Activity")]
  dt_joined <- dt_joined[!is.na(emission)]
  
  if(identical(sum(dt_joined$emission, na.rm=T), sum(dt_naei$emission, na.rm=T))==T){
    NULL
  }else{
    print("Emissions total has changed since joining to profile table - CHECK")
    break
  }
    
  return(dt_joined)
  
}

##########################################################################################################
#### function to profile emissions into classification system, using the profile IDs per NFR sector.  ####
####  Result: a weighted mean temporal profile for the whole Sector, based on emissions per profile   ####

TempProfileBySector <- function(year, species, timestep, classification = c("SNAP","GNFR"), emis, yr_spec_NFR = NULL){
  
  ####################################################
  
  # year           = *numeric* year to process. Will determine calendar structure plus year specific profiles.
  if(!is.numeric(year)) stop ("Year is not numeric")
  # species        = *character* name of air pollutant or GHG or metal etc. Needs to conform to a list of options.
  if(species %!in% c("NOx","SOx","CH4","CO2","N2O") ) stop ("Species must be in: 
                                            AP:    BaP, CO, NH3, NMVOC, NOx, SO2
                                            PM:    PM2.5, PM10
                                            GHG:   CH4, CO2, N2O
                                            Metal: Cd, Cu, Hg, Ni, Pb, Zn")
  # timestep       = *vector* timesteps to recreate profiles for
  if(timestep %!in% c("yday","month","wday","hour","hourwday")) stop("Not currently a valid timestep to run profiles for")
  # classification = *character* GNFR or SNAP
  classification <- match.arg(classification)
  # emis           = *data.table* from JoinNAEItoProfiles; total emissions per NFR with matching profile IDs
  # yr_spec_NFR    = *character* optional vector of NFR codes to be profiled by year specific data
  if(sum(yr_spec_NFR %!in% unique(dt_sect_to_prof$NFR19))) stop(paste0("The following are not NFR codes, check: ",yr_spec_NFR[yr_spec_NFR %!in% unique(dt_sect_to_prof$NFR19)]))
  
  ####################################################
  
  ## loop through the timesteps and create new coefficients, for aggregated Sector level
  ## (using a calendar method smooths over the individual temporal profiles, it's not right)
  
  print(paste0(Sys.time(),": Profiling ",classification," by ",timestep))
  
  # multiplier to adjust the coefficents to a fraction of 1. 
  c_mult <- ifelse(timestep == "yday", 1/365, ifelse(timestep == "month", 1/12, ifelse(timestep == "wday", 1/7, ifelse(timestep == "hour", 1/24, 1/7/24))))
  
  # set sectors to work through
  if(classification == "SNAP"){
    sectors <- 1:11 
  }else{
    sectors <- c("ff") # FILL IN GNFR HERE
  }
  
  # blank list for sector level profiles
  l_sec_by_time <- list()
  
  # loop through sectors set above and do temporal processing
  for(s in sectors){
    
  # extract the sector relevant NFR codes from NAEI data. Aggregate by the profile. 
  dt_sect_specific <- emis[SNAP == s]
  dt_sect_specific_agg <- dt_sect_specific[, .(emission = sum(emission, na.rm=T)), by = .(Sector = get(classification), Profile_ID)]
  
  # sum of sector emissions (for checking at end)
  total_emis_checker <- sum(dt_sect_specific$emission)
  units(total_emis_checker) <- "Gg /yr"
  
  l_prof_profiles <- list()
      
  ## loop through the unique temporal profiles in sector and disperse emissions over timestep.
  # timesteps are according to the Profile name and year of choice. 
  for(prof in unique(dt_sect_specific_agg[,Profile_ID])){
    
    if(prof %in% yr_spec_NFR){
      dt_p <- get(paste0("dt_prof_",timestep))[Profile_ID == prof & Year == year]    
    }else{
      dt_p <- get(paste0("dt_prof_",timestep))[Profile_ID == prof & Year == 0]    
    }
      
    # back to a fraction of 1
    dt_p[, coeff := get(paste0(timestep,"_coeff")) * c_mult]
    
    # total emissions for the sector and profile ID
    Ekt <- dt_sect_specific_agg[Profile_ID == prof, emission]
    dt_p[, ann_emis_kt := Ekt]
      
    dt_p[,frac_emis_kt := ann_emis_kt * coeff]
      
    l_prof_profiles[[paste0(prof,"_",timestep)]] <- dt_p
      
  } # end of profile ID loop
    
    # combine the data, emissions now distributed over time, at Profile_ID level
    dt_prof_profiles <- rbindlist(l_prof_profiles, use.names = T)
    
    # sum all the emissions for a sector total over the timesteps
    if(timestep=="hourwday"){
      dt_p_emis <- dt_prof_profiles[,.(timestep_emis_kt = sum(frac_emis_kt)), by=.(wday,hour)]
    }else{
      dt_p_emis <- dt_prof_profiles[,.(timestep_emis_kt = sum(frac_emis_kt)), by=get(timestep)]
      setnames(dt_p_emis, "get", timestep) 
    }
    
    # check emissions totals still match
    if((as.numeric(sum(dt_p_emis$timestep_emis_kt) / total_emis_checker)) < 0.999 | 
       (as.numeric(sum(dt_p_emis$timestep_emis_kt) / total_emis_checker)) > 1.001){
      print(paste0("Sector ",sector," hourly emissions do not add up to NAEI sector total. Check yday." ))
    }else{
      NULL
    }
    
    # add in columns etc., ready to list. Centre the coefficient around 1 again
    dt_p_emis[, c("Sector", "Pollutant", "coeff") := list(s, species, (timestep_emis_kt/sum(timestep_emis_kt)) / c_mult)]
    dt_p_emis[, coeff := set_units(coeff, NULL)]
    
    # add to master list
    l_sec_by_time[[as.character(s)]] <- dt_p_emis
      
 } # end of sector loop
  
  # combine, write & return
  dt_sec_by_time <- rbindlist(l_sec_by_time, use.names = T)
  
  if(timestep=="hourwday"){
    cols <- c("Pollutant", "Sector", "wday", "hour", "coeff")
  }else{
    cols <- c("Pollutant", "Sector", timestep, "coeff")
  }
  
  dt_return <- dt_sec_by_time[,..cols]
  
  fwrite(dt_return,paste0("./doc/coefficients_sector/",species,"_",timestep,"_profiles_",classification,"_",year,"_",v_num,".csv"))
  
  return(dt_return)
  
} # end of function


#########################################################################################################
#### function to visualise the above sectoral level profiles, as emissions graphs. (Sectors & total) ####

PlotEmissionsOverTime <- function(year, species, classification, sec_profs, emis){
  
  y_spec_NFR <- NULL
  
  if(classification=="SNAP"){
    emis <- emis[!(get(classification) %in% c("","avi.cruise")) ] # need to get rid of aviation cruise
  }else{
    NULL  # might be some things to put here in future
  }
  
  # aggregate the emissions to sector level, to distribute onto sector profiles. 
  dt_sect_emis <- emis[, .(emission = sum(emission, na.rm=T)), by = .(Sector = get(classification))]
  if(classification=="SNAP") dt_sect_emis[,Sector := as.numeric(as.character(Sector))]
  
  # aggregate the emissions to sector level, to distribute onto sector profiles. 
  dt_tot_emis <- emis[, .(emission = sum(emission, na.rm=T))]
  
  l_plot <- list()
  l_ymwhhw <- list()
  
  # cycle through individual timesteps and create graphs of emissions, per sector. 
  for(p in c("yday","month","wday","hour","hourwday")){
    
    # select multipler to make fraction of 1
    c_mult <- ifelse(p == "yday", 1/365, ifelse(p == "month", 1/12, ifelse(p == "wday", 1/7, ifelse(p == "hour", 1/24, 1/7/24))))
    
    # select correct coefficient table and join on the sector emission total
    dt_p <- sec_profs[[p]]
    dt_p <- dt_sect_emis[dt_p,on="Sector"]
    dt_p[, emission := set_units(emission, NULL)]
    dt_p[, emis_adj := emission * (coeff * c_mult)]
    dt_p[, Sector := factor(Sector)]
    dt_p[,Time_unit := p]
    
    #if(p=="hourwday"){
    #  setnames(dt_p, "hour", "Time_value")
    #}else{
    #  setnames(dt_p, paste0(p), "Time_value")
    #}
    
    
    # check emissions totals still match
    if((as.numeric(sum(dt_p$emis_adj) / dt_tot_emis)) < 0.999 | 
       (as.numeric(sum(dt_p$emis_adj) / dt_tot_emis)) > 1.001){
      print(paste0("Sector ",sector," hourly emissions do not add up to NAEI sector total. Check yday." ))
    }else{
      NULL
    }
    
    l_ymwhhw[[p]] <- dt_p
    
    if(p %in% c("month","hour")){
     
      g1 <- ggplot()+
        geom_line(data=dt_p, aes(x=get(p) , y=emis_adj, group=Sector, colour=Sector))+
        labs(x=p,y="Emission (kt)")+
        #facet_wrap(~Sector, nrow=1)+
        theme_bw()+
                theme(axis.text.x = element_text(size=12),
                axis.text.y = element_text(size=12),
                axis.title = element_text(size=16),
                legend.title=element_blank(),
                legend.position = "none")
      
      l_plot[[p]] <- ggplotGrob(g1)
      
    }else if(p %in% c("yday","wday")){
      
      g1 <- ggplot()+
        geom_line(data=dt_p, aes(x=get(p) , y=emis_adj, group=Sector, colour=Sector))+
        labs(x=p,y=NULL)+
        #facet_wrap(~Sector, nrow=1)+
        theme_bw()+
        theme(axis.text.x = element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.title = element_text(size=16),
              legend.title=element_blank(),
              legend.text= element_text(size=12),
              legend.position = NULL)
      
      l_plot[[p]] <- ggplotGrob(g1)
    
    }else{
      
      g1 <- ggplot()+
        geom_line(data=dt_p, aes(x=hour, y=emis_adj, group=Sector, colour=Sector))+
        labs(x="hour (by wday)",y="Emission (kt)")+
        facet_wrap(~wday, nrow=1)+
        theme_bw()+
        theme(axis.text.x = element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.title = element_text(size=16),
              legend.text= element_text(size=12),
              legend.title=element_blank())
      
      l_plot[[p]] <- ggplotGrob(g1)
      
    } # if else hourwday
  
  } # end of time loop
  
  # total emissions plot - do this on a yday and hourwday basis, adjusting the yday for wday
  ## create a blank calendar for the given year to attach profiled emissions to. 
  #startDate <- dmy_hm(paste0("01/01/",year," 00:00"))
  #endDate <- dmy_hm(paste0("31/12/",year," 23:00"))
  #
  #dt_calendar <- data.table(DateTime = seq(from = startDate, to = endDate, by = 3600))
  #dt_calendar[,n_mon_days := days_in_month(DateTime)]
  #
  #dt_calendar[ , c("yday","month", "wday", "hour") := list(yday(DateTime), month(DateTime), lubridate::wday(DateTime, week_start = getOption("lubridate.week.start", 1)), hour(DateTime))]
  #
  # sum the emissions to temporal steps for total emissions. Hourwday must remain separate. 
  #l1 <- lapply(l_ymwhhw[1:4], function(x) x[, .(p_emis = sum(emis_adj)),by=get(names(x)[4])])
  #sapply(1:4, function(x) setnames(l1[[x]], "get", names(l1)[x]))
  #sapply(1:4, function(x) setnames(l1[[x]], "p_emis", paste0(names(l1)[x],"_emis")))
  #l2 <- l_ymwhhw[[5]][,.(p_emis = sum(emis_adj)),by=.(wday,hour)]
  #setnames(l2,"p_emis","hourwday_emis")
  #
  # join emissions by temporal fraction to the calendar
  #dt_calendar_emis <- dt_calendar[l1[["yday"]], on = "yday"][l1[["month"]], on = "month"][l1[["wday"]], on = "wday"][l2, on = c("wday","hour")]
  
  # adjust the yday to take into account the actual day of the week (differs year to year)
  #dt_calendar_emis[,yday_emis_adj := yday_emis]
  
  ## save data
  d1 <- plot_grid(l_plot[["hour"]], l_plot[["wday"]], l_plot[["month"]], l_plot[["yday"]], ncol=2, nrow=2,rel_heights = c(1,1,1,1), rel_widths = c(1,1,1,1), align="h", axis="l")
  d2 <- plot_grid(d1,l_plot[["hourwday"]], ncol=1, nrow=2, rel_heights = c(1.1,0.9), align="h", axis="l")
    
  save_plot(paste0("./doc/coefficients_sector/",species,"_",classification,"_",y_emis,"_EmisPlots_",v_num,".png"), d2, base_height = 14, base_width = 16)
  
  
  ## metadata
  v_file <- paste0("./doc/coefficients_sector/",species,"_",classification,"_",y_emis,"_",v_num,".info")
  file.create(v_file,overwrite=T)
  
  line1 <- paste0("Sector level coefficients created for version ",v_num)
  line2 <- paste0("Species: ", species)
  line3 <- paste0("Emissions Year: ", year)
  line4 <- paste0("Classification system: ", classification)
  line5 <- ifelse(is.null(y_spec_NFR),"Year specific profiles used for: No NFR sectors", paste0("Year specific profiles used for: ", y_spec_NFR))
  line6 <- paste0("Coefficients created for: ", paste(names(sec_profs), collapse = " "))
  
  lines = paste(line1,line2,line3,line4,line5,line6,
             sep="\n")
  
  write(lines,file=v_file,append=TRUE)
  
}

###############################################################################################
###############################################################################################

#### Sentinel 5p comparisons ####

######################################################################################################
#### function to process NAEI UK emissions data to AOI ####
AOIEmissionsData <- function(aoi_name, species, year){
  
  ### extent of London study area + an expanded area for plotting interest (to ID large nearby sources?)
  aoi <- suppressWarnings(suppressMessages(st_read("./AOI_extents",paste0("aoi_",aoi_name), quiet = T)))
  aoi_bng <- st_transform(aoi, BNG)
  
  aoi_expand <- raster::extent(extent(aoi_bng)[1] - (extent(aoi_bng)[2] - extent(aoi_bng)[1])/2,
                               extent(aoi_bng)[2] + (extent(aoi_bng)[2] - extent(aoi_bng)[1])/2,
                               extent(aoi_bng)[3] - (extent(aoi_bng)[4] - extent(aoi_bng)[3])/2,
                               extent(aoi_bng)[4] + (extent(aoi_bng)[4] - extent(aoi_bng)[3])/2)
  
  ### for each SNAP, read in year/spec and crop/stack
  NAEI_diff_files <- list.files(paste0("//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/NAEI_data_and_SNAPS/Emissions_grids_plain/BNG/nox/diffuse/",year,"/rasters_SNAP"), pattern = paste0("^",species,"_diff_",year,"_uk_SNAP.*.t_1km_BNG_2019NAEImap.tif$"), full.names = T)
  
  NAEI_pt_files <- fread(paste0("//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/NAEI_data_and_SNAPS/Emissions_grids_plain/BNG/nox/point/",year,"/",species,"_pt_",year,"_uk_SNAP_t_BNG.csv"))
  
  
  
  l_SNAP_tots <- list()
  
  st_dif   <- stack()
  st_difE  <- stack()
  st_pt    <- stack()
  st_ptE   <- stack()
  st_SNAP  <- stack()
  st_SNAPE <- stack()
  
  for(SN in paste0("S",1:11)){
    
    print(paste0(SN))
    
    r_file <- NAEI_diff_files[grep(paste0("_",SN,"_"), NAEI_diff_files)]
    r <- suppressWarnings(raster(r_file))
    
    pts <- NAEI_pt_files[SNAP==as.numeric(substr(SN,2,nchar(SN)))]
    units(pts$Emission) <- "Mg/yr"
    
    if(nrow(pts)>0){
      
      pt_xy <- copy(pts)
      coordinates(pt_xy) <- ~Easting + Northing
      crs(pt_xy) <- BNG
      
      r_pt <- rasterize(pt_xy, r_uk_BNG, field="Emission", fun=`sum`)
      
    }else{
      r_pt <- r_uk_BNG
    }
    
    # crop to AOI and expanded AOI
    rc_dif <- crop(r, aoi_bng)
    rc_difE <- crop(r, aoi_expand)
    
    rc_pt <- crop(r_pt, aoi_bng)
    rc_ptE <- crop(r_pt, aoi_expand)
    
    rc_SNtot <- suppressWarnings(calc(stack(rc_dif, rc_pt), sum, na.rm=T))
    names(rc_SNtot) <- paste0(species,"_",SN)
    rc_SNtotE <- suppressWarnings(calc(stack(rc_difE, rc_ptE), sum, na.rm=T))
    names(rc_SNtotE) <- paste0(species,"_",SN)
    
    # stacks
    st_dif  <- stack(st_dif, rc_dif)
    st_difE  <- stack(st_difE, rc_difE)
    
    st_pt   <- stack(st_pt, rc_pt)
    st_ptE   <- stack(st_ptE, rc_ptE)
    
    st_SNAP <- stack(st_SNAP, rc_SNtot)
    st_SNAPE <- stack(st_SNAPE, rc_SNtotE)
    
    # totals
    tot_SNdif <- cellStats(rc_dif, sum)
    tot_SNdifE <- cellStats(rc_difE, sum)
    tot_SNpt  <- cellStats(rc_pt, sum)
    tot_SNptE  <- cellStats(rc_ptE, sum)
    tot_SNAP  <- cellStats(rc_SNtot, sum)
    tot_SNAPE  <- cellStats(rc_SNtotE, sum)
    
    # stats table
    dt_row <- data.table(AOI = aoi_name, Year=year, Pollutant=species, SNAP = substr(SN, 2, nchar(SN)), emission_diff = tot_SNdif, emission_diff_exp = tot_SNdifE, emission_pt = tot_SNpt, emission_pt_exp = tot_SNptE, emission_total = tot_SNAP, emission_total_exp = tot_SNAPE)
    
    l_SNAP_tots[[SN]] <- dt_row
    
  }
  
  
  # bind the SNAP AOI totals and make a sum raster (tonnes)
  dt_SNAP_tots <- rbindlist(l_SNAP_tots, use.names = T)
  dt_SNAP_tots[, SNAP := as.numeric(SNAP)]
  
  # biggest point sources
  top5_pts <- NAEI_pt_files[Easting <= aoi_expand[2] & Easting >= aoi_expand[1] & Northing <= aoi_expand[4] & Northing >= aoi_expand[3]][order(-Emission)][1:5]
  
  top5_pts <- NAEI_pt_names[Pollutant=="NOx"][top5_pts, on=c("Easting","Northing")]
  
  # total snaps
  r_diffuse <- calc(st_dif, sum, na.rm=T)
  r_diffuseE <- calc(st_difE, sum, na.rm=T)
  r_point <- calc(st_pt, sum, na.rm=T)
  r_pointE <- calc(st_ptE, sum, na.rm=T)
  r_total <- calc(st_SNAP, sum, na.rm=T)
  r_totalE <- calc(st_SNAPE, sum, na.rm=T)
  
  st_aoi <- stack(r_diffuse, r_point, r_total)
  names(st_aoi) <- c("diffuse","point","total")
  st_exp <- stack(r_diffuseE, r_pointE, r_totalE)
  names(st_exp) <- c("diffuse","point","total")
  
  return(list("dt_aoi_emis" = dt_SNAP_tots, "st_aoi" = st_aoi, "st_expanded" = st_exp, "top5_pts" = top5_pts))
  
}


###########################################################################################
#### function to process AOI emissions data into one total temporal profile (SNAP weighted) ####

AOITemporalProfile <- function(aoi_emis, species, year, write){
  
  units(aoi_emis$emission_total) <- "Mg/yr"
  units(aoi_emis$emission_total_exp) <- "Mg/yr"
  
  ## bring in profile data
  dt_spec_SNAP <- readRDS(paste0("./",paste0(toupper(substr(species,1,2)),substr(species,3,3)),"_SNAP_",year,"_tp.rds"))
  dt_spec_hourly <- rbindlist(lapply(dt_spec_SNAP,`[[`,4), use.names = T)
  dt_spec_hourly[, doy := yday(DateTime)]
  
  dt_spec_daily <- dt_spec_hourly[,.(sector_doy_emission = sum(sector_emission)), by=.(doy, Sector)]
  dt_spec_daily[,doy_coeff := sector_doy_emission/sum(sector_doy_emission, na.rm=T), by=Sector]
  units(dt_spec_daily$doy_coeff) <- NULL
  setnames(dt_spec_daily, "Sector", "SNAP")
  
  ## join the SNAP totals from NAEI map to the profiles
  dt_spec_daily <- aoi_emis[dt_spec_daily, on = "SNAP"]
  
  dt_spec_daily[, aoi_doy_emission := emission_total * doy_coeff]
  dt_spec_daily[, exp_doy_emission := emission_total_exp * doy_coeff]
  
  ## sum the aoi emissions back to one total. split this out to a new total temporal profile (aoi & exp).
  dt_aoi_tp_sum <- dt_spec_daily[,.(doy_emission = sum(aoi_doy_emission)), by=doy]
  dt_aoi_tp_sum[, doy_coeff := doy_emission/sum(doy_emission)]
  units(dt_aoi_tp_sum$doy_coeff) <- NULL
  
  dt_exp_tp_sum <- dt_spec_daily[,.(doy_emission = sum(exp_doy_emission)), by=doy]
  dt_exp_tp_sum[, doy_coeff := doy_emission/sum(doy_emission)]
  units(dt_exp_tp_sum$doy_coeff) <- NULL
  
  ## also make temporal profiles for the SNAPs.
  dt_aoi_tp_SN <- dt_spec_daily[,.(doy_emission = sum(aoi_doy_emission)), by=.(doy, SNAP)]
  dt_aoi_tp_SN[, doy_coeff := doy_emission/sum(doy_emission), by=SNAP]
  units(dt_aoi_tp_SN$doy_coeff) <- NULL
  
  dt_exp_tp_SN <- dt_spec_daily[,.(doy_emission = sum(exp_doy_emission)), by=.(doy, SNAP)]
  dt_exp_tp_SN[, doy_coeff := doy_emission/sum(doy_emission), by=SNAP]
  units(dt_exp_tp_SN$doy_coeff) <- NULL
  
  # check #
  identical(sum(aoi_emis$emission_total), sum(dt_spec_daily$aoi_doy_emission), sum(dt_aoi_tp_sum$doy_emission), sum(dt_aoi_tp_SN$doy_emission))
  identical(sum(aoi_emis$emission_total_exp), sum(dt_spec_daily$exp_doy_emission), sum(dt_exp_tp_sum$doy_emission), sum(dt_exp_tp_SN$doy_emission))
  
  # to return
  rt <- list(dt_aoi_tp_sum, dt_exp_tp_sum, dt_aoi_tp_SN, dt_exp_tp_SN)
  names(rt) <- c("DOY_emis_tp_AOI","DOY_emis_tp_EXP", "DOY_SN_emis_tp_AOI", "DOY_SN_emis_tp_EXP")
  
  ## write if needed
  if(write==T){
    
    fwrite(rt[["DOY_emis_tp_AOI"]], paste0("./",aoi_name,"/",aoi_name,"_AOI_",species,"_emis_tp_doy_",y_emis,".csv"))
    fwrite(rt[["DOY_emis_tp_EXP"]], paste0("./",aoi_name,"/",aoi_name,"_EXP_",species,"_emis_tp_doy_",y_emis,".csv"))
    fwrite(aoi_emis, paste0("./",aoi_name,"/",aoi_name,"_",species,"_emissions_t_",y_emis,".csv"))
    
  }else{
    print("data files not written")
  }
  
  return(rt)
  
  
}


###########################################################################################
#### function to plot the profile and emissions data

PlottingAOIdata <- function(species, year, aoi_name, r_aoi_tot, r_exp_tot, aoi_snap_tot, tps, pt_labels){
  
  tps[["DOY_SN_emis_tp_AOI"]][,SNAP := factor(SNAP, levels=1:11)]
  tps[["DOY_SN_emis_tp_EXP"]][,SNAP := factor(SNAP, levels=1:11)]
  
  aoi_snap_tot[, aoi_diff_share := sum(emission_diff) / sum(emission_total)]
  aoi_snap_tot[, aoi_pt_share := sum(emission_pt) / sum(emission_total)]
  
  aoi_snap_tot[, aoi_dif_SNAP_share := (emission_diff/sum(emission_diff)) * aoi_diff_share]
  aoi_snap_tot[, aoi_pt_SNAP_share := (emission_pt/sum(emission_pt)) * aoi_pt_share]
  
  aoi_snap_tot[, aoi_tot_SNAP_share := emission_total/sum(emission_total)]
  aoi_snap_tot[, aoi_tot_SNAP_share_rank := frank(aoi_tot_SNAP_share)]
  aoi_snap_tot[, exp_tot_SNAP_share := emission_total_exp/sum(emission_total_exp)]
  
  aoi <- suppressWarnings(suppressMessages(st_read("./AOI_extents",paste0("aoi_",aoi_name), quiet = T)))
  aoi_bng <- st_transform(aoi, BNG)
  
  ## plot temporal profile for total emissions in the AOI & in the EXP
  gg <- ggplot(data=tps[["DOY_emis_tp_AOI"]], aes(x=doy,y=doy_coeff))+
    geom_line()+
    geom_smooth()+
    stat_smooth(data=tps[["DOY_SN_emis_tp_AOI"]][SNAP %in% aoi_snap_tot[aoi_tot_SNAP_share_rank >= 7, SNAP]], aes(group=SNAP, colour=SNAP),geom='line', se=F, alpha=0.6)+
    geom_hline(yintercept = 1/365, linetype="dashed")+
    #geom_smooth(data=tps[["DOY_emis_tp_EXP"]], linetype="dashed", se=F)+
    #geom_text_repel(data=tps[["DOY_SN_emis_tp_AOI"]][doy==365 & SNAP %in% aoi_snap_tot[aoi_SNAP_share_rank >= 7, SNAP]], aes(label=SNAP), hjust=0, nudge_x = 2, size=3.3, colour="black" , segment.color = "grey70" )+
    scale_x_continuous(name="Day of Year", breaks=c(0,50,100,150,200,250,300,350))+
    scale_y_continuous(name = "Coefficient", breaks=c(0.001,0.002,0.003,0.004, 0.005), limits=c(0, 0.005))+
    theme_bw()+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=16),
          legend.text = element_text(hjust = 0, size = 14),axis.ticks = element_blank(),
          legend.title = element_text(size = 16, hjust = 1))
  
  
  #ggsave(paste0("London_",species,"_tp_doy_total_",y_emis,".png"), gg, width=14, height=8)
  
  ggSNprop <- ggplot()+
    #geom_bar(data=aoi_snap_tot, aes(x=SNAP,y=aoi_tot_SNAP_share),stat="identity", fill="#F07F69")+
    geom_bar(data = melt(aoi_snap_tot[,c("SNAP","aoi_dif_SNAP_share","aoi_pt_SNAP_share")], id.vars = "SNAP"), aes(x=SNAP, y=value, group=variable,fill=variable), position="stack", stat = "identity")+
    geom_bar(data=aoi_snap_tot, aes(x=SNAP,y=exp_tot_SNAP_share),stat="identity", fill=NA, colour="black", linetype="dashed")+
    scale_fill_manual(values=c("#F07F69","#6FE5F1"))+
    scale_x_continuous(breaks=c(1:11))+
    scale_y_continuous(name = "Proportion")+
    annotate(geom="text", x=1,y=0.5, label = "Red    = AOI (diffuse)", size=3.5, hjust=0)+
    annotate(geom="text", x=1,y=0.47, label = "Blue   = AOI (points)", size=3.5, hjust=0)+
    annotate(geom="text", x=1,y=0.44, label = "Black = Expanded", size=3.5, hjust=0)+
    ggtitle("Proportion emissions by SNAP")+
    theme_bw()+
    theme(legend.position = "none",
          axis.text = element_text(size=10),
          axis.title = element_text(size=11))
  
  
  #gginset <- ggdraw(gg)+
  #  draw_plot(ggSNprop, 0.075,0.11,0.35,0.31, hjust = 0, vjust=0)
  
  #ggsave(paste0("London_",species,"_tp_doy_total_withinset_",y_emis,".png"), gginset, width=14, height=8)
  
  dt_exp_toplot <- as.data.table(as.data.frame(r_exp_tot, xy=T))
  dt_exp_dif <- dt_exp_toplot[,c("x","y","diffuse")]
  dt_exp_dif[diffuse<0.001, diffuse:=NA]
  dt_exp_pt <- dt_exp_toplot[,c("x","y","point")]
  dt_exp_pt[point==0, point:=NA]
  
  midp_dif <- (min(dt_exp_dif$diffuse,na.rm=T) + max(dt_exp_dif$diffuse,na.rm=T))/2
  midp_pt <- (min(dt_exp_pt$point,na.rm=T) + max(dt_exp_pt$point,na.rm=T))/2
  
  gmap_dif <- ggplot() +
    geom_raster(data = dt_exp_dif , aes(x = x, y = y, fill = diffuse))+
    # coord_fixed()+
    geom_rect(aes(xmin = extent(aoi_bng)[1], xmax = extent(aoi_bng)[2],   ymin = extent(aoi_bng)[3], ymax = extent(aoi_bng)[4]),   fill = NA, colour="black") +
    #geom_sf(data = aoi_bng, fill = NA, colour = "black", size = 0.3)+
    #coord_sf(crs=st_crs(27700))+
    scale_fill_gradient2(name=bquote(tonnes), midpoint = midp_dif, low = "#fde0dd", mid = "#fa9fb5",high = "#c51b8a",na.value = "white")+
    ggtitle("Diffuse emissions")+
    guides(fill = guide_colourbar(barwidth = 1.5, barheight = 10))+
    theme(axis.line = element_blank(),
          axis.text = element_text(size=10),
          axis.title = element_blank(),
          panel.grid.major = element_line(colour = "transparent"),
          panel.background = element_rect(fill = NA, colour = NA),
          #strip.text.x = element_blank(),
          #strip.background = element_blank(),
          legend.text = element_text(hjust = 0, size = 12),axis.ticks = element_blank(),
          legend.title = element_text(size = 14, hjust = 0))
  
  gmap_pt <-  ggplot() +
    geom_raster(data = dt_exp_pt , aes(x = x, y = y, fill = point))+
    #coord_fixed()+
    geom_rect(aes(xmin = extent(aoi_bng)[1], xmax = extent(aoi_bng)[2],   ymin = extent(aoi_bng)[3], ymax = extent(aoi_bng)[4]),   fill = NA, colour="black") +
    #geom_sf(data = aoi_bng, fill = NA, colour = "black", size = 0.3)+
    #coord_sf(crs=st_crs(27700))+
    scale_fill_gradient2(name=bquote(tonnes), midpoint = midp_pt, low = "#fde0dd", mid = "#fa9fb5",high = "#c51b8a",na.value = "white")+
    geom_text_repel(data=pt_labels, aes(x=Easting,y= Northing, label=paste0(Site," (",SectorID,", ",SNAP,")")), hjust=0, size=2, colour="black", min.segment.length = 0)+
    ggtitle("Point emissions")+
    guides(fill = guide_colourbar(barwidth = 1.5, barheight = 10))+
    theme(axis.line = element_blank(),
          axis.text = element_text(size=10),
          axis.title = element_blank(),
          panel.grid.major = element_line(colour = "transparent"),
          panel.background = element_rect(fill = NA, colour = NA),
          #strip.text.x = element_blank(),
          #strip.background = element_blank(),
          legend.text = element_text(hjust = 0, size = 12),axis.ticks = element_blank(),
          legend.title = element_text(size = 14, hjust = 0))
  
  
  
  d1 <- plot_grid(ggSNprop, gmap_dif, gmap_pt,labels=c("B","C","D"), ncol=3, nrow=1, rel_widths = 1, rel_heights = 1, align="h", axis="l")
  d2 <- plot_grid(gg,NULL, d1,labels=c("A","",""), ncol=1,rel_heights = c(1.8,-.02,1), align="h", axis="l")
  
  save_plot(paste0("./",aoi_name,"/",aoi_name,"_",species,"_plots_",year,".png"), d2, base_height = 15, base_width = 15)
  
  print(paste0("plot for ",aoi_name," complete and saved"))
  
  
}



###############################################################################################
###############################################################################################
### BELOW IS CALENDAR METHOD EMISSIONS PROFILE - DONT DELETE! ###

## create a blank calendar for the given year to attach profiled emissions to. 
#startDate <- dmy_hm(paste0("01/01/",year," 00:00"))
#endDate <- dmy_hm(paste0("31/12/",year," 23:00"))

#dt_calendar <- data.table(DateTime = seq(from = startDate, to = endDate, by = 3600))
#dt_calendar[,n_mon_days := days_in_month(DateTime)]

#dt_calendar[ , c("yday","month", "wday", "hour") := list(yday(DateTime), month(DateTime), lubridate::wday(DateTime, week_start = getOption("lubridate.week.start", 1)), hour(DateTime))]


## empty list for the hourly emissions for the profile, per SNAP/GNFR (might be > 1 sector in a profile ID)
#l_sector_profiles_yday <- list()
#l_sector_profiles_month <- list()

## loop through the unique temporal profiles and disaggregate emissions to hourly
#for(prof in unique(dt_sect_specific_agg[,Profile_ID])){

## all data has access to yday/month/wday/hour/hourwday. 
## using yday, month is not needed. wday is needed to know which hourwday to use. 

### method 'yday': using yday generated GAM and attaching hour, determined by wday ###
### the month is not used here. The yday is used and the month cab be calculated ###
# subset the temporal profile tables to the specific profile ID in the loop.
#if(prof %in% yr_spec_NFR){
#  dt_yday_ID <- dt_prof_yday[Profile_ID == prof & Year == year]
#  dt_hourwday_ID <- dt_prof_hourwday[Profile_ID == prof & Year == year]
#}else{
#  dt_yday_ID <- dt_prof_yday[Profile_ID == prof & Year == 0]
#  dt_hourwday_ID <- dt_prof_hourwday[Profile_ID == prof & Year == 0]
#}
#dt_yday_ID[,c("Pollutant","Year","Profile_ID") := NULL]
#dt_hourwday_ID[,c("Pollutant","Year","Profile_ID") := NULL]

# join the profiles to the hourly calendar
#dt_hours_coeffs <- dt_calendar[dt_yday_ID, on = "yday"][dt_hourwday_ID, on = c("wday","hour")]

# total emissions for the sector and profile ID
#Ekt <- dt_sect_specific_agg[Sector == sector & Profile_ID == prof, emission]
#dt_hours_coeffs[, ann_emis_kt := Ekt]
#  
#if(species == "BaP"){
#  units(dt_hours_coeffs$ann_emis_kt) <- "Kg /yr"
#}else{
#  units(dt_hours_coeffs$ann_emis_kt) <- "Gg /yr"
#}

# split emissions out into hours and day of the year
# Readjust to total (sometimes very slightly out due to proportion of hod/dow in year)
#dt_hours_coeffs[,     hour_emis := ann_emis_kt * (yday_coeff/365) * (hrwd_coeff/24)]
#dt_hours_coeffs[, hour_emis_adj := ( ((yday_coeff/365) * ann_emis_kt) / sum(hour_emis, na.rm=T)) * hour_emis, by = yday]

# add some information
#dt_hours_coeffs[,Sector := sector]
#dt_hours_coeffs[,Profile_ID := prof]
#setkey(dt_hours_coeffs, DateTime)
#dt_hours_coeffs[,Method := "yday"]

#l_sector_profiles_yday[[paste0("Sector_",sector," Profile_",prof)]] <- dt_hours_coeffs

### method 'month': reflecting EMEP, using month/wday/hourwday ###
# subset the temporal profile tables to the specific profile ID in the loop.
#if(prof %in% yr_spec_NFR){
#  dt_month_ID <- dt_prof_month[Profile_ID == prof & Year == year]
#  dt_wday_ID <- dt_prof_wday[Profile_ID == prof & Year == year]
#  dt_hourwday_ID <- dt_prof_hourwday[Profile_ID == prof & Year == year]
#}else{
#  dt_month_ID <- dt_prof_month[Profile_ID == prof & Year == 0]
#  dt_wday_ID <- dt_prof_wday[Profile_ID == prof & Year == 0]
#  dt_hourwday_ID <- dt_prof_hourwday[Profile_ID == prof & Year == 0]
#}

#dt_month_ID[,c("Pollutant","Year","Profile_ID") := NULL]
#dt_wday_ID[,c("Pollutant","Year","Profile_ID") := NULL]
#dt_hourwday_ID[,c("Pollutant","Year","Profile_ID") := NULL]

# join the profiles to the hourly calendar
#dt_hours_coeffs <- dt_calendar[dt_month_ID, on = "month"][dt_wday_ID, on = "wday"][dt_hourwday_ID, on = c("wday","hour")]

# total emissions for the sector and profile ID
#Ekt <- dt_sect_specific_agg[Sector == sector & Profile_ID == prof, emission]
#dt_hours_coeffs[, ann_emis_kt := Ekt]

#if(species == "BaP"){
#  units(dt_hours_coeffs$ann_emis_kt) <- "Kg /yr"
#}else{
#  units(dt_hours_coeffs$ann_emis_kt) <- "Gg /yr"
#}

# split emissions out into hours, wday, month.
# Readjust to total (very slightly out due to proportion of days not in full week in a month)
#dt_hours_coeffs[,     hour_emis := ann_emis_kt * (month_coeff/12) * (7/n_mon_days) * (wday_coeff/7) * (hrwd_coeff/24)]
#dt_hours_coeffs[, hour_emis_adj := ( ((month_coeff/12) * ann_emis_kt) / sum(hour_emis, na.rm=T)) * hour_emis, by = month]

# add some information
#dt_hours_coeffs[,Sector := sector]
#dt_hours_coeffs[,Profile_ID := prof]
#setkey(dt_hours_coeffs, DateTime)
#dt_hours_coeffs[,Method := "month"]

#l_sector_profiles_month[[paste0("Sector_",sector," Profile_",prof)]] <- dt_hours_coeffs

#} # end of unique profile ID loop (within broader sector)

## collapse all profiles in the list into one table (8760 rows per profile ID)
#dt_sector_profiles_yday  <- rbindlist(l_sector_profiles_yday,  use.names = T, fill = T)
#dt_sector_profiles_month <- rbindlist(l_sector_profiles_month, use.names = T, fill = T)

#################
## summarise the hourly emissions by the sector. This can be returned later and also check against NAEI total.
#dt_sector_yday_totals <- dt_sector_profiles_yday[, .(sector_emission = sum(hour_emis_adj, na.rm=T)), 
#                                            by= .(DateTime, Sector)]

#dt_sector_month_totals <- dt_sector_profiles_month[, .(sector_emission = sum(hour_emis_adj, na.rm=T)), 
#                                                 by= .(DateTime, Sector)]

# set to tonnes
#dt_sector_yday_totals$sector_emission <- units::set_units(dt_sector_yday_totals$sector_emission, Mg/yr)
#dt_sector_month_totals$sector_emission <- units::set_units(dt_sector_month_totals$sector_emission, Mg/yr)

# checker
#if((as.numeric(sum(dt_sector_yday_totals$sector_emission) / total_emis_checker)) < 0.999 | 
#   (as.numeric(sum(dt_sector_yday_totals$sector_emission) / total_emis_checker)) > 1.001){
#  print(paste0("Sector ",sector," hourly emissions do not add up to NAEI sector total. Check yday." ))
#}else{
#  NULL
#}

#if((as.numeric(sum(dt_sector_month_totals$sector_emission) / total_emis_checker)) < 0.999 | 
#   (as.numeric(sum(dt_sector_month_totals$sector_emission) / total_emis_checker)) > 1.001){
#  print(paste0("Sector ",sector," hourly emissions do not add up to NAEI sector total. Check month." ))
#}else{
#  NULL
#}
#################

# produce a weighted mean temporal profile for the entire sector. Format up for discussion. 
# e.g. can the ACTM take a hod profile that is day specific? or a dow profile that is month specific?
# this summary table will summarise HDD/DOY data into the MDH format. 
# 17/11/2021 : using same structure as the temporal profile data coming in. 
# 18/11/2021 : function allows for hod by dow to be written
#dt_full_sector_emis_yday <- dt_sector_profiles_yday[, .(hour_emis_adj = (sum(hour_emis_adj, na.rm=T))), 
#                                                    by= .(DateTime, yday, month, wday, hour, Sector)]
#
#if((as.numeric(sum(dt_full_sector_emis_yday$hour_emis_adj) / total_emis_checker)) < 0.999 | 
#   (as.numeric(sum(dt_full_sector_emis_yday$hour_emis_adj) / total_emis_checker)) > 1.001){
#  print(paste0("Sector ",sector," hourly emissions do not add up to NAEI sector total. Check yday." ))
#}else{
#  NULL
#}

#dt_full_sector_emis_month <- dt_sector_profiles_month[, .(hour_emis_adj = (sum(hour_emis_adj, na.rm=T))), 
#                                                      by= .(DateTime, month, wday, hour, Sector)]

### return the data for the sector for both methods. 
#l_complete <- list()

#l_complete[["yday_hourwday"]] <- dt_full_sector_emis_yday
#l_complete[["month_wday_hour"]] <- dt_full_sector_emis_month

#return(l_complete)



######################################################################################################
## function to format the profiled sector emissions into new coefficient tables, at sector level

#FormatToSectorCoeffs <- function(year, species, classification = c("SNAP","GNFR")){
#  
#  # read in the relevant .rds file
#  l_emis <- readRDS(paste0("./doc/profiled_emissions_sector/",species,"_",classification,"_",year#,"_profemis2.rds"))
#  
#  dt_emis <- rbindlist(unlist(l_emis, recursive = FALSE), use.names = T)
#  
#  l_ymwh <- list()
#  l_plot <- list()
#  
#  # go through time steps and create coefficients (hourwday is separate)
#  for(p in c("yday","month","wday","hour","hourwday")){
#    
#    if(p %in% c("yday","month","hour")){
#      
#      # multiplier to retain the coefficients centred on 1
#      c_mult <- ifelse(p == "yday", 365, ifelse(p == "month", 12, ifelse(p == "wday", 7, 24)))
#      
#      dt_p <- dt_emis[, .(p_sum = sum(hour_emis_adj, na.rm=T)), by=.(get(p), Sector)]
#      setnames(dt_p, "get", p)
#      
#      if(sum(dt_p$p_sum) != sum(dt_emis$hour_emis_adj)){
#        print(paste0("The ",p," emissions are not adding to original total. CHECK"))
#        break
#      }
#      
#      dt_p[, coeff := (p_sum/sum(p_sum)) * c_mult, by = Sector]
#      dt_p[, coeff := set_units(dt_p$coeff, NULL)]
#      
#      l_ymwh[[p]] <- dt_p
#      
#      g1 <- ggplot()+
#        geom_line(data=dt_p, aes(x=get(p), y=coeff), colour="red")+
#        labs(x=p,y="Coefficient")+
#        facet_wrap(~Sector, nrow=2)+
#        theme_bw()+
#        theme(axis.text.x = element_text(size=12),
#              axis.text.y = element_text(size=12),
#              axis.title = element_text(size=16),
#              legend.title=element_blank())
#      
#      l_plot[[p]] <- ggplotGrob(g1)
#      
#    }else if(p=="wday"){
#      
#      # multiplier to retain the coefficients centred on 1
#      c_mult <- ifelse(p == "yday", 365, ifelse(p == "month", 12, ifelse(p == "wday", 7, 24)))
#      
#      dt_p <- dt_emis[, .(p_sum = sum(hour_emis_adj, na.rm=T)), by=.(get(p), Sector)]
#      setnames(dt_p, "get", p)
#      
#      if(sum(dt_p$p_sum) != sum(dt_emis$hour_emis_adj)){
#        print(paste0("The ",p," emissions are not adding to original total. CHECK"))
#        break
#      }
#      
#      dt_p[, coeff := (p_sum/sum(p_sum)) * c_mult, by = Sector]
#      dt_p[, coeff := set_units(dt_p$coeff, NULL)]
#      
#      l_ymwh[[p]] <- dt_p
#      
#      g1 <- ggplot()+
#        geom_line(data=dt_p, aes(x=get(p), y=coeff), colour="red")+
#        labs(x=p,y="Coefficient")+
#        facet_wrap(~Sector, nrow=2)+
#        theme_bw()+
#        theme(axis.text.x = element_text(size=12),
#              axis.text.y = element_text(size=12),
#              axis.title = element_text(size=16),
#              legend.title=element_blank())
#      
#      l_plot[[p]] <- ggplotGrob(g1)
#      
#    }else{
#      
#      # create coefficients for hourwday
#      dt_p <- dt_emis[, .(p_sum = sum(hour_emis_adj, na.rm=T)), by=.(hour, wday, Sector)]
#      
#      if(sum(dt_p$p_sum) != sum(dt_emis$hour_emis_adj)){
#        print(paste0("The ",p," emissions are not adding to original total. CHECK"))
#        break
#      }
#      
#      dt_p[, coeff := (p_sum/sum(p_sum)) * 24, by = .(Sector, wday)]
#      dt_p[, coeff := set_units(dt_p$coeff, NULL)]
#      
#      l_ymwh[["hourwday"]] <- dt_p
#      
#      g1 <- ggplot()+
#        geom_line(data=dt_p, aes(x=hour, y=coeff), colour="red")+
#        labs(x="hour (by wday)",y="Coefficient")+
#        facet_grid(wday~Sector)+
#        theme_bw()+
#        theme(axis.text.x = element_text(size=12),
#              axis.text.y = element_text(size=12),
#              axis.title = element_text(size=16),
#              legend.title=element_blank())
#      
#      l_plot[[p]] <- ggplotGrob(g1)
#    }
#    
#  }
#  
#  ## save data
#  saveRDS(l_ymwh, paste0("./doc/coefficients_sector/",species,"_",classification,"_",y_emis,"_coefficients.rds"))
#  
#  ## plot data
#  d1 <- plot_grid(l_plot[["hour"]], l_plot[["wday"]], l_plot[["month"]], l_plot[["yday"]], l_plot[["hourwday"]], ncol=1, nrow=5,rel_heights = c(1,1,1,1,3), align="h", axis="l")
#  
#  save_plot(paste0("./doc/coefficients_sector/",species,"_",classification,"_",y_emis,"_plots_month.png"), d1, base_height = 18, base_width = 14)
#  
#  
#}
