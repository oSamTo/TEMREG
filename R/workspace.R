## initialise the workspace for splitting annual emissions into hourly profiles ##

## packages and workspace
packs <- c("sp","raster","stringr","gdalUtils","rgeos","rgdal","grid","plyr","ggplot2","ggrepel","data.table","stats","readr","ggplot2","sf","lubridate","units", "cowplot")
lapply(packs, require, character.only = TRUE)

"%!in%" <- Negate("%in%")

## lookup table ##
dt_sect_to_prof <<- fread("./Data/Sectors.csv")

## profile tables ##
dt_prof_hour <<- melt(fread("./Data/hour_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year"), variable.name = "hour", value.name = "hour_coeff", variable.factor = F, value.factor = F ) %>% .[, hour := as.numeric(hour)]

dt_prof_hourwday <<- melt(fread("./Data/hourwday_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year","wday"), variable.name = "hour", value.name = "hrwd_coeff", variable.factor = F, value.factor = F ) %>% .[, hour := as.numeric(hour)]

dt_prof_wday <<- melt(fread("./Data/wday_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year"), variable.name = "wday", value.name = "wday_coeff", variable.factor = F, value.factor = F ) %>% .[, wday := as.numeric(wday)]

dt_prof_month <<- melt(fread("./Data/month_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year"), variable.name = "month", value.name = "month_coeff", variable.factor = F, value.factor = F ) %>% .[, month := as.numeric(month)]

dt_prof_yday <<- melt(fread("./Data/yday_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year"), variable.name = "yday", value.name = "yday_coeff", variable.factor = F, value.factor = F ) %>% .[, yday := as.numeric(yday)]

################################################################################################
## function to 1. return the NAEI emissions data, formatted, for the given year & Species.    ##
##             2. match the annual NAEI data to the NFR codes and the temporal profile codes. ##

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

######################################################################################################
## function to profile emissions into classification system, using the profile IDs per NFR sector.  ##
##  Result: a weighted mean temporal profile for the whole Sector, based on emissions per profile   ##

EmissionsProfileBySector <- function(year, species, sector, classification = c("SNAP","GNFR"), emis, yr_spec_NFR = NULL, hod_by_dow=F, hour_emis=F){
  
  ####################################################
  
  # year           = *numeric* year to process. Will determine calendar structure plus year specific profiles.
  if(!is.numeric(year)) stop ("Year is not numeric")
  # species        = *character* name of air pollutant or GHG or metal etc. Needs to conform to a list of options.
  if(species %!in% c("NOx","SOx","CH4","CO2","N2O") ) stop ("Species must be in: 
                                            AP:    BaP, CO, NH3, NMVOC, NOx, SO2
                                            PM:    PM2.5, PM10
                                            GHG:   CH4, CO2, N2O
                                            Metal: Cd, Cu, Hg, Ni, Pb, Zn")
  # sector         = *list* sectors to run, e.g. some SNAPS or some GNFR codes
  if(sector %in% 1:11 & classification != "SNAP") stop("Sector choices do not match classification system")
  if(sector %!in% 1:11 & classification == "SNAP") stop("Sector choices do not match classification system")
  if(sector %in% LETTERS[1:16] & classification != "GNFR") stop("Sector choices do not match classification system")
  if(sector %!in% LETTERS[1:16] & classification == "GNFR") stop("Sector choices do not match classification system")
  # classification = *character* GNFR or SNAP
  classification <- match.arg(classification)
  # emis           = *data.table* from JoinNAEItoProfiles; total emissions per NFR with matching profile IDs
  # yr_spec_NFR    = *character* optional vector of NFR codes to be profiled by year specific data, not mean
  if(sum(yr_spec_NFR %!in% unique(dt_sect_to_prof$NFR19))) stop(paste0("The following are not NFR codes, check: ",yr_spec_NFR[yr_spec_NFR %!in% unique(dt_sect_to_prof$NFR19)]))
  # hod_by_dow     = *logical* should hour of day coeffs be specific to the day of the week. Default = False. 
  # hour_emis      = *logical* should a table of emissions (t) by hour in year be provided. Default = False.
  # save_tp        = *logical* should the list of profiles (and emissions) be saved as csvs
  
  ####################################################
  
  print(paste0(Sys.time(),": Profiling ",classification," ",sector))
  
  ## create a blank calendar for the given year to attach profiled emissions to. 
  startDate <- dmy_hm(paste0("01/01/",year," 00:00"))
  endDate <- dmy_hm(paste0("31/12/",year," 23:00"))
  
  dt_calendar <- data.table(DateTime = seq(from = startDate, to = endDate, by = 3600))
  dt_calendar[,n_mon_days := days_in_month(DateTime)]
  
  dt_calendar[ , c("yday","month", "wday", "hour") := list(yday(DateTime), month(DateTime), lubridate::wday(DateTime, week_start = getOption("lubridate.week.start", 1)), hour(DateTime))]
  
  ## extract the sector relevant NFR codes from NAEI data. Aggregate by the profile. 
  dt_sect_specific <- emis[SNAP == sector]
  dt_sect_specific_agg <- dt_sect_specific[, .(emission = sum(emission, na.rm=T)), by = .(Sector = get(classification), Profile_ID)]
  
  # sum of sector emissions (for checking at end)
  total_emis_checker <- sum(dt_sect_specific$emission)
  units(total_emis_checker) <- "Gg /yr"
  
  ## empty list for the hourly emissions for the profile, per SNAP/GNFR (might be > 1 sector in a profile ID)
  l_sector_profiles_yday <- list()
  #l_sector_profiles_month <- list()
    
  ## loop through the unique temporal profiles and disaggregate emissions to hourly
  for(prof in unique(dt_sect_specific_agg[,Profile_ID])){
    
    ## all data has access to yday/month/wday/hour/hourwday. 
    ## using yday, month is not needed. wday is needed to know which hourwday to use. 
    
    ### method 'yday': using yday generated GAM and attaching hour, determined by wday ###
    ### the month is not used here. The yday is used and the month cab be calculated ###
    # subset the temporal profile tables to the specific profile ID in the loop.
    if(prof %in% yr_spec_NFR){
      dt_yday_ID <- dt_prof_yday[Profile_ID == prof & Year == year]
      dt_hourwday_ID <- dt_prof_hourwday[Profile_ID == prof & Year == year]
    }else{
      dt_yday_ID <- dt_prof_yday[Profile_ID == prof & Year == 0]
      dt_hourwday_ID <- dt_prof_hourwday[Profile_ID == prof & Year == 0]
    }
    dt_yday_ID[,c("Pollutant","Year","Profile_ID") := NULL]
    dt_hourwday_ID[,c("Pollutant","Year","Profile_ID") := NULL]
    
    # join the profiles to the hourly calendar
    dt_hours_coeffs <- dt_calendar[dt_yday_ID, on = "yday"][dt_hourwday_ID, on = c("wday","hour")]
    
    # total emissions for the sector and profile ID
    Ekt <- dt_sect_specific_agg[Sector == sector & Profile_ID == prof, emission]
    dt_hours_coeffs[, ann_emis_kt := Ekt]
    
    if(species == "BaP"){
      units(dt_hours_coeffs$ann_emis_kt) <- "Kg /yr"
    }else{
      units(dt_hours_coeffs$ann_emis_kt) <- "Gg /yr"
    }
    
    # split emissions out into hours and day of the year
    # Readjust to total (sometimes very slightly out due to proportion of hod/dow in year)
    dt_hours_coeffs[,     hour_emis := ann_emis_kt * (yday_coeff/365) * (hrwd_coeff/24)]
    dt_hours_coeffs[, hour_emis_adj := ( ((yday_coeff/365) * ann_emis_kt) / sum(hour_emis, na.rm=T)) * hour_emis, by = yday]
    
    # add some information
    dt_hours_coeffs[,Sector := sector]
    dt_hours_coeffs[,Profile_ID := prof]
    setkey(dt_hours_coeffs, DateTime)
    #dt_hours_coeffs[,Method := "yday"]
    
    l_sector_profiles_yday[[paste0("Sector_",sector," Profile_",prof)]] <- dt_hours_coeffs
    
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
    
  } # end of unique profile ID loop (within broader sector)
  
  ## collapse all profiles in the list into one table (8760 rows per profile ID)
  dt_sector_profiles_yday  <- rbindlist(l_sector_profiles_yday,  use.names = T, fill = T)
  #dt_sector_profiles_month <- rbindlist(l_sector_profiles_month, use.names = T, fill = T)
  
  #################
  ## summarise the hourly emissions by the sector. This can be returned later and also check against NAEI total.
  dt_sector_yday_totals <- dt_sector_profiles_yday[, .(sector_emission = sum(hour_emis_adj, na.rm=T)), 
                                              by= .(DateTime, Sector)]
  
  #dt_sector_month_totals <- dt_sector_profiles_month[, .(sector_emission = sum(hour_emis_adj, na.rm=T)), 
  #                                                 by= .(DateTime, Sector)]
  
  # set to tonnes
  dt_sector_yday_totals$sector_emission <- units::set_units(dt_sector_yday_totals$sector_emission, Mg/yr)
  #dt_sector_month_totals$sector_emission <- units::set_units(dt_sector_month_totals$sector_emission, Mg/yr)
  
  # checker
  if((as.numeric(sum(dt_sector_yday_totals$sector_emission) / total_emis_checker)) < 0.999 | 
     (as.numeric(sum(dt_sector_yday_totals$sector_emission) / total_emis_checker)) > 1.001){
    print(paste0("Sector ",sector," hourly emissions do not add up to NAEI sector total. Check yday." ))
  }else{
    NULL
  }
  
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
  dt_full_sector_emis_yday <- dt_sector_profiles_yday[, .(hour_emis_adj = (sum(hour_emis_adj, na.rm=T))), 
                                                      by= .(DateTime, yday, month, wday, hour, Sector)]
  
  if((as.numeric(sum(dt_full_sector_emis_yday$hour_emis_adj) / total_emis_checker)) < 0.999 | 
     (as.numeric(sum(dt_full_sector_emis_yday$hour_emis_adj) / total_emis_checker)) > 1.001){
    print(paste0("Sector ",sector," hourly emissions do not add up to NAEI sector total. Check yday." ))
  }else{
    NULL
  }
  
  #dt_full_sector_emis_month <- dt_sector_profiles_month[, .(hour_emis_adj = (sum(hour_emis_adj, na.rm=T))), 
  #                                                      by= .(DateTime, month, wday, hour, Sector)]
  
  ### return the data for the sector for both methods. 
  l_complete <- list()
  
  l_complete[["yday_hourwday"]] <- dt_full_sector_emis_yday
  #l_complete[["month_wday_hour"]] <- dt_full_sector_emis_month
  
  return(l_complete)
  
  
} # end of function


######################################################################################################
## function to format the profiled sector emissions into new coefficient tables, at sector level

FormatToSectorCoeffs <- function(year, species, classification = c("SNAP","GNFR")){
  
  # read in the relevant .rds file
  l_emis <- readRDS(paste0("./doc/profiled_emissions_sector/",species,"_",classification,"_",year,"_profemis.rds"))
  
  dt_emis <- rbindlist(unlist(l_emis, recursive = FALSE), use.names = T)
  
  l_ymwh <- list()
  
  # go through time steps and create coefficients (hourwday is separate)
  for(p in c("yday","month","wday","hour")){
    
    dt_emis[, .(p_sum = sum(hour_emis_adj, na.rm=T)), by=.(yday, sector)]
    
    
  }
  
  dt_sector_moy_profile <- dt_full_sector_emis[, .(moy_sum = sum(hod_emis_adj, na.rm=T)), by=.(moy)]
  #dt_sector_moy_profile[, c("Sector","coeff") := list(sector, moy_sum / sum(moy_sum))]
  #dt_sector_moy_profile <- dt_sector_moy_profile[, c("Sector","moy","coeff")]
  
  
  
}

######################################################################################################
## function to create new GAMs to the newly aggregated and profiled sector emissions

FitGAMsToSectorTotals <- function(){
  
  
  
  
}













## return the data as a list of MOY, DOW, HOD. Multiple sectors gives list of lists. 
## choices for format in the function. 
#l_mdh <- list()

#dt_sector_moy_profile <- dt_full_sector_emis[, .(moy_sum = sum(hod_emis_adj, na.rm=T)), by=.(moy)]
#dt_sector_moy_profile[, c("Sector","coeff") := list(sector, moy_sum / sum(moy_sum))]
#dt_sector_moy_profile <- dt_sector_moy_profile[, c("Sector","moy","coeff")]
#
#dt_sector_dow_profile <- dt_full_sector_emis[, .(dow_sum = sum(hod_emis_adj, na.rm=T)), by=.(dow)]
#dt_sector_dow_profile[, c("Sector","coeff") := list(sector, dow_sum / sum(dow_sum))]
#dt_sector_dow_profile <- dt_sector_dow_profile[, c("Sector","dow","coeff")]


#if(hod_by_dow == T){
#  dt_sector_hod_profile <- dt_full_sector_emis[, .(hod_sum = sum(hod_emis_adj, na.rm=T)), by=.(hod, dow)]
#  dt_sector_hod_profile[, c("Sector","coeff") := list(sector, hod_sum / sum(hod_sum))]
#  dt_sector_hod_profile <- dt_sector_hod_profile[, c("Sector","hod","dow","coeff")]
#}else{
#  dt_sector_hod_profile <- dt_full_sector_emis[, .(hod_sum = sum(hod_emis_adj, na.rm=T)), by=.(hod)]
#  dt_sector_hod_profile[, c("Sector","coeff") := list(sector, hod_sum / sum(hod_sum))]
#  dt_sector_hod_profile <- dt_sector_hod_profile[, c("Sector","hod","coeff")]
#}

#l_mdh[["moy"]] <- dt_sector_moy_profile
#l_mdh[["dow"]] <- dt_sector_dow_profile
#l_mdh[["hod"]] <- dt_sector_hod_profile

## add in the hourly emissions to the list if that is requested
#if(hour_emis==T){
#  l_mdh[["hour_emis"]] <- dt_sector_hour_emis
#}else{
#  NULL
#}

#return(l_mdh)

