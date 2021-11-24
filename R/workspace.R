## initialise the workspace for splitting annual emissions into hourly profiles ##

## packages and workspace
packs <- c("sp","raster","stringr","gdalUtils","rgeos","rgdal","grid","plyr","ggplot2","ggrepel","data.table","stats","readr","ggplot2","sf","lubridate","units")
require(cowplot)
lapply(packs, require, character.only = TRUE)

"%!in%" <- Negate("%in%")

## lookup table ##
dt_sect_to_prof <<- fread("./Data/Sectors.csv")

## profile tables ##
dt_prof_moy <<- melt(fread("./Data/moy_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year"), variable.name = "moy", value.name = "moy_coeff", variable.factor = F, value.factor = F ) %>% .[, moy := as.numeric(moy)]

dt_prof_dow <<- melt(fread("./Data/dow_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year"), variable.name = "dow", value.name = "dow_coeff", variable.factor = F, value.factor = F ) %>% .[, dow := as.numeric(dow)]

dt_prof_hod <<- melt(fread("./Data/hod_profiles.csv"), id.vars = c("Profile_ID","Pollutant","Year","dow"), variable.name = "hod", value.name = "hod_coeff", variable.factor = F, value.factor = F ) %>% .[, hod := as.numeric(hod)]

dt_prof_hdd <<- fread("./Data/hdd_profiles.csv")

################################################################################################
## function to 1. return the NAEI emissions data, formatted, for the given year & Species.    ##
##             2. match the annual NAEI data to the NFR codes and the temporal profile codes. ##

JoinNAEItoProfiles <- function(year, species){
  
  ########################################################
  
  # year           = *numeric* year to process. Will determine calendar structure plus year specific profiles.
  if(!is.numeric(year)) stop ("Year is not numeric")
  # species        = *character* name of air pollutant or GHG or metal etc. Needs to conform to a list of options.
  if(species %!in% c("NOx","SO2") ) stop ("Species must be in: 
                                            AP:    BaP, CO, NH3, NMVOC, NOx, SO2
                                            PM:    PM2.5, PM10
                                            GHG:   CH4, CO2, N2O
                                            Metal: Cd, Cu, Hg, Ni, Pb, Zn")
  
  #########################################################
  
  
  colskeep <- c("NFR/CRF Group","Source","Activity",year)
  
  # read NAEI emissions data - currently on latest year = 2019. way to automate this?
  dt_naei <- fread(paste0("//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/NAEI_data_and_SNAPS/NAEI_data/diffuse/",tolower(species),"/NFC_time_series/naei_",tolower(species),"_1970-2019.csv"), header=T)
  
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
  if(species %!in% c("NOx","SO2") ) stop ("Species must be in: 
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
  
  dt_calendar[ , c("hdd","doy", "moy", "dow", "hod") := list(yday(DateTime), yday(DateTime), month(DateTime), lubridate::wday(DateTime, week_start = getOption("lubridate.week.start", 1)), hour(DateTime))]
  
  ## extract the sector relevant NFR codes from NAEI data. Aggregate by the profile. 
  dt_sect_specific <- emis[SNAP == sector]
  dt_sect_specific_agg <- dt_sect_specific[, .(emission = sum(emission, na.rm=T)), by = .(Sector = get(classification), Profile_ID)]
  
  # sum of sector emissions (for checking at end)
  total_emis_checker <- sum(dt_sect_specific$emission)
  units(total_emis_checker) <- "Gg /yr"
  
  ## empty list for the hourly emissions for the profile, per SNAP/GNFR (might be > 1 sector in a profile ID)
  l_sector_profiles <- list()
    
  ## loop through the unique temporal profiles and disaggregate emissions to hourly
  for(prof in unique(dt_sect_specific_agg[,Profile_ID])){
    
    # data needs to be processed differently if it requires HDD information
    if(prof %in% c(unique(emis[HDD==1,Profile_ID]))){
      
      # HDD processing is currently limited to domestic and commercial type combustion
      dt_prof_hdd_pollyear <- copy(dt_prof_hdd)
      if(prof %in% yr_spec_NFR){dt_prof_hdd_pollyear <- dt_prof_hdd_pollyear[Year == year]}else{dt_prof_hdd_pollyear <- dt_prof_hdd_pollyear[Year == 0]}
      dt_prof_hdd_pollyear[,c("Pollutant","Year","Profile_ID") := NULL]
      
      dt_prof_hod_pollyear <- copy(dt_prof_hod)
      dt_prof_hod_pollyear <- dt_prof_hod_pollyear[Profile_ID==prof]
      if(prof %in% yr_spec_NFR){dt_prof_hod_pollyear <- dt_prof_hod_pollyear[Year == year]}else{dt_prof_hod_pollyear <- dt_prof_hod_pollyear[Year == 0]}
      dt_prof_hod_pollyear[,c("Pollutant","Year","Profile_ID") := NULL]
      
      # join the profiles to the hourly calendar
      dt_hours_coeffs <- dt_calendar[dt_prof_hdd_pollyear, on = "hdd"][dt_prof_hod_pollyear, on = c("dow","hod")]
      
      # total emissions for the sector and profile ID
      Ekt <- dt_sect_specific_agg[Sector == sector & Profile_ID == prof, emission]
      dt_hours_coeffs[, ann_emis_kt := Ekt]
      
      if(species == "BaP"){
        units(dt_hours_coeffs$ann_emis_kt) <- "Kg /yr"
      }else{
        units(dt_hours_coeffs$ann_emis_kt) <- "Gg /yr"
      }
      
      # split emissions out into hours and day of the year
      # Readjust to total (very slightly out due to proportion of hod/dow in year)
      dt_hours_coeffs[,     hod_emis := ann_emis_kt * hdd_coeff * hod_coeff]
      dt_hours_coeffs[, hod_emis_adj := ( (hdd_coeff * ann_emis_kt) / sum(hod_emis, na.rm=T)) * hod_emis, by = hdd]
      
      # add some information
      dt_hours_coeffs[,Sector := sector]
      dt_hours_coeffs[,Profile_ID := prof]
      setkey(dt_hours_coeffs, DateTime)
      
      l_sector_profiles[[paste0("Sector_",sector," Profile_",prof)]] <- dt_hours_coeffs
      
            
    }else{
      
      # subset the temporal profile tables to the specific profile ID in the loop.
      # use this bit to subset to species and year specific profiles. 
      dt_prof_moy_pollyear <- copy(dt_prof_moy)
      dt_prof_moy_pollyear <- dt_prof_moy_pollyear[Profile_ID==prof]
      if(prof %in% yr_spec_NFR){dt_prof_moy_pollyear <- dt_prof_moy_pollyear[Year == year]}else{dt_prof_moy_pollyear <- dt_prof_moy_pollyear[Year == 0]}
      dt_prof_moy_pollyear[,c("Pollutant","Year","Profile_ID") := NULL]
      
      dt_prof_dow_pollyear <- copy(dt_prof_dow)
      dt_prof_dow_pollyear <- dt_prof_dow_pollyear[Profile_ID==prof]
      if(prof %in% yr_spec_NFR){dt_prof_dow_pollyear <- dt_prof_dow_pollyear[Year == year]}else{dt_prof_dow_pollyear <- dt_prof_dow_pollyear[Year == 0]}
      dt_prof_dow_pollyear[,c("Pollutant","Year","Profile_ID") := NULL]
      
      dt_prof_hod_pollyear <- copy(dt_prof_hod)
      dt_prof_hod_pollyear <- dt_prof_hod_pollyear[Profile_ID==prof]
      if(prof %in% yr_spec_NFR){dt_prof_hod_pollyear <- dt_prof_hod_pollyear[Year == year]}else{dt_prof_hod_pollyear <- dt_prof_hod_pollyear[Year == 0]}
      dt_prof_hod_pollyear[,c("Pollutant","Year","Profile_ID") := NULL]
      
      # join the profiles to the hourly calendar
      dt_hours_coeffs <- dt_calendar[dt_prof_moy_pollyear, on = "moy"][dt_prof_dow_pollyear, on = "dow"][dt_prof_hod_pollyear, on = c("dow","hod")]
      
      # total emissions for the sector and profile ID
      Ekt <- dt_sect_specific_agg[Sector == sector & Profile_ID == prof, emission]
      dt_hours_coeffs[, ann_emis_kt := Ekt]
      
      if(species == "BaP"){
        units(dt_hours_coeffs$ann_emis_kt) <- "Kg /yr"
      }else{
        units(dt_hours_coeffs$ann_emis_kt) <- "Gg /yr"
      }
      
      # split emissions out into hours, days, months.
      # Readjust to total (very slightly out due to proportion of days not in full week in a month)
      dt_hours_coeffs[,     hod_emis := ann_emis_kt * moy_coeff * (7/n_mon_days) * dow_coeff * hod_coeff]
      dt_hours_coeffs[, hod_emis_adj := ( (moy_coeff * ann_emis_kt) / sum(hod_emis, na.rm=T)) * hod_emis, by = moy]
      
      # add some information
      dt_hours_coeffs[,Sector := sector]
      dt_hours_coeffs[,Profile_ID := prof]
      setkey(dt_hours_coeffs, DateTime)
      
      l_sector_profiles[[paste0("Sector_",sector," Profile_",prof)]] <- dt_hours_coeffs
      
    } # end of ifelse statement for HDD or M/D/H
    
  } # end of unique profile ID loop (within broader sector)
  
  ## collapse all profiles in the list into one table (8760 rows per profile ID)
  dt_sector_profiles <- rbindlist(l_sector_profiles, use.names = T, fill = T)
  
  ######
  ## summarise the hourly emissions by the sector. This can be returned later and also check against NAEI total.
  dt_sector_hour_emis <- dt_sector_profiles[, .(sector_emission = sum(hod_emis_adj, na.rm=T)), 
                                              by= .(DateTime, Sector)]
  # set to tonnes
  dt_sector_hour_emis$sector_emission <- units::set_units(dt_sector_hour_emis$sector_emission, Mg/yr)
  
  # checker
  if((as.numeric(sum(dt_sector_hour_emis$sector_emission) / total_emis_checker)) < 0.999 | 
     (as.numeric(sum(dt_sector_hour_emis$sector_emission) / total_emis_checker)) > 1.001){
    print(paste0("Sector ",sector," hourly emissions do not add up to NAEI sector total. Check." ))
  }else{
    NULL
  }
  ######
  
  # produce a weighted mean temporal profile for the entire sector. Format up for discussion. 
  # e.g. can the ACTM take a hod profile that is day specific? or a dow profile that is month specific?
  # this summary table will summarise HDD/DOY data into the MDH format. 
  # 17/11/2021 : using same structure as the temporal profile data coming in. 
  # 18/11/2021 : function allows for hod by dow to be written
  dt_full_sector_emis <- dt_sector_profiles[, .(hod_emis_adj = (sum(hod_emis_adj, na.rm=T))), 
                                            by= .(DateTime, moy, dow, hod, Sector)]
  
  ## return the data as a list of MOY, DOW, HOD. Multiple sectors gives list of lists. 
  ## choices for format in the function. 
  l_mdh <- list()
  
  dt_sector_moy_profile <- dt_full_sector_emis[, .(moy_sum = sum(hod_emis_adj, na.rm=T)), by=.(moy)]
  dt_sector_moy_profile[, c("Sector","coeff") := list(sector, moy_sum / sum(moy_sum))]
  dt_sector_moy_profile <- dt_sector_moy_profile[, c("Sector","moy","coeff")]
  
  dt_sector_dow_profile <- dt_full_sector_emis[, .(dow_sum = sum(hod_emis_adj, na.rm=T)), by=.(dow)]
  dt_sector_dow_profile[, c("Sector","coeff") := list(sector, dow_sum / sum(dow_sum))]
  dt_sector_dow_profile <- dt_sector_dow_profile[, c("Sector","dow","coeff")]
  
  
  if(hod_by_dow == T){
    dt_sector_hod_profile <- dt_full_sector_emis[, .(hod_sum = sum(hod_emis_adj, na.rm=T)), by=.(hod, dow)]
    dt_sector_hod_profile[, c("Sector","coeff") := list(sector, hod_sum / sum(hod_sum))]
    dt_sector_hod_profile <- dt_sector_hod_profile[, c("Sector","hod","dow","coeff")]
  }else{
    dt_sector_hod_profile <- dt_full_sector_emis[, .(hod_sum = sum(hod_emis_adj, na.rm=T)), by=.(hod)]
    dt_sector_hod_profile[, c("Sector","coeff") := list(sector, hod_sum / sum(hod_sum))]
    dt_sector_hod_profile <- dt_sector_hod_profile[, c("Sector","hod","coeff")]
  }
  
  l_mdh[["moy"]] <- dt_sector_moy_profile
  l_mdh[["dow"]] <- dt_sector_dow_profile
  l_mdh[["hod"]] <- dt_sector_hod_profile
  
  ## add in the hourly emissions to the list if that is requested
  if(hour_emis==T){
    l_mdh[["hour_emis"]] <- dt_sector_hour_emis
  }else{
    NULL
  }
  
  
  return(l_mdh)
  
} # end of function


#############################################################################################

