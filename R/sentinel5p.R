
###############################################################################################

#### Sentinel 5p comparisons #### WORKING TO SNAP SECTORS ONLY 08/02/22

######################################################################################################
#### function to process NAEI UK emissions data to AOI ####
AOIEmissionsData <- function(aoi_name, species, year){
  
  ### extent of London study area + an expanded area for plotting interest (to ID large nearby sources?)
  aoi <- suppressWarnings(suppressMessages(st_read("./data/AOI_extents",paste0("aoi_",aoi_name), quiet = T)))
  aoi_bng <- st_transform(aoi, crs="EPSG:27700")
  
  aoi_expand <-            ext(ext(aoi_bng)[1] - (ext(aoi_bng)[2] - ext(aoi_bng)[1])/2,
                               ext(aoi_bng)[2] + (ext(aoi_bng)[2] - ext(aoi_bng)[1])/2,
                               ext(aoi_bng)[3] - (ext(aoi_bng)[4] - ext(aoi_bng)[3])/2,
                               ext(aoi_bng)[4] + (ext(aoi_bng)[4] - ext(aoi_bng)[3])/2)
  
  ### for each SNAP, read in year/spec and crop/stack
  NAEI_diff_files <- list.files(paste0("//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/NAEI_data_and_SNAPS/Emissions_grids_plain/BNG/nox/diffuse/",year,"/rasters_SNAP"), pattern = paste0("^",dt_pollutants[upper_name==species,file_name],"_diff_",year,"_uk_SNAP.*.t_1km_BNG_2019NAEImap.tif$"), full.names = T)
  
  NAEI_pt_files <- fread(paste0("//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/NAEI_data_and_SNAPS/Emissions_grids_plain/BNG/nox/point/",year,"/",dt_pollutants[upper_name==species,file_name],"_pt_",year,"_uk_SNAP_t_BNG.csv"))
  
  # run function for all SNAPS to return cropped rasters and stats etc
  
  l_AOI_emis <- lapply(paste0("S",1:11) , NAEItoAOI, NAEI_diff_files, NAEI_pt_files, aoi_bng, aoi_expand, species)
  
  # bind the SNAP AOI totals and make a sum raster (tonnes)
  dt_SNAP_tots <- rbindlist(lapply(l_AOI_emis, '[[',"stats"), use.names = T)
  dt_SNAP_tots[, SNAP := as.numeric(SNAP)]
  
  # biggest point sources
  top5_pts <- NAEI_pt_files[Easting <= aoi_expand[2] & Easting >= aoi_expand[1] & Northing <= aoi_expand[4] & Northing >= aoi_expand[3]][order(-Emission)][1:5]
  
  #top5_pts <- NAEI_pt_names[Pollutant=="NOx"][top5_pts, on=c("Easting","Northing")]
  
  # total snaps
  r_diffuse <- app(rast(lapply(l_AOI_emis, '[[',"aoi_diff")), sum, na.rm=T)
  r_diffuseE <- app(rast(lapply(l_AOI_emis, '[[',"aoiE_diff")), sum, na.rm=T)
  r_point <- app(rast(lapply(l_AOI_emis, '[[',"aoi_pt")), sum, na.rm=T)
  r_pointE <- app(rast(lapply(l_AOI_emis, '[[',"aoiE_pt")), sum, na.rm=T)
  r_total <- app(rast(lapply(l_AOI_emis, '[[',"aoi_total")), sum, na.rm=T)
  r_totalE <- app(rast(lapply(l_AOI_emis, '[[',"aoiE_total")), sum, na.rm=T)
  
  st_aoi <- c(r_diffuse, r_point, r_total)
  names(st_aoi) <- c("diffuse","point","total")
  st_exp <- c(r_diffuseE, r_pointE, r_totalE)
  names(st_exp) <- c("diffuse","point","total")
  
  return(list("dt_aoi_emis" = dt_SNAP_tots, "st_aoi" = st_aoi, "st_expanded" = st_exp, "top5_pts" = top5_pts))
  
}

###########################################################################################
#### function to take NAEI diffuse and point data and crop it to the AOI (+ expanded)  ####

NAEItoAOI <- function(SN, NAEI_diff_files, NAEI_pt_files, aoi_bng, aoi_expand, species){

  print(paste0(SN))
  
  r_file <- NAEI_diff_files[grep(paste0("_",SN,"_"), NAEI_diff_files)]
  r <- suppressWarnings(rast(r_file))
  
  pts <- NAEI_pt_files[SNAP==as.numeric(substr(SN,2,nchar(SN)))]
  
  if(nrow(pts)>0){
    
    pt_xy <- vect(pts, geom=c("Easting", "Northing"), crs = "EPSG:27700")
    r_pt <- terra::rasterize(pt_xy, r_uk_BNG, field="Emission", fun=`sum`)
    
  }else{
    r_pt <- r_uk_BNG
  }
  
  # crop to AOI and expanded AOI
  rc_dif <- crop(r, aoi_bng)
  rc_difE <- crop(r, aoi_expand)
  
  rc_pt <- crop(r_pt, aoi_bng)
  rc_ptE <- crop(r_pt, aoi_expand)
  
  rc_SNtot <- app(c(rc_dif, rc_pt), sum, na.rm=T)
  names(rc_SNtot) <- paste0(species,"_",SN)
  rc_SNtotE <- app(c(rc_difE, rc_ptE), sum, na.rm=T)
  names(rc_SNtotE) <- paste0(species,"_",SN)
  
  # list the results
  l <- list(rc_dif, rc_difE, rc_pt, rc_ptE, rc_SNtot, rc_SNtotE)
  names(l) <- c("aoi_diff","aoiE_diff","aoi_pt","aoiE_pt","aoi_total","aoiE_total")
  
  # totals
  tot_SNdif <- global(rc_dif, sum, na.rm=T)
  tot_SNdifE <- global(rc_difE, sum, na.rm=T)
  tot_SNpt  <- global(rc_pt, sum, na.rm=T)
  tot_SNptE  <- global(rc_ptE, sum, na.rm=T)
  tot_SNAP  <- global(rc_SNtot, sum, na.rm=T)
  tot_SNAPE  <- global(rc_SNtotE, sum, na.rm=T)
  
  # stats table
  dt_row <- data.table(AOI = aoi_name, Year=year, Pollutant=species, SNAP = substr(SN, 2, nchar(SN)), emission_diff = tot_SNdif, emission_diff_exp = tot_SNdifE, emission_pt = tot_SNpt, emission_pt_exp = tot_SNptE, emission_total = tot_SNAP, emission_total_exp = tot_SNAPE)
  
  l[[paste0("stats")]] <- dt_row
  
  return(l)
  
}



###########################################################################################
#### function to process AOI emissions data into one temporal profile (SNAP weighted)  ####
#### This is only done on the yday/wday GAMs, as sentinel5p does not cover hourly res  ####
#### all of this is just working to SNAP sector at the moment. to be developed         ####

AOITemporalProfile <- function(aoi_emis, species, year, write){
  
  ## bring in profile data - GAMs generated at the sector level
  ## either use specific-species GAMs or the 'general-pollutant' GAMs
  
  # make a blank calendar
  #dt_calendar <- data.table(datetime = seq(from=as.POSIXct(paste0(year,"-1-1 0:00"), tz="UTC"),to=as.POSIXct(paste0(year,"-12-31 23:00"), tz="UTC"), by="hour" )  )
  dt_calendar <- data.table(date = seq(from=as.POSIXct(paste0(year,"-1-1")),to=as.POSIXct(paste0(year,"-12-31")), by="day" ))
  dt_calendar[, c("yday","wday") := list(yday(date),wday(date, week_start = getOption("lubridate.week.start", 1)))]
  
  # call coefficient tables and add to calendar, work out emission split. Across all SNAPS into list. 
  l_SNAP_calendar <- lapply(1:11, SNAP_TP, aoi_emis, species, year, dt = dt_calendar)
  
  ## for plotting, the hour/wday fluctuation is hideous to look at
  
  dt_SNAP_calendar <- rbindlist(l_SNAP_calendar, use.names = T)
  dt_SNAP_calendar[, sector := factor(sector, levels = 1:11)] 
  dt_total_calendar <- dt_SNAP_calendar[, .(emis_t = sum(E_adjusted, na.rm=T)), by=.(date, yday, wday)]
  dt_total_calendar[, emis_alpha := emis_t/mean(emis_t, na.rm=T)]
  
  
                    
  ggplot()+
    geom_smooth(data=dt_total_calendar, aes(x=yday, y=emis_alpha), size=0.8)+
    #geom_line(data=dt_total_calendar, aes(x=yday, y=emis_alpha), size=0.8)+
    #geom_line(data=dt_SNAP_calendar, aes(x=yday, y=E_alpha, group=sector, colour=sector), alpha=0.6)+
    stat_smooth(data=dt_SNAP_calendar, aes(x=yday, y=E_alpha, group=sector, colour=sector),geom='line', se=F, alpha=0.6)+
    geom_hline(yintercept = 1, linetype="dashed")+
    theme_bw()
  
  
  
  
  
}



############################################################################################
#### function to split out SNAP level emissions onto a calendar, using GAM coefficients ####

SNAP_TP <- function(SNAPsec, aoi_emis, species, year, dt = dt_calendar){
  
  # call coefficient tables and add to calendar
  for(ts in c("wday","yday")){
    
    dt_temporal <- fread(paste0("./output/coeffs_sector/SNAP/",species,"/GAM_",ts,"_SNAP_",species,"_allYr_LIST.csv"))
    dt_temporal <- dt_temporal[sector == SNAPsec]
    dt_temporal[,c("coeff_l","coeff_u") := NULL]
    setnames(dt_temporal, "coeff", paste0(ts,"_coeff"))
    
    #if(ts == "hourwday"){
    #  dt <- dt_temporal[dt, on = c("hour","wday")] 
    #}else{
      dt <- dt_temporal[dt, on = paste0(ts)] 
    #}
    
  }
  
  dt <- dt[,c("date","sector","yday","wday","yday_coeff","wday_coeff")]
  # get the sector emissions and divide into one standard day. 
  E_t_hr <- aoi_emis[SNAP==SNAPsec & Year == year, emission_total.sum]/365
  dt[,E_t_hr :=  E_t_hr]
  
  # now apply the scaling alphas to get exact hour emission. 
  dt[, E_adjusted := E_t_hr * yday_coeff * wday_coeff]
  
  dt[, E_alpha := E_adjusted/mean(E_adjusted, na.rm=T)]
  
  return(dt)
  
} # end of function


############################################################################################
#### function to assess the raw Sentinel5p data and fit a GAM to ALL data points/days   ####

funcName <- function(){
  
}

## all point data
dt <- fread("C:/Users/samtom/OneDrive - UKCEH/DUKEMS/Sentinel5P/london/london_all_pts_values_nitrogendioxide_tropospheric_column.csv")
dt <- dt[year(date) == 2019]
dt <- dt[pt_val>=0]

#ggplot(data = dt[,.N,by=yday], aes(x=yday,y=N))+
#  geom_point()

dt_GAM <- dt[,c("yday","pt_val","qa_val")]
setnames(dt_GAM,"yday","time")
dt_GAM <- dt_GAM[qa_val >= 0.75]
dt_GAM[, coeff := pt_val/mean(pt_val, na.rm=T)]
dt_GAM[,name := "All measurements"]

gam_S5P <- gam(pt_val ~ s(time, bs="cc", k=15), data = dt_GAM, method="REML")
#gam_S5P <- gam(pt_val ~ s(time, bs="cc"), data = dt_GAM, method="REML")
dt_fit <- blankResponse(v_sectors="Sentinel5p", sectorColName="sector", timescale="yday")
dt_coeffs <- SectorCoeffs(gam = gam_S5P, dt = dt_fit, timescale = "yday")

## Claire's mean of means data
dtMean <- fread("C:/Users/samtom/OneDrive - UKCEH/DUKEMS/Sentinel5P/london/london_day_agg_nitrogendioxide_tropospheric_column_2019.csv")
setnames(dtMean,"yday","time")

l_ydaymean_gams <- list()

for(i in c("yday_mean_of_pts_summed_NO2_value","yday_mean_of_pts_mean_NO2_value")){
  
  cols <- c("time",i)
  dt_GAMmn <- dtMean[, ..cols]
  setnames(dt_GAMmn, paste0(i), "N")
  dt_GAMmn <- dt_GAMmn[N >= 0]
  dt_GAMmn[, coeff := N/mean(N, na.rm=T)]
  
  gam_Mns <- gam(N ~ s(time, bs="cc"), data = dt_GAMmn, method="REML")
  dt_fit <- blankResponse(v_sectors=i, sectorColName="sector", timescale="yday")
  dt_coeffsMns <- SectorCoeffs(gam = gam_Mns, dt = dt_fit, timescale = "yday")
  
  l_ydaymean_gams[[i]] <- dt_coeffsMns
  
}

dt_ydaymean_gams <- rbindlist(l_ydaymean_gams, use.names = T)

dt_ydaymean_pts <- melt(dtMean[,c("time","yday_mean_of_pts_summed_NO2_value","yday_mean_of_pts_mean_NO2_value")], id.vars = "time", variable.name = "sector", value.name = "N")
dt_ydaymean_pts <- dt_ydaymean_pts[N >= 0]
dt_ydaymean_pts[, coeff := N/mean(N, na.rm=T), by=sector]


## data per image
dt_image <- fread("C:/Users/samtom/OneDrive - UKCEH/DUKEMS/Sentinel5P/london/london_all_image_nitrogendioxide_tropospheric_column.csv")
setnames(dt_image,"yday","time")
dt_image <- dt_image[year(date) == 2019]

l_imgmean_gams <- list()

for(i in c("pts_sum_NO2_value","pts_mean_NO2_value")){
  
  cols <- c("time",i)
  dt_GAMimg <- dt_image[, ..cols]
  setnames(dt_GAMimg, paste0(i), "N")
  dt_GAMimg <- dt_GAMimg[N >= 0]
  dt_GAMimg[, coeff := N/mean(N, na.rm=T)]
  
  gam_img <- gam(N ~ s(time, bs="cc"), data = dt_GAMimg, method="REML")
  dt_fit <- blankResponse(v_sectors=i, sectorColName="sector", timescale="yday")
  dt_coeffsimg <- SectorCoeffs(gam = gam_img, dt = dt_fit, timescale = "yday")
  
  l_imgmean_gams[[i]] <- dt_coeffsimg
  
}

dt_imgmean_gams <- rbindlist(l_imgmean_gams, use.names = T)

dt_imgmean_pts <- melt(dt_image[,c("time","pts_sum_NO2_value","pts_mean_NO2_value")], id.vars = "time", variable.name = "sector", value.name = "N")
dt_imgmean_pts <- dt_imgmean_pts[N >= 0]
dt_imgmean_pts[, coeff := N/mean(N, na.rm=T), by=sector]

###########


## plot of modelled TP, mean data GAMs and all point data GAM
g1 <- ggplot()+
  geom_point(data = dt_ydaymean_pts, aes(x=time,y=coeff,group=sector,colour=sector))+
  geom_line(data=dt_ydaymean_gams, aes(x=yday,y=coeff,group=sector,colour=sector))+
  geom_smooth(data=dt_total_calendar, aes(x=yday, y=emis_alpha), size=0.8)+
  geom_line(data=dt_coeffs, aes(x=yday,y=coeff), colour="black")+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_vline(xintercept = 90, linetype="dashed")+
  geom_vline(xintercept = 181, linetype="dashed")+
  geom_vline(xintercept = 273, linetype="dashed")+
  labs(x="yday",y="Coefficient")+
  ylim(0,5)+
  theme_bw()+
  theme(legend.position = "bottom")

g2 <- ggplot()+
  geom_point(data = dt_imgmean_pts, aes(x=time,y=coeff,group=sector,colour=sector))+
  geom_line(data=dt_imgmean_gams, aes(x=yday,y=coeff,group=sector,colour=sector))+
  geom_smooth(data=dt_total_calendar, aes(x=yday, y=emis_alpha), size=0.8)+
  geom_line(data=dt_coeffs, aes(x=yday,y=coeff), colour="black")+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_vline(xintercept = 90, linetype="dashed")+
  geom_vline(xintercept = 181, linetype="dashed")+
  geom_vline(xintercept = 273, linetype="dashed")+
  labs(x="yday",y="")+
  ylim(0,5)+
  theme_bw()+
  theme(legend.position = "bottom")


g3 <- ggplot()+
  geom_point(data=dt_GAM, aes(x=time,y=coeff,group=name, colour=name), alpha=0.3)+
  geom_line(data=dt_coeffs, aes(x=yday,y=coeff))+
  geom_line(data=dt_total_calendar, aes(x=yday, y=emis_alpha),  alpha=0.5)+
  geom_smooth(data=dt_total_calendar, aes(x=yday, y=emis_alpha), size=0.8)+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_vline(xintercept = 90, linetype="dashed")+
  geom_vline(xintercept = 181, linetype="dashed")+
  geom_vline(xintercept = 273, linetype="dashed")+
  labs(x="yday",y="")+
  ylim(0,5)+
  theme_bw()+
  theme(legend.position = "bottom")
  

d1 <- plot_grid(g1,g2,g3, ncol=3, nrow=1,rel_heights = 1, rel_widths = 1, align="h", axis="l")

save_plot(paste0("C:/FastProcessingSam/DUKEMS/Data/RoadTransport/differing_S5P_GAMs.png"), d1, base_height = 12, base_width = 16)

ggplot()+
  geom_smooth(data=dt_total_calendar, aes(x=yday, y=emis_alpha), size=0.8)+
  #geom_line(data=dt_total_calendar, aes(x=yday, y=emis_alpha), size=0.8)+
  #geom_line(data=dt_SNAP_calendar, aes(x=yday, y=E_alpha, group=sector, colour=sector), alpha=0.6)+
  stat_smooth(data=dt_SNAP_calendar, aes(x=yday, y=E_alpha, group=sector, colour=sector),geom='line', se=F, alpha=0.6)+
  geom_line(data=dt_coeffs, aes(x=yday,y=coeff), colour="red")+
  geom_line(data=dt_coeffsMns, aes(x=yday,y=coeff), colour="red", linetype="dashed")+
  geom_point(data = dtMean, aes(x=yday,y=(yday_mean_of_rast_mean_NO2_value/mean(yday_mean_of_rast_mean_NO2_value,na.rm=T))), alpha=0.5)+
  geom_hline(yintercept = 1, linetype="dashed")+
  theme_bw()



  
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
