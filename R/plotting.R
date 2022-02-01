

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
  
  for(timescale in c("hour","hourwday","wday","yday")){
    
    #l_plot_data <- list()
    
    if(((year %in% 2010:2019) & (species %in% dt_pollutants[,upper_name]))){
      filepath <- paste0("./output/coeffs_sector/",classification,"/",species,"/")
      filename <- paste0("GAM_",timescale,"_",classification,"_",species,"_",year,"_LIST")
    }else if(is.na(year) & (species %in% dt_pollutants[,upper_name])){
      filepath <- paste0("./output/coeffs_sector/",classification,"/",species,"/")
      filename <- paste0("GAM_",timescale,"_",classification,"_",species,"_allYr_LIST")
    }else if((year %in% 2010:2019) & is.na(species)){
      filepath <- paste0("./output/coeffs_sector/",classification,"/allGas/")
      filename <- paste0("GAM_",timescale,"_",classification,"_allGas_",year,"_LIST")
    }else{
      filepath <- paste0("./output/coeffs_sector/",classification,"/allGas/")
      filename <- paste0("GAM_",timescale,"_",classification,"_allGas_allYr_LIST")
    }
    
    # read in the data
    dt <- fread(paste0(filepath,filename,".csv"))
    #  dt[, time := timescale]
    # dt[, Year := year]
    # dt[, Gas:= s]
    if(timescale != "hourwday") setnames(dt, timescale, "time")
    if(timescale == "hourwday") dt[,wday := factor(wday, levels = 1:7)]
    
    # l_plot_data[[paste0(s,"_",year,"_",timescale)]] <- dt
    
    
    g1 <- ggplot()+
      {if(timescale == "hourwday")geom_line(data=dt, aes(x=hour,y=coeff, group=wday, colour=wday))}+
      {if(timescale != "hourwday")geom_line(data=dt, aes(x=time,y=coeff), colour="red")}+
      {if(timescale == "hourwday")geom_ribbon(data=dt, aes(x=hour, ymin=coeff_l, ymax=coeff_u, group=wday, fill=wday), alpha=0.3)}+
      {if(timescale != "hourwday")geom_ribbon(data=dt, aes(x=time, ymin=coeff_l, ymax=coeff_u), alpha=0.3)}+
      labs(x=timescale,y="Coefficient")+
      facet_wrap(~sector, nrow=2)+
      theme_bw()+
      theme(axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=12),
            axis.title = element_text(size=16),
            legend.title=element_blank())
    
    l_plot[[timescale]] <- ggplotGrob(g1)
    
  } # time loop
  
  # stitch plots together
  
  d1 <- plot_grid(l_plot[["hourwday"]], l_plot[["hour"]], l_plot[["wday"]], l_plot[["yday"]], ncol=1, nrow=4,rel_heights = c(1.5,1,1,1), rel_widths = c(1,1,1,1), align="h", axis="l")
  
  dir.create(file.path(paste0("./output/plots/",classification)), showWarnings = FALSE)
  
  save_plot(paste0("./output/plots/",classification,"/",filename,"_PLOT.png"), d1, base_height = 15, base_width = 15)
  
  print(paste0(Sys.time(),": DONE."))
  
} # end of function
