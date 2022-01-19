source("C:/Users/samtom/OneDrive - UKCEH/DUKEMS/Sentinel5P/R/workspace.R")
#########################################
#### NOx totals for Clare P Sentinel ####
#########################################

####################################################
#### Compare Sentinel profile to NAEI emissions ####
#### profile, using DUKEMs tp's applied to NAEI ####
####    emissions data - look at NOx column     ####
####################################################

### Take an extent area (AOI) and crop down the UK NAEI emissions, by SNAP. 
### Then apply the pre-made temporal profiles calc'd for each SNAP under DUKEMs
### Sum the resulting hourly emissions (by SNAP) into one total NOx profile for the AOI
### The temporal profiles are sub-SNAP weighted, profiles are applied to groups of NFRs
### This means SNAP level profile might not necessarily represent the AOI, e.g. SNAP8 ships/inland etc
###      to do this, you would need sub-SNAP maps and apply the required profile ID. 

species <- "nox"
y_emis <- 2019
aoi_name <- "london" # "london", "leeds_brad", "notts_coal"

## return the NAEI emissions data in the AOI as SNAP totals and a raster
l_emis <- AOIEmissionsData(aoi_name = aoi_name, species = species, year = y_emis)

## The profiles need to be re-applied at the SNAP level
## The annual general profile (against UK total) does not apply everwhere as the AOI dictates
## the contribution of each SNAP and therefore how the temporal profile for the total looks

## work the NAEI emissions data for the AOI into total emissions profile. 
l_tp <- AOITemporalProfile(aoi_emis = l_emis[["dt_aoi_emis"]], species = species, year = y_emis, write=T)

## plot the profile and emissions data into a composite plot

PlottingAOIdata(species, year = y_emis, aoi_name, r_aoi_tot = l_emis[["st_aoi"]], r_exp_tot = l_emis[["st_expanded"]], aoi_snap_tot = l_emis[["dt_aoi_emis"]], tps = l_tp, pt_labels = l_emis[["top5_pts"]])

###################################################################

## generic SNAP tps from UK totals. 

ggSN <- ggplot(data=dt_nox_daily, aes(x=doy,y=doy_coeff))+
  geom_line()+
  #geom_point()+
  #scale_x_continuous(name="Day of Year", breaks=c(0,50,100,150,200,250,300,350))+
  #scale_y_continuous(name = "Coefficient", breaks=c(0.001,0.002,0.003,0.004), limits=c(0.001, 0.004))+
  facet_wrap(~SNAP)+
  geom_smooth()+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16))

ggsave(paste0(species,"_tp_doy_SNAP_",y_emis,".png"), ggSN, width=14, height=8)




###################################################################
###################################################################
###################################################################
## uk air data

setwd("C:/Users/samtom/OneDrive - UKCEH/DUKEMS/Sentinel5P")

library(openair)
require(data.table)
require(sf)
require(mgcv)
require(lmtest)
require(dtw)

BNG <- suppressWarnings(CRS("+init=epsg:27700"))
LL <- suppressWarnings(CRS("+init=epsg:4326"))

#### choose AOI
aoi_name <- "leeds_brad" # "london", "leeds_brad", "notts_coal"
aoi <- suppressWarnings(suppressMessages(st_read("./AOI_extents",paste0("aoi_",aoi_name), quiet = T)))
aoi_bng <- st_transform(aoi, BNG)

#### for 3 AQ networks, call all sites in the AOI and all data. Some extra variables.
l_sites_aoi <- list()

for(nw in c("AURN","AQE","KCL")){
  
  print(nw)
  
  if(nw=="KCL"){
    dt_meta <- as.data.table(importMeta(source = nw, all = TRUE))
  }else{
    dt_meta <- as.data.table(importMeta(source = nw, all = TRUE))
    dt_meta <- dt_meta[variable=="NOx"]
  }
  
  #v_sites <- unique(dt_meta[,code])
  v_sites_aoi <- dt_meta[longitude >= extent(aoi)[1] & longitude <= extent(aoi)[2] & latitude >= extent(aoi)[3] & latitude <= extent(aoi)[4], unique(code)]
  
  if(length(v_sites_aoi)==0){
    print(paste0("No sites in AOI in the ",nw," network. Next."))
    next
  }
  
  dt_sites_aoi <- dt_meta[code %in% v_sites_aoi]
  
  
  dt_aurn_nox <- as.data.table(get(paste0("import",nw))(site = v_sites_aoi, year = 2019, pollutant = "nox"))
  
  dt_aurn_nox[, code_entries := .N, by=code]
  dt_aurn_nox[, code_NAs := sum(is.na(nox)), by=code]
  dt_aurn_nox[, prop_NAs := code_NAs/code_entries]
  
  dt_aurn_nox[, network := nw]
  
  dt_aurn_nox[,year:=year(date)]
  dt_aurn_nox[,yday:=yday(date)]
  dt_aurn_nox[,hour:=hour(date)]
  
  dt_aurn_nox_xy <- dt_sites_aoi[dt_aurn_nox, on="code"]
  
  l_sites_aoi[[nw]] <- dt_aurn_nox_xy
  
}

#### bind all data together, subset cols, keep sites with >= 75% temporal coverage, create coeff
dt_aoi_meas <- rbindlist(l_sites_aoi, use.names = T, fill=T)
cols_keep <- c("code","site","network","site_type","latitude","longitude","code_entries","code_NAs","prop_NAs","date","year","yday","hour","nox")
dt_aoi_meas <- dt_aoi_meas[,..cols_keep]
dt_aoi_meas <- dt_aoi_meas[prop_NAs < 0.25]
dt_aoi_meas[,coeff := (nox/sum(nox, na.rm=T)), by=.(code,network)]

## in London remove 2 sites warping data
if(aoi_name == "london"){
  site_rem <- c("LHRBR","RHH")
  dt_aoi_meas <- dt_aoi_meas[!(code %in% site_rem)]
  
}

#### make one gam for AOI from all data (yday for now)
#gam_hr <- gam(nox ~ s(hour, bs="cc", k=23), data=dt_aoi_meas)
#gam_wd <- gam(N ~ s(wday, bs="cc", k=6),  data=dt_aoi_meas)
gam_yd <- gam(nox ~ s(yday, bs="cr"), data=dt_aoi_meas)

dt_gc <- data.table(yday=1:365)

fits = predict(gam_yd, newdata=dt_gc, type='response', se=T)
predicts = as.data.table(data.frame(dt_gc, fits) %>% 
                           mutate(lower = fit - 1.96*se.fit,
                                  upper = fit + 1.96*se.fit))

predicts[, coeff   := (fit/sum(fit))/24]
predicts[, coeff_l := (1-(se.fit/fit)) * coeff ]
predicts[, coeff_u := ((se.fit/fit)+1) * coeff ]
predicts[,date := as.Date(yday, origin = "2018-12-31")]
predicts[, date := as.POSIXct(date)]

#### Sentinel5p returns ####

# (currently mean value per image)
## why has Clare's `PROD_London_s5p_all_images_2019.png` not got all the points in the file below? is this from quality flag?

dt_sent5p <- fread(paste0("./",aoi_name,"/",aoi_name,"_all_image_nitrogendioxide_tropospheric_column.csv"))
dt_sent5p <- dt_sent5p[pts_mean_NO2_value>=0]

gam_yd <- gam(pts_mean_NO2_value ~ s(yday, bs="cc"), data=dt_sent5p)

dt_gc <- data.table(yday=1:365)

sen5p_fits = predict(gam_yd, newdata=dt_gc, type='response', se=T)
sen5p_predicts = as.data.table(data.frame(dt_gc, sen5p_fits) %>% 
                           mutate(lower = fit - 1.96*se.fit,
                                  upper = fit + 1.96*se.fit))

sen5p_predicts[, coeff   := (fit/sum(fit))/24]
sen5p_predicts[, coeff_l := (1-(se.fit/fit)) * coeff ]
sen5p_predicts[, coeff_u := ((se.fit/fit)+1) * coeff ]
sen5p_predicts[, date := as.Date(yday, origin = "2018-12-31")]
sen5p_predicts[, date := as.POSIXct(date)]

#### set some graphical params
if(aoi_name=="london"){
  ph <- 15
  pw <- 20
  pcol <- 3
}else{
  ph <- 15
  pw <- 15
  pcol <- 1
}

#### plot the site data smooth lines, facet by site type, with the overall GAM from above
gl <- ggplot()+
  #geom_line(data=dt_aoi_meas, aes(x=date,y=coeff,group=code,colour=code))+
  geom_hline(yintercept = 1/8760, linetype="dashed")+
  geom_smooth(data=dt_aoi_meas, aes(x=date,y=coeff,group=code,colour=code),se=F)+
  geom_line(data=predicts, aes(x=date, y=coeff), colour="black")+
  geom_ribbon(data=predicts, aes(x=date, ymin=coeff_l, ymax=coeff_u), alpha=0.5)+
  geom_line(data=sen5p_predicts, aes(x=date, y=coeff), colour="black", linetype="dashed", size=1)+
  geom_ribbon(data=sen5p_predicts, aes(x=date, ymin=coeff_l, ymax=coeff_u), alpha=0.5)+
  facet_wrap(~site_type, ncol = pcol)+
  theme_bw()+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 14))

ggsave(paste0("./",aoi_name,"/",aoi_name,"_measurements_gam_2019.png"), gl, height=ph, width=pw)

#### compare the gams 1:1

# https://stats.stackexchange.com/questions/19103/how-to-statistically-compare-two-time-series

dt_meas_gam <- predicts[,c("yday","fit")]
dt_meas_gam[, meas_coeff  := (fit/sum(fit))]
dt_sent_gam <- sen5p_predicts[,c("yday","fit")]
dt_sent_gam[, sent_coeff  := (fit/sum(fit))]

dt_gams <- dt_sent_gam[dt_meas_gam, on="yday"]

dt_gams[,sent_min_meas := sent_coeff - meas_coeff]

ggplot()+
  geom_point(data=dt_gams, aes(x=sent_coeff,y=meas_coeff,colour=yday),size=2)+
  geom_smooth(data=dt_gams, method = lm, formula=y~x, se=FALSE,aes(x=sent_coeff,y=meas_coeff),colour="red", linetype="dashed")+
  geom_abline(slope=1, intercept=0)+
  theme_bw()

grangertest(meas_coeff ~ sent_coeff, order = 3, data = dt_gams)

alignment<-dtw(dt_gams[,meas_coeff],dt_gams[,sent_coeff],keep=TRUE)

plot(alignment,type="threeway")

## Align and plot with the Rabiner-Juang type VI-c unsmoothed recursion
plot(
  dtw(dt_gams[,meas_coeff],dt_gams[,sent_coeff],keep=TRUE,
      step=rabinerJuangStepPattern(1,"c")),
  type="twoway",offset=-2);
dtwPlotDensity(alignment)

ggplot()+
  geom_line(data=dt_gams, aes(x=yday, y=sent_coeff), colour="red")+
  geom_line(data=dt_gams, aes(x=yday, y=meas_coeff), colour="blue")+
  geom_line(data=dt_gams, aes(x=yday, y=sent_min_meas))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 1/365, linetype="dashed")+
  theme_bw()


#### extract sites and make point surface - plot as v simple map in AOI
dt_sites <- dt_aoi_meas[,c("code","site","network","site_type","latitude","longitude")]
dt_sites <- dt_sites[!duplicated(dt_sites)]
sf_sites <- sf::st_as_sf(dt_sites, coords = c("longitude", "latitude"), crs=LL)

gm <- ggplot()+
  geom_sf(data = aoi, fill = NA, colour = "black", size = 0.3)+
  geom_sf(data = sf_sites,aes(colour=network), fill = NA, size = 1.5)+
  facet_wrap(~site_type)+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        strip.text.x = element_text(size = 14))

ggsave(paste0("./",aoi_name,"/",aoi_name,"_sites_map_2019.png"), gm, height=10, width=10)

###################################################################
###################################################################
###################################################################

