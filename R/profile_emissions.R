source("./R/workspace.R")

#############################################
#### TEMPORAL PROFILING ANNUAL EMISSIONS ####
#############################################

## Taking NAEI annual totals for a species and splitting into hourly fractions.     ##
## This is done with the MOY/DOW/HOD profiles. Match a profile to an NFR and split. ##
## Then sum up the NFRs into SNAPs and create a temporal profile fo the SNAP        ##
## Also create a total emissions annual profile, for comparison with Sentinel-5p.   ##

## set the year to model and the species
y_emis <- 2018 # the emissions year
## !! Need a way of choosing NFR/fuel combos to be modelled on specific year data, as opposed to means !! ##
y_spec_NFR <- NULL # the NFR codes that are to have year-specific profiles (otherwise `average` year = 0)
species <- "NOx"

## read and match the NAEI data to the temporal profile in the lookup table, along with SNAP sector
dt_naei_profs <- JoinNAEItoProfiles(year = y_emis, species = species) 
  
## create SNAP-wide temporal profiles, a weighted profile from all NFR & Profile combinations
v_sectors <- list(1)

l_DUKEMs_profiles <- lapply(X=setNames(v_sectors, v_sectors), FUN = EmissionsProfileBySector, year = y_emis, species = species, classification = "SNAP", emis = dt_naei_profs, yr_spec_NFR=NULL, hod_by_dow=F, hour_emis=T)

# enact year specific profile choice
# enact pollutant specific choice (L~140)


#!!# The sum of hourly emissions != original annual emissions:
##   This is because nearly all months are not wholly divisible by complete weeks. 
##   The week's worth of emissions, as defined by E x moy x 7/31(/30), is split into days (the sum of which = a week)
##   Full weeks are fine. The parts/days left over use the daily coeff value against a full week value. 
##   Those day coeff values are (or can be) different to a standard 1/7 value
##   Therefore the sum of the daily parts != month, as the extra days are not behaving exactly as 1/7
##   Some months this is more (if the mean extra days coeff > 1/7), some is less (if the mean extra days coeff < 1/7)

#### suggest scaling the hourly emission by the monthly coefficient total ####



tp <- "moy"
x_max <- ifelse(tp=="hod", 23, ifelse(tp=="dow", 7, 12))
plotting_profs <- get(paste0("dt_prof_",tp))[Profile_ID %in% c(dt_naei_profs[SNAP == v_sectors, unique(Profile_ID)])]
emis_profs <- dt_naei_profs[SNAP == v_sectors, .(emis_kt = sum(emis_kt)) , by=Profile_ID ]

plotting_profs <- emis_profs[plotting_profs, on=c("Profile_ID")]



ggplot()+
  geom_line(data=plotting_profs, aes(x=get(tp), y=get(paste0(tp,"_coeff")), group=Profile_ID, colour=Profile_ID, size=emis_kt))+
  geom_line(data=profiles[Unit==tp], aes(x=Period, y=coeff), colour="black", linetype="dashed")+
  #facet_wrap(~dow)+
  geom_text_repel(data=plotting_profs[get(tp)==max(get(tp))], aes(x = x_max, y = get(paste0(tp,"_coeff")), group = Profile_ID,label=Profile_ID), color="black", cex = 3, direction="y", hjust=0, nudge_x = 0.5, segment.curvature = 0, segment.size = 0.2)+
  xlim(0,x_max+1)+
  theme_bw()+
  theme(legend.position = "none")




