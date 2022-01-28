# DUKEMs_TP

######################################################################
#### **Processing UK annual emissions to sector-level temporal profiles.**
######################################################################

# Structural change coming: 
* Create raw-data GAMs at Profile_ID level, not aggregated level
    * this will allow for direct NFR to GAM relationship
    * skips the 'grouping' sector, e.g. Shipping, which needs updating if raw data changes. Mistakes can happen easily here. 
    * don't group `by=` in raw-date GAM anymore, so faster per GAM but many more GAMs produced from raw data
    * the sector GAM calls in many more GAMs to weight. Yet to be seen if this causes processing issues. 
* The NAEI data should be collated for all pollutants/years and output to one table.
* a generic NFR to sector contribution (i.e. weighting) should be made across all years and pollutant to represent 'generic'.
    * consider only using more recent (10?) years to year of interest
* Maybe an option to continue to produce pollutant specific? yes. Put in ./output/GAM_specific (?) and set on .gitignore
* Read and weight GAMs in a sector but build one huge table and create one GAM per timescale
    * read sector, read GAMs from NFR_to_GAM, sample them and weight by emissions contribution. 
    * add to large table. repeat for all sectors. 
    * take timescale GAM, group `by = sector` 
* Output = 1 GAM for each timescale: contains n_sector groups, representing generic pollutant in generic year. 

* All of this to go in targets?
    * if you updated the raw data GAMs, you need to update NFR_to_GAM table. Targets will spot this and re-run
    * if you update the pollutants to work over, targets should spot and re-run


*Info:*
----------------

* Raw data collected for activity data in the UK, processed to standardized output 
    * Processing of data in C:/...../DUKEMS/Database/R/gam_generation.R
    * Standard input files uploaded to ./data/GAM_input  (not domestic, too large)
* Each NFR code maps to a small grouping (Profile_ID), which maps to a GAM_input data category. 
* GAMs are produced from the raw data (./R/xxxxx.R), grouped by Profile_ID, and placed in ./data/GAM_output
    * *Need to move GAM creation function over to this Git as currently done on C:/*
* They are produced for hour, hourwday (hour for specific day of week), wday, month and yday
* These GAM outputs are the basis for sector-level groupings;
    * Need to know the NFR to sector look-up (e.g. to SNAP)
    * The NFRs are grouped by Profile_ID and the total emissions represented by the Profile_ID, as a proportion of sector total, defines a Profile weight
    * The weights influence the construction of the final sector level GAM (./output/GAM_sector)
* As well as GAMs, standard csv outputs of temporal profiles are produced; coefficients are relative to 1

* When further raw data is found and ingested, it must either map to current Profile_IDs or new IDs need to be added and mapped to NFR codes. 

----------------

**Emissions contribution by Profile_ID per sector here; ./doc/profile_contributions_sector/**

**IMPORTANT:** When making GAM of GAMs, the weighted raw data is having an effect when the scales of data are very different, which is annoying. Need to re-scale to 1 (div by mean). 

**Thought:** The year actually does determine the shape of the profile, even when not using year specific activity data, as the emissions contributions to a sector change per year (not a lot but they do). Do i need to take an average of NFR contribution? another GAM???

**Thought:** The shape of the sector GAM will differ from pollutant to pollutant as emissions are more pertinent to certain NFRs and therefore certain shapes. NH3 in SNAP 10 is the most obvious example. Yet EMEP (e.g.) just uses one profile per SNAP, all pollutants.

**IDEA to develop (19/01/22):** The method could be changed by giving the original raw data a weighting based on a look-up table. Going through each SNAP, you would calculate what each Profile_ID contributes to that SNAP, via the NFRs the Profile_ID relates to, *first*. Then weight the raw data accordingly (dependent on the Profile_ID it belongs to), and then produce the sector GAM for the SNAP on the fly. I think there would be a much bigger processing overhead here as the sector GAMs are created from raw data every time.

**IDEA to develop (28/01/22):** Create one large table of emissions for all years and pollutants. When creating the sector GAM, use the mean grouped-NFR contribution to that sector across *all years and all pollutants*. That way it becomes generic. This could be a simple mean (total per sector would need to be re-adjusted to 1), a weighted mean to give more recent years more influence (as old stuff is out of date?) or another GAM. It could be fine that one % contribution figure is used to make the GAM or the range across all years/pollutants could feed into a MCMC?

### To Do:
* convert emissions to emissions per hour over a year. 
* have added Sentinel5p comparison script - need to adjust functions
* Test integration of year-specific data into generic profiles (e.g. coal-energy into larger SNAP 1)

**Future DUKEMs work:**\
    * Integrate wood burning hour of day data; Gary Fuller (email 25/10/2021)\
    * continue to add data & detail around sub-sector profiles\
    * Pollutant specific profiles\
    * Agricultural information re activity AND temperature/climate related profiles\
    * More year-specific data\
    * Spatially variable temporal profiles\
    * Point data temporal profiles\
    * (near) Real time data interface

