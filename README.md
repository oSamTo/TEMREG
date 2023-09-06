# TEMREG: Temporally Resolved Emissions Generator

######################################################################
#### **Processing UK annual emissions to sector-level temporal profiles.**
######################################################################

# Notes on structural change: 
* Code for timescale GAM with `by=` grouping for sectors is ready, but very computationally expensive. Takes a long time. 
* Therefore code is also available to construct sector GAMs one by one and place into list. List is saved as .rds

* All of this to go in targets?
    * if you updated the raw data GAMs, you need to update NFR_to_GAM table. Targets will spot this and re-run
    * if you update the pollutants to work over, targets should spot and re-run


*Info:*
----------------

* UK activity data is used - not emissions data - to match to activity sectors in the NAEI
* Raw data processed to standardized output 
    * Processing of data in /DUKEMS/Database/R/gam_generation.R
    * Standard input files uploaded to ./data/input_formatted  (all on .gitignore)
* Each NFR code translates to an activity grouping (Profile_ID), and the grouping's total emissions will determine the its weight in the final aggregated sector.
* GAMs for every Profile_ID are produced from the raw data (data/GAM_profileID/*timescale*)
* They are produced for hour, hourwday (hour for specific day of week), wday, month and yday
* These GAM outputs are the basis for sector-level groupings;
    * Need to know the NFR to sector look-up (e.g. to SNAP, GNFR, custom)
    * The weights (contribution of Profile_ID to aggregated sector) influence the construction of the final sector level GAM (./output/GAM_sector)
    * Final sector level GAMs are output (./output/GAM_sector)
* As well as GAMs, standard csv outputs of temporal profiles are produced; coefficients are relative to 1

* GAMs for sectors can be constructed from year/pollutant averages or year/pollutant specific.  
    * This allows for longer/shorter-term representations, multiple or single pollutants

* When further raw data is found and ingested, it must either map to current Profile_IDs or new IDs need to be added and mapped to NFR codes. 

----------------

**Emissions contribution by Profile_ID per sector here; ./doc/profile_contributions_sector/**

**IMPORTANT:** Raw data should be normalised (units removed) to enable combination of disparate data.

**Aggregating GAMs** When aggregating Profile_ID GAMs to aggregated GAMs, each Profile_ID GAM is sampled across the prediction interval (n times). This uses the uncertainty range properly, instead of sampling the SE (which creates false confidence).   

**IDEA to develop (28/01/22):** Create one large table of emissions for all years and pollutants. When creating the sector GAM, use the mean grouped-NFR contribution to that sector across *all years and all pollutants*. That way it becomes generic. This could be a simple mean (total per sector would need to be re-adjusted to 1), a weighted mean to give more recent years more influence (as old stuff is out of date?) or another GAM. It could be fine that one % contribution figure is used to make the GAM or the range across all years/pollutants could feed into a MCMC?

**Future DUKEMs work:**\
    * Integrate wood burning hour of day data; Gary Fuller (email 25/10/2021)\
    * continue to add data & detail around sub-sector profiles\
    * Pollutant specific profiles\
    * Agricultural information re activity AND temperature/climate related profiles\
    * More year-specific data\
    * Spatially variable temporal profiles\
    * Point data temporal profiles\
    * (near) Real time data interface

