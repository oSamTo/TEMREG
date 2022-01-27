# DUKEMs_TP

######################################################################
#### **Processing UK annual emissions to sector-level temporal profiles.**
######################################################################

*Info:*
----------------

* Raw data collected for activity data in the UK, processed to standardized output 
    * Processing of data in C:/...../DUKEMS/Database/R/xxxxxxx.R
    * Standard input files uploaded to ./data/GAM_input  (not domestic, too large)
* Each NFR code maps to a small grouping (Profile_ID), which maps to a GAM_input data category. 
* GAMs are produced from the raw data (./R/xxxxx.R), grouped by Profile_ID, and placed in ./data/GAM_output
* They are produced for hour, hourwday (hour for specific day of week), wday, month and yday
* These GAM outputs are the basis for sector-level groupings;
    * Need to know the NFR to sector look-up (e.g. to SNAP)
    * The NFRs are grouped by Profile_ID and the total emissions represented by the Profile_ID, as a proportion of sector total, defines a Profile weight
    * The weights influence the construction of the final sector level GAM (./output/GAM_sector)
* As well as GAMs, standard csv outputs of temporal profiles are produced; coefficients are relative to 1

* When further raw data is found and ingested, it must either map to current Profile_IDs or new IDs need to be added and mapped to NFR codes. 

----------------

**Emissions contribution by Profile_ID per sector here; ./doc/profile_contributions_sector/**

Re GAM of GAMs: There is an example script (./R/weighted_gams.R) of taking a GAM by group from raw data and then producing a GAM from those two GAMs, with a weighting assigned to the two GAMs produced from raw data.

This is to mimic emissions; the raw data groups represent emissions grouped by NFR code choices (Profile_IDs) and the weights on those GAMs for the next level GAM represents the overall emissions contribution to a sector (e.g. SNAP). Therefore you end up taking into account the relative contribution of a sub-sector in each main sector.

**IDEA (19/01/22):** The above could be changed by giving the original raw data a weighting based on a lookup table. For example, at the moment there are ~90 Profile_IDs across all SNAPs. Going through each SNAP, you calculate what each Profile_ID contributes to that SNAP *first*, weight the raw data accordingly (dependent on the Profile_ID it belongs to), and then produce the sector GAM for the SNAP on the fly. This allows for **any** combination of NFR to Profile_ID to be read in but requires every NFR to be linked to a datasource and for a separate NFR to Profile_ID lookup table. There would be a much bigger processing overhead here.

This totally rewrites what I have done and would require some effort. Discuss. 

### Currently: NOx, SOx, CH4, CO2, N2O

### To Do:
* have added Sentinel5p comparison script - need to adjust functions
* Talk about the method for sector GAM outlined above
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

