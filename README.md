# DUKEMs_TP

######################################################################
#### **Processing UK annual emissions to sector-level temporal profiles.**
######################################################################

*Info:*
----------------

* Temporal profiles created for ~90 sub-sectors (currently; see DUKEMs project work) - "Profile_IDs"
* These are produced for hour, hourwday (hour for specific day of week), wday, month and yday

* Functions here will;
    * Join NFR NAEI emissions to DUKEMs Profile_IDs
    * For a Sector Classification, aggregate emissions by Profile_ID and create new sector temporal profiles
    * Plot emissions and write some metadata

* Currently produces profiles for SNAP or GNFR classes (16/01/2022 : GNFR not ready)

* run from ./R/profile_emissions.R

* coefficients are produced relative to 1


**Emissions contribution by Profile_ID per sector here; ./doc/profile_contributions_sector/**

Re GAM of GAMs: There is an example script (./R/weighted_gams.R) of taking a GAM by group from raw data and then producing a GAM from those two GAMs, with a weighting assigned to the two GAMs produced from raw data.

This is to mimic emissions; the raw data groups represent emissions grouped by NFR code choices (Profile_IDs) and the weights on those GAMs for the next level GAM represents the overall emissions contribution to a sector (e.g. SNAP). Therefore you end up taking into account the relative contribution of a sub-sector in each main sector.

**IDEA (19/01/22):** The above could be changed by giving the original raw data a weighting based on a lookup table. For example, at the moment there are ~90 Profile_IDs across all SNAPs. Going through each SNAP, you calculate what each Profile_ID contributes to that SNAP *first*, weight the raw data accordingly (dependent on the Profile_ID it belongs to), and then produce the sector GAM for the SNAP on the fly. This allows for **any** combination of NFR to Profile_ID to be read in but requires every NFR to be linked to a datasource and for a separate NFR to Profile_ID lookup table. 

This totally rewrites what I have done below. 

## currently: NOx, SOx, CH4, CO2, N2O

**Future DUKEMs work:**\
    * Integrate wood burning hour of day data; Gary Fuller (email 25/10/2021)\
    * continue to add data & detail around sub-sector profiles\
    * Pollutant specific profiles\
    * Agricultural information re activity AND temperature/climate related profiles\
    * More year-specific data\
    * Test integration of year-specific data into generic profiles (e.g. coal-energy into larger SNAP 1)
    * Spatially variable temporal profiles\
    * Point data temporal profiles\
    * (near) Real time data interface

