# DUKEMs_TP

######################################################################
#### **Processing UK annual emissions to sector temporal profiles.**
######################################################################

*Info:*
----------------

* Temporal profiles created for ~90 (currently) sub-sectors (see DUKEMs project work) - "Profile_IDs"
* These are produced for hour, hourwday (hour for specific day of week), wday, month and yday

* Functions here will;
    * Join NFR NAEI emissions to DUKEMs Profile_IDs
    * For a Sector Classification, aggregate emissions by Profile_ID and create new temporal profiles
    * Plot emissions and write some metadata

* Currently produces profiles for SNAP or GNFR classes (16/01/2022 : GNFR not ready)

* run from ./R/profile_emissions.R

* coefficients are produced to 1

## currently: NOx, SOx, CH4, CO2, N2O

**Future work:**\
    * Integrate wood burning hour of day data; Gary Fuller (email 25/10/2021)\
    * continue to add data & detail around sub-sector profiles\
    * Pollutant specific profiles\
    * Agricultural information re activity AND temperature/climate related profiles\
    * More year-specific data\
    * HDD year-specific\
    * DOY patterns as option?\
    * Spatially variable temporal profiles\
    * Point data temporal profiles\
    * (near) Real time data interface

