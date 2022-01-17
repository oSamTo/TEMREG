# DUKEMs_TP

########################################################################
#### **Processing UK annual emissions to hourly temporal profiles.**
########################################################################

*Info:*
----------------

* Currently produces profiles for SNAP or GNFR classes
* Currently produces separate hour of day / day of week / month of year tables for each sector
* Option to produce hourly emissions table (`hour_emis`)

* run from ./R/profile_emissions.R

* coefficients are produced to 1

## currently: NOx, SOx, CH4, CO2, N2O

**Future work:**\
    * Integrate wood burning hour of day data; Gary Fuller (email 25/10/2021)\
    * Increase data & detail around sub-sector profiles\
    * Pollutant specific profiles\
    * Agricultural information re activity AND temperature/climate related profiles\
    * More year-specific data\
    * HDD year-specific\
    * DOY patterns as option?\
    * Spatially variable temporal profiles\
    * Point data temporal profiles\
    * (near) Real time data interface

