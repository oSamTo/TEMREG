packs <- c("sp","raster","stringr","gdalUtils","rgeos","rgdal","grid","plyr","car","reshape2","ggplot2","ggrepel","data.table","stats","readr","ggplot2","sf","rasterVis","ncdf4","readxl","lubridate")
lapply(packs, require, character.only = TRUE)

setwd("C:/FastProcessingSam/DUKEMS")

## read in many NAEI emissions files to get the full picture of Source/Activity sectors. Match to SNAP etc. 


all_dirs <- list.dirs(path = "//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/NAEI_data_and_SNAPS/NAEI_data/diffuse",  full.names = TRUE, recursive = TRUE)

sub_dirs <- all_dirs[grep("time_series", all_dirs)]

ems_files <- list.files(sub_dirs, recursive = T, pattern = "-2019.csv$", full.names = T)

l_dt <- lapply(ems_files, fread, header=T)

l_dt_sectors <- lapply(l_dt, "[", , c("Gas","NFR/CRF Group","Source","Activity"))

dt_sectors <- rbindlist(l_dt_sectors, use.names = T)
dt_sectors <- dt_sectors[Source != ""]

dt_sectors_uni <- dt_sectors[!duplicated(dt_sectors)]

setnames(dt_sectors_uni, "NFR/CRF Group", "NFR19")

#dt_sectors_uni[grepl("Field burn", Source)]
#dt_sectors_uni[Source=="Digestate"]
#dt_sectors_uni[grepl("managed as", Activity)]
## lookup

dt_lu <- fread("//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/NAEI_data_and_SNAPS/lookups/NFR19_SNAP_lookup.csv")
dt_lu <- dt_lu[,c("NFR19","SNAP","GNFR19")]
dt_lu <- dt_lu[NFR19 != "" ]
dt_lu <- dt_lu[GNFR19 != "" ]
dt_lu <- dt_lu[!duplicated(dt_lu)]

gferf <- dt_lu[dt_sectors_uni, on= "NFR19"]


fwrite(gferf,"./Database/Sectors.csv")




