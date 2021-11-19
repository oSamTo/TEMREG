source("./R/workspace.R")

## create time data from 2000 and attach bank holiday information

bankhols <- fread("./Data/ukbankholidays.csv", header=T)
bankhols[,V2:=NULL]
bankhols[, Date := as_date(`UK BANK HOLIDAYS`, format="%d-%b-%Y")]
bankhols[,pub_hol := T]
bankhols <- bankhols[Date >= "2000-01-01"]
bankhols[,`UK BANK HOLIDAYS` := NULL]

dt_time <- data.table(Date=seq(as.Date("2000-01-01"),as.Date("2021-12-31"),by="days"))

dt_time[, year_int := year(Date)]
dt_time[, mon_int := month(Date)]
dt_time[, mon_char := month.abb[mon_int]]
dt_time[, woy_int := week(Date)]
dt_time[, dow_int := wday(Date)]
dt_time[, dow_char := c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")[dow_int]]
dt_time[, jul_day := yday(Date)]
dt_time[, timezone := "UTC"]

dt_with_hols <- bankhols[dt_time, on="Date"]
dt_with_hols[is.na(pub_hol), pub_hol := F]

fwrite(dt_with_hols,"./Data/Dates.csv")


