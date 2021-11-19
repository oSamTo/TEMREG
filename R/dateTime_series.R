require(lubridate)
require(data.table)

## create time data from 2000 and attach bank holiday information

bnakhols <- fread("C:/FastProcessingSam/DUKEMS/Database/ukbankholidays.csv", header=T)
bnakhols[,V2:=NULL]
bnakhols[, Date := as_date(`UK BANK HOLIDAYS`, format="%d-%b-%Y")]
bnakhols[,pub_hol := T]
bnakhols <- bnakhols[Date >= "2000-01-01"]
bnakhols[,`UK BANK HOLIDAYS` := NULL]

dt_time <- data.table(Date=seq(as.Date("2000-01-01"),as.Date("2021-12-31"),by="days"))

dt_time[, year_int := year(Date)]
dt_time[, mon_int := month(Date)]
dt_time[, mon_char := month.abb[mon_int]]
dt_time[, woy_int := week(Date)]
dt_time[, dow_int := wday(Date)]
dt_time[, dow_char := c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")[dow_int]]
dt_time[, jul_day := yday(Date)]
dt_time[, timezone := "UTC"]

dt_with_hols <- bnakhols[dt_time, on="Date"]
dt_with_hols[is.na(pub_hol), pub_hol := F]

fwrite(dt_with_hols,"C:/FastProcessingSam/DUKEMS/Database/temporal.csv")


