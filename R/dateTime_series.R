source("./R/dukem.R")

## create time data from 2000 and attach bank holiday information

bankhols <- fread("./Data/ukbankholidays.csv")
bankhols[, Date := lubridate::dmy(BANK_HOLIDAYS)]
bankhols[,pub_hol := T]
bankhols <- bankhols[Date >= "2000-01-01"]
bankhols[,BANK_HOLIDAYS :=  NULL]

dt_time <- data.table(Date=seq(as.Date("2000-01-01"),as.Date("2023-12-31"),by="days"))

dt_time[, year := year(Date)]
dt_time[, month := month(Date)]
#dt_time[, month_char := month.abb[month_int]]
dt_time[, week := week(Date)]
dt_time[, wday := wday(Date, week_start = getOption("lubridate.week.start", 1))]
#dt_time[, wday_char := c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")[wday_int]]
dt_time[, yday := yday(Date)]
dt_time[, timezone := "UTC"]

dt_with_hols <- bankhols[dt_time, on="Date"]
dt_with_hols[is.na(pub_hol), pub_hol := F]

fwrite(dt_with_hols,"./Data/Dates.csv")


