#### Vincent Fugere, 2018 - 2020
#### LEAP 2017 experiment
#### Code to format data, to be used in Rmd file

rm(list=ls())

library(tidyverse)
library(scales)
source('/Users/vincentfugere/Google Drive/Recherche/PhD/R/functions/utils.R')
library(scales)
library(mgcv)
library(itsadug)
library(magrittr)
library(vegan)
library(readxl)
library(chron)
library(brms)

cols<-c('firebrick2','gold2','forestgreen','darkblue')
treat <- read.csv('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/LEAP2017treatments.csv', stringsAsFactors = F)
to.rm <- c('S1','S2','S3','S4','LAKE','P2C1','P2C2','P2C3','P2C4') #useless for this project

#relevant experimental period (week 0 to week 7)
phase1.dates <- as.Date(c('06-06-2017','07-27-2017'),'%m-%d-%Y')
phase1.dates.julian <- as.numeric(format(phase1.dates, '%j'))
date.range <- 157:208
rm(phase1.dates.julian,phase1.dates)

#### Load and format data ####

## depth

depth <- read.csv('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/depth.csv') %>%
  select(month:depth) %>% unite(pond, sub.array, pond.number, sep='') %>%
  unite(date, month, day, sep='_') %>% filter(pond %!in% to.rm)
depth$date <- depth$date %>% paste0(.,'_17') %>% as.Date(., format = '%m_%d_%y')
depth$date <- as.numeric(format(depth$date, '%j'))
depth %<>% filter(date %in% date.range) %>% arrange(date,pond)
depth$week <- rep(0:7,each=96)
depth %<>% select(pond,week,depth)

## YSI data

ysi <- read.csv('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/ysi.csv') %>%
  select(month:ph.after) %>% unite(pond, sub.array, pond.number, sep='') %>% unite(date, month, day, sep='_') %>%
  filter(pond %!in% to.rm) %>% filter(pond != 'NA')
ysi$date <- ysi$date %>% paste0(.,'_17') %>% as.Date(., format = '%m_%d_%y')
ysi %<>% filter(date != '2017-06-14')
ysi$date <- as.numeric(format(ysi$date, '%j'))
ysi %<>% filter(date %in% date.range) %>% arrange(date,pond)
ysi$week <- rep(0:7,each=96)
ysi %<>% select(pond, week, SPC:ph.after)

## phytoplankton

tp0 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Jun07.txt', skip=2, stringsAsFactors = F)
tp1 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Jun14.txt', skip=2, stringsAsFactors = F)
tp2 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Jun21.txt', skip=2, stringsAsFactors = F)
tp3 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Jun28.txt', skip=2, stringsAsFactors = F)
tp4 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Jul05.txt', skip=2, stringsAsFactors = F)
tp5 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Jul12.txt', skip=2, stringsAsFactors = F)
tp6 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Jul19.txt', skip=2, stringsAsFactors = F)
tp7 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Jul26.txt', skip=2, stringsAsFactors = F)
tp8 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Aug02.txt', skip=2, stringsAsFactors = F)
tp9 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Aug09.txt', skip=2, stringsAsFactors = F)
tp10 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Aug16.txt', skip=2, stringsAsFactors = F)
tp11 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Aug30.txt', skip=2, stringsAsFactors = F)
tp12 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/fluoroprobe/Sep25.txt', skip=2, stringsAsFactors = F)
phyto <- rbind(tp0,tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9,tp10,tp11,tp12)
rm(tp0,tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9,tp10,tp11,tp12)

colnames(phyto)[1:10] <- c('date','time1','time2','pond','greens','cyanos','diatoms','cryptos','ys','total')
phyto <- phyto %>% unite(time, time1, time2, remove=T) %>% select(date:total) %>% select(-ys) %>%
  filter(pond %!in% to.rm) %>% group_by(date, pond) %>% filter(nchar(pond) == 2) %>%
  summarise_if(is.numeric, mean) %>% ungroup
phyto$date <- as.Date(phyto$date, format = '%d/%m/%Y')
phyto$date <- as.numeric(format(phyto$date, '%j'))
phyto %<>% filter(date %in% date.range) %>% arrange(date,pond)
phyto$week <- rep(0:7,each=96)
phyto %<>% select(pond, week, greens:total)

## zooplankton

zoops_tot <- read.csv('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/zooplankton/total_counts_undergrads.csv', stringsAsFactors = F)
zoops_tot$Chidorus[is.na(zoops_tot$Chidorus)] <- 0
zoops_tot <- zoops_tot %>% rename('sample' = Date) %>%
  filter(Pond %!in% to.rm) %>%
  filter(Pond != 'NA') %>%
  mutate('Clad.perL' = (Cladocerans + Chidorus)/2, 'Cop.perL' = Copepods/2) %>%
  select(sample, Clad.perL, Cop.perL, density.indperL) %>%
  distinct(sample, .keep_all=T) %>%
  separate(sample, c('pond','day','month')) %>%
  add_column(year = 2017) %>%
  unite(date, day, month, year, sep='_') %>%
  rename('zd' = density.indperL)
zoops_tot$date <- as.Date(zoops_tot$date, format = '%d_%m_%Y')

zoops_com <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/zooplankton/community_composition_Lynne.xlsx') %>% 
  select(pond:density.indperL) %>% as.data.frame
zoops_com[is.na(zoops_com)] <- 0
zoops2 <- zoops_com %>% filter(pond %!in% to.rm) %>%
  filter(pond != 'NA') %>%
  mutate('Clad.perL' = Cladocerans/2, 'Cop.perL' = Copepods/2) %>%
  select(pond,date,Clad.perL,Cop.perL,density.indperL) %>%
  rename('zd' = density.indperL) %>% 
  filter(date != '08.06.17')
zoops2$date <- as.Date(zoops2$date, format = '%d.%m.%y')
zoo <- bind_rows(zoops_tot,zoops2)
rm(zoops2, zoops_tot)

zoo$date <- as.numeric(format(zoo$date, '%j'))
zoo %<>% filter(date %in% date.range) %>% arrange(date,pond)
zoo$week <- rep(0:7,each=96)
zoo %<>% select(pond, week, Clad.perL:zd)

zoops_com <- filter(zoops_com, date == '27.07.17', pond %!in% to.rm) %>% 
  select(pond, `Bosmina longirostris`:copepodids)

## metabolism

metab <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/LEAP/leap2017/data/LEAP2017_metabolism.xlsx') %>%
  unite(pond, sub.array, pond.number, sep='', remove=T) %>%
  unite(date, month.day1, date.day1, sep='_') %>%
  filter(pond %!in% to.rm) %>% filter(pond != 'NA')

metab$date <- as.Date(metab$date, format = '%m_%d')
metab$date <- as.numeric(format(metab$date, '%j'))
metab %<>% filter(date %in% date.range) %>% arrange(date,pond)
metab$week <- rep(c(0,2,4,6,7),each=96)

#calculating the duration (in hours) of nightime and daytime periods
time1 <- strftime(metab$sunset.day1, format = "%H:%M:%S", tz='UTC')
time2 <- strftime(metab$sunrise.day2, format = "%H:%M:%S", tz='UTC')
time3 <- strftime(metab$sunset.day2, format = "%H:%M:%S", tz='UTC')
day1 <- rep('01/01/2019', length(time1))
day2 <- rep('01/02/2019', length(time1))
x1 <- chron(dates = day1, times = time1)
x2 <- chron(dates = day2, times = time2)
x3 <- chron(dates = day2, times = time3)
period1 <- as.numeric(difftime(x2, x1, unit="hours"))
period2 <- as.numeric(difftime(x3, x2, unit="hours"))
metab$period1 <- period1
metab$period2 <- period2
rm(time1,time2,time3,day1,day2,x1,x2,x3,period1,period2)

#converting mg to DO surplus or deficit (compared to mg at 100% sat) 
#to account for large temp changes during the night. ER and NEP computed from these.
metab$DO.mg.1 <- metab$DO.mg.1-(metab$DO.mg.1*(100/metab$DO.prctsat.1)) #term in parenthesis is mg/L at 100% sat
metab$DO.mg.2 <- metab$DO.mg.2-(metab$DO.mg.2*(100/metab$DO.prctsat.2)) #term in parenthesis is mg/L at 100% sat
metab$DO.mg.3 <- metab$DO.mg.3-(metab$DO.mg.3*(100/metab$DO.prctsat.3)) #term in parenthesis is mg/L at 100% sat

#adding nightime respiration and daytime NEP, both in mg O2 / L / hr
metab$ER <- with(metab, (DO.mg.2-DO.mg.1)/period1)
metab$NEP <- with(metab, (DO.mg.3-DO.mg.2)/period2)

# hist(metab$ER) #some positive values
# hist(metab$NEP) #some negative values

#excel sheet suggests data for 06-20 (week 2) might be garbage due to a probe problem
metab <- filter(metab, week != 2)

# hist(metab$ER)
# hist(metab$NEP)
# #yes that fixes it!

metab <- select(metab, pond, week, ER, NEP)

############################################################
## MISSING FOR PAPER: PERIPHYTON, NUTRIENTS, WATER COLOUR ##
############################################################

## binding and adding treatments

data <- left_join(depth, ysi, by = c('pond','week')) %>% 
  left_join(phyto, by = c('pond','week')) %>% 
  left_join(zoo, by = c('pond','week')) %>% 
  left_join(metab, by = c('pond','week'))

data <- treat %>% select(pond.ID:pH.var) %>%
  right_join(data, by = c('pond.ID' = 'pond'))

data <- select(data, pond.ID,week,MC.ID:pH.var,depth:NEP)

colnames(data)[c(1,3,4,5,7,8)] <- c('pond','MC','disp','pH.trt','MC.pH','str')


#### Optional filters ####

#getting rid of two ponds which were accidently acidified to a lower pH than intended

#getting rid of homogeneous metacommunities for now for simplicity (fix later)

#adding measured pH

data$pH <- data$ph.after
data$pH[1:96] <- data$ph.before[1:96]
data <- select(data, -ph.before, -ph.after)
data <- select(data, pond, week, pH, everything())

data <- data %>% mutate_at(vars(Clad.perL,Cop.perL,zd,greens,diatoms,total), list(log = ~log1p(.)))

#end of data formatting