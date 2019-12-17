#### Vincent Fugere, 2018 - 2019
#### LEAP 2017 experiment
#### Code to examine effect of dispersal and pH heterogeneity on zooplankton metacommunities

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
treat <- read.csv('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/LEAP2017treatments.csv', stringsAsFactors = F)
to.rm <- c('S1','S2','S3','S4','LAKE','P2C1','P2C2','P2C3','P2C4') #useless for this project

#relevant experimental period (week 0 to week 7)
date.range <- 157:208

#### Load and format data ####

## depth

depth <- read.csv('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/depth.csv') %>%
  select(month:depth) %>% unite(pond, sub.array, pond.number, sep='') %>%
  unite(date, month, day, sep='_') %>% filter(pond %!in% to.rm)
depth$date <- as.Date(depth$date, format = '%m_%d')
depth$date <- as.numeric(format(depth$date, '%j'))
depth %<>% filter(date %in% date.range) %>% arrange(date,pond)
depth$week <- rep(0:7,each=96)
depth %<>% select(pond,week,depth)

## YSI data

ysi <- read.csv('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/ysi.csv') %>%
  select(month:ph.after) %>% unite(pond, sub.array, pond.number, sep='') %>% unite(date, month, day, sep='_') %>%
  filter(pond %!in% to.rm) %>% filter(pond != 'NA')
ysi$date <- as.Date(ysi$date, format = '%m_%d')
ysi %<>% filter(date != '2019-06-14')
ysi$date <- as.numeric(format(ysi$date, '%j'))
ysi %<>% filter(date %in% date.range) %>% arrange(date,pond)
ysi$week <- rep(0:7,each=96)
ysi %<>% select(pond, week, SPC:ph.after)

## phytoplankton

tp0 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Jun07.txt', skip=2, stringsAsFactors = F)
tp1 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Jun14.txt', skip=2, stringsAsFactors = F)
tp2 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Jun21.txt', skip=2, stringsAsFactors = F)
tp3 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Jun28.txt', skip=2, stringsAsFactors = F)
tp4 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Jul05.txt', skip=2, stringsAsFactors = F)
tp5 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Jul12.txt', skip=2, stringsAsFactors = F)
tp6 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Jul19.txt', skip=2, stringsAsFactors = F)
tp7 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Jul26.txt', skip=2, stringsAsFactors = F)
tp8 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Aug02.txt', skip=2, stringsAsFactors = F)
tp9 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Aug09.txt', skip=2, stringsAsFactors = F)
tp10 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Aug16.txt', skip=2, stringsAsFactors = F)
tp11 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Aug30.txt', skip=2, stringsAsFactors = F)
tp12 <- read.table('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/fluoroprobe/Sep25.txt', skip=2, stringsAsFactors = F)
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

zoops_tot <- read.csv('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/zoops.csv', stringsAsFactors = F)
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

zoops_com <- read.csv('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/zoops_community.csv', stringsAsFactors = F) 
zoops2 <- zoops_com %>% filter(pond %!in% to.rm) %>%
  filter(pond != 'NA') %>%
  mutate('Clad.perL' = (Bosmina.longispina + Diaphanosoma.sp + Sida.crystallina + Chydorus.sphaericus + Daphnia.pulex + Simocephalus.sp + Daphnia.ambigua + Ceriodaphnia.dubia)/2, 'Cop.perL' = Cyclops.scutifer/2) %>%
  select(pond,date,Clad.perL,Cop.perL,density) %>%
  rename('zd' = density)
zoops2$date <- as.Date(zoops2$date, format = '%d.%m.%y')
zoo <- bind_rows(zoops_tot,zoops2)
rm(zoops2, zoops_tot)

zoo$date <- as.numeric(format(zoo$date, '%j'))
zoo %<>% filter(date %in% date.range) %>% arrange(date,pond)
zoo$week <- rep(0:7,each=96)
zoo %<>% select(pond, week, Clad.perL:zd)

zoops_com <- filter(zoops_com, date == '27.07.17', pond %!in% to.rm) %>% 
  select(pond, Chaoborus:Ceriodaphnia.dubia)

## metabolism

metab <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/LEAP2017_metabolism.xlsx') %>%
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

#### Exploring pH, disp, and effects on local communities ####

data$pH <- data$ph.after
data$pH[1:96] <- data$ph.before[1:96]
data <- select(data, -ph.before, -ph.after)
data <- select(data, pond, week, pH, everything())

data <- data %>% mutate_at(vars(Clad.perL,Cop.perL,zd,greens,diatoms,total), list(log = ~log1p(.)))

data <- filter(data, week < 7)

homo <- filter(data, str == 'homogeneous', week > 0)
hete <- filter(data, str == 'heterogeneous', week > 0)

colfunc <- colorRampPalette(cols)
colfunc(100) -> cols.plot

#####

pdf('~/Desktop/explo.pdf',width=8.5,height=11,onefile = T)
par(mfrow=c(2,1),cex=1)

for(i in 10:27){

  real.var.name <- colnames(data)[i]

  tempdata <- homo[,c(1,2,3,5,6,i)] %>% drop_na %>% as.data.frame
  colnames(tempdata)[6] <- 'y'
  tempdata$disp <- as.factor(tempdata$disp)
  tempdata$pH.trt <- as.factor(tempdata$pH.trt)
  tempdata$pond <- as.factor(tempdata$pond)
  # gam.homo <- bam(y ~ disp + s(week, k = 4) + s(pH, k = 4) + ti(week, pH, k=4) + s(week, by = disp, k=4) + s(week, pond, bs='fs',m=1,k=4), data=tempdata, method='fREML')
  # summary(gam.homo)
  tempdata$ph.col.idx <- round(rescale(tempdata$pH, c(1,100)))
  plot(y~week,tempdata,type='n',yaxt='n',xaxt='n',xlim=c(1,6),ylim=range(tempdata$y),ann=F,bty='l')
  title(ylab=real.var.name)
  title(xlab="weeks of acid treatment")
  axis(2,lwd=0,lwd.ticks=1)
  axis(1,lwd=0,lwd.ticks=1)
  for(p in 1:48){
    data.tmp <- subset(tempdata, pond == levels(tempdata$pond)[p])
    points(y ~ week, type='l', data.tmp, pch=16, cex=0.5, lwd=0.2, col = alpha(cols.plot[data.tmp$ph.col.idx],0.5))
  }
  means <- aggregate(y ~ disp*pH.trt*week, tempdata, FUN = 'mean')
  for(pH in 1:4){
    for(d in 1:3){
      data.tmp <- filter(means, disp == levels(means$disp)[d], pH.trt == levels(means$pH.trt)[pH])
      plot.data <- as.data.frame(approx(x = data.tmp$week, y = data.tmp$y, n = 1000))
      lines(y ~ x, plot.data, lwd=3, col = alpha(cols[pH],0.8), lty=c(1,2,3)[d])
    }
  }
  rm(data.tmp,p,pH,d,plot.data)
  title(main = 'homo', cex = 0.5)
  # legend('topright',bty='n',legend = '(a)',inset = c(0,-0.05), cex =1.2)
  # legend('bottomright',bty='n',legend = c('pH 8.5','pH 7','pH 5.5','pH 4'),pch=16,col=cols[4:1],y.intersp=1.2)

  tempdata <- hete[,c(1,2,3,5,6,i)] %>% drop_na %>% as.data.frame
  colnames(tempdata)[6] <- 'y'
  tempdata$disp <- as.factor(tempdata$disp)
  tempdata$pH.trt <- as.factor(tempdata$pH.trt)
  tempdata$pond <- as.factor(tempdata$pond)
  tempdata$ph.col.idx <- round(rescale(tempdata$pH, c(1,100)))
  plot(y~week,tempdata,type='n',yaxt='n',xaxt='n',xlim=c(1,6),ylim=range(tempdata$y),ann=F,bty='l')
  title(ylab=real.var.name)
  title(xlab="weeks of acid treatment")
  axis(2,lwd=0,lwd.ticks=1)
  axis(1,lwd=0,lwd.ticks=1)
  for(p in 1:48){
    data.tmp <- subset(tempdata, pond == levels(tempdata$pond)[p])
    points(y ~ week, type='l', data.tmp, pch=16, cex=0.5, lwd=0.2, col = alpha(cols.plot[data.tmp$ph.col.idx],0.5))
  }
  means <- aggregate(y ~ disp*pH.trt*week, tempdata, FUN = 'mean')
  for(pH in 1:4){
    for(d in 1:3){
      data.tmp <- filter(means, disp == levels(means$disp)[d], pH.trt == levels(means$pH.trt)[pH])
      plot.data <- as.data.frame(approx(x = data.tmp$week, y = data.tmp$y, n = 1000))
      lines(y ~ x, plot.data, lwd=3, col = alpha(cols[pH],0.8), lty=c(1,2,3)[d])
    }
  }
  rm(data.tmp,p,pH,d,plot.data)
  title(main = 'hetero', cex = 0.5)
}

dev.off()

#####

# #pdf('~/Desktop/hists.pdf',width=8,height=10.5,pointsize = 8)
# par(mfrow=c(4,5),cex=1)
# for(i in 10:27){
#   hist(data[,i],main=NULL,xlab=colnames(data[i]),breaks=30)
# }
# #dev.off()

##### 

# MODELS

homo <- mutate_at(homo, vars(pond,disp,pH.trt,MC), list(f = ~as.factor(.)))
hete <- mutate_at(hete, vars(pond,disp,pH.trt,MC), list(f = ~as.factor(.)))
data <- mutate_at(data, vars(pond,disp,pH.trt,MC,str), list(f = ~as.factor(.)))

# library(lme4)
# m1 <- lmer(zd_log ~ pH.trt_f + disp_f + str_f + pH.trt_f:disp_f + disp_f:str_f + pH.trt_f:disp_f:str_f + (1|MC_f/pond_f), data)
# #summary(m1)
# anova(m1)
# plot(m1)
# 
# dd <- hete
# dd$var <- dd$total_log
# m1 <- update(m1, newdata = dd)
# #m1 <- brm(var ~ pH.trt_f * disp_f + (1|MC_f/pond_f), dd, cores=4)
# summary(m1)
# plot(m1)
# pp_check(m1)
# marginal_effects(m1)
# 
# m1 -> m.total
# 
#save(m.total,m.diatoms,m.greens,m.zd,m.cop,m.clad,m.NEP,m.ER, file = '~/Desktop/LMMs_6_hete.RData')
# save(chla_hete,chla_homo,clad_hete,clad_homo,cop_hete,cop_homo,depth_hete,depth_homo,diatoms_hete,diatoms_homo,ER_hete,ER_homo,greens_hete,greens_homo,NEP_hete,NEP_homo,SPC_hete,SPC_homo,zd_hete,zd_homo, file = '~/Desktop/LMMs.RData')
# rm(dd,m1,means,plot,tempdata,i,real.var.name)

#### stability ####

stability <- data.frame()

for(i in 1:96){
  tmp <- filter(data, pond == levels(data$pond_f)[i])
  greens <- mean(tmp$greens_log)/sd(tmp$greens_log)
  diatoms <- mean(tmp$diatoms_log)/sd(tmp$diatoms_log)
  total <- mean(tmp$total_log)/sd(tmp$total_log)
  clad <- mean(tmp$Clad.perL_log)/sd(tmp$Clad.perL_log)
  cop <- mean(tmp$Cop.perL_log)/sd(tmp$Cop.perL_log)
  zoo <- mean(tmp$zd_log)/sd(tmp$zd_log)
  results <- data.frame(c(tmp[1,c('pond','MC','disp','pH.trt','str')],'greens'=greens,'diatoms'=diatoms,'total'=total,'clad'=clad,'cop'=cop,'zoo'=zoo))
  stability <- rbind(stability,results)
}

stability <- filter(stability, str == 'heterogeneous')

pdf('~/Desktop/stability.pdf',width=6,height = 4,onefile = T,pointsize = 8)

boxplot(greens ~ disp * pH.trt, stability, ylab = 'greens')
anova(lmer(greens ~ disp * pH.trt + (1|MC), data = stability))

boxplot(diatoms ~ disp * pH.trt, stability, ylab = 'diatoms')
anova(lmer(diatoms ~ disp * pH.trt + (1|MC), data = stability))

boxplot(total ~ disp * pH.trt, stability, ylab = 'total')
anova(lmer(total ~ disp * pH.trt + (1|MC), data = stability))

boxplot(clad ~ disp * pH.trt, stability, ylab = 'clad')
anova(lmer(clad ~ disp * pH.trt + (1|MC), data = stability))

boxplot(cop ~ disp * pH.trt, stability, ylab = 'cop')
anova(lmer(cop ~ disp * pH.trt + (1|MC), data = stability))

boxplot(zoo ~ disp * pH.trt, stability, ylab = 'zoo')
anova(lmer(zoo ~ disp * pH.trt + (1|MC), data = stability))

dev.off()

#### zoops composition ####

com <- zoops_com %>% left_join(treat, by = c('pond' = 'pond.ID')) %>%
  mutate('zd' = rowSums(.[3:11])) %>%
  filter(zd > 0)
trt.sub <- select(com, pond, MC.ID:pH.var)
com <- select(com, Bosmina.longispina:Ceriodaphnia.dubia)
row.names(com) <- trt.sub$pond

dm <- vegdist(com,method = 'jaccard')
adonis(dm~trt.sub$pH.local*trt.sub$dispersal*trt.sub$pH.var+trt.sub$upstream)
dm <- vegdist(com,method = 'bray')
adonis(dm~trt.sub$pH.local*trt.sub$dispersal*trt.sub$pH.var)
dm <- vegdist(log1p(com),method = 'bray')
adonis(dm~trt.sub$pH.local*trt.sub$dispersal*trt.sub$pH.var)

ordi <- metaMDS(com, distance = 'bray', k = 2, autotransform = FALSE, trymax = 500)
g<-ordi$points[,1:2]

pdf('~/Desktop/nmds_plot.pdf',width = 5,height = 5,pointsize = 12)
plot(g[,2] ~ g[,1], type = "n",yaxt='n',xaxt='n',ann=F)
title(xlab='NMDS dimension 1',cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='NMDS dimension 2',cex.lab=1,line = 2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# pch.vec <- as.numeric(as.factor(trt.sub$dispersal))-1
# pch.vec[trt.sub$pH.var == 'heterogeneous'] <- pch.vec[trt.sub$pH.var == 'heterogeneous'] + 15
# points(g[,2] ~ g[,1],pch=pch.vec,col=cols[as.numeric(as.factor(trt.sub$pH.local))])
points(g[,2] ~ g[,1],pch=16,col=cols[as.numeric(as.factor(trt.sub$pH.local))])
labels <- ordi$species[,]
names <- rownames(labels) %>% str_replace('\\.',' ') %>%  make.italic
text(labels, names, cex = 0.7, col = 1)
legend('topright',bty='n',legend=bquote('Stress ='~.(round(ordi$stress,2))))
dev.off()

#####

pdf('~/Desktop/Figure.pdf',width=9,height=9,pointsize = 8)

par(mfrow=c(4,3),cex=1)
col2 <- c('#1b9e77','#d95f02','#7570b3')

gam.zoo.disp.tot <- bam(zd_l ~ pH*disp + s(stime, by = pH, k=7) + s(stime, by = disp, k=7) + s(stime, by = interaction(pH,disp), k=7) + s(stime, pond, bs='re',m=1,k=5), data=zoo, method='fREML')
gam.zoo.disp.clad <- bam(log_clad ~ pH*disp + s(stime, by = pH, k=7) + s(stime, by = disp, k=7) + s(stime, by = interaction(pH,disp), k=7) + s(stime, pond, bs='re',m=1,k=5), data=zoo, method='fREML')
gam.zoo.disp.cop <- bam(log_cop ~ pH*disp + s(stime, by = pH, k=7) + s(stime, by = disp, k=7) + s(stime, by = interaction(pH,disp), k=7) + s(stime, pond, bs='re',m=1,k=5), data=zoo, method='fREML')

#pH 8.5
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='8.5',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(crustacean density)'),main='pH 8.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='8.5',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='8.5',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)
legend('bottomleft',bty='n',legend = c('no dispersal','stepping stone','global dispersal'),pch=16,col=col2,y.intersp=1)

plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='8.5',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(cladoceran density)'),main='pH 8.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='8.5',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='8.5',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='8.5',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(copepod density)'),main='pH 8.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='8.5',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='8.5',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)


#pH 7
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='7',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(crustacean density)'),main='pH 7',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='7',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='7',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='7',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(cladoceran density)'),main='pH 7',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='7',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='7',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='7',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(copepod density)'),main='pH 7',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='7',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='7',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

#pH 5.5
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='5.5',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(crustacean density)'),main='pH 5.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='5.5',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='5.5',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='5.5',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(cladoceran density)'),main='pH 5.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='5.5',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='5.5',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='5.5',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(copepod density)'),main='pH 5.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='5.5',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='5.5',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

#pH 4
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='4',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(crustacean density)'),main='pH 4',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='4',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH='4',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='4',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(cladoceran density)'),main='pH 4',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='4',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH='4',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='4',disp='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(copepod density)'),main='pH 4',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='4',disp='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH='4',disp='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

dev.off()

#### ####

pdf('~/Desktop/boxplot.pdf',width=10,height=4,pointsize = 7.5)
par(mfrow=c(1,2))
zoo.last <- filter(zoo, week == 6)
boxplot(log_cop~disp*pH,zoo.last,main='copepod')
boxplot(log_clad~disp*pH,zoo.last,main='cladoceran')
dev.off()

