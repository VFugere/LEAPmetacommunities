#### Vincent Fugere, 2018 - 2019
#### LEAP 2017 experiment
#### Code to examine effect of dispersal and pH heterogenety on zooplankton metacommunities

rm(list=ls())

library(tidyverse)
library(scales)
source('/Users/vincentfugere/Google Drive/Recherche/PhD/R/functions/utils.R')
library(scales)
library(mgcv)
library(itsadug)
library(magrittr)
library(vegan)

cols<-c('firebrick2','gold2','forestgreen','darkblue')
treat <- read.csv('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/LEAP2017treatments.csv')
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
  mutate('Clad.perL' = (Eubosmina.longispina + Diaphanosoma + Sida.crystallina + Chydorus.sphaericus + Daphnia.pulex + Simocephalus + Daphnia.ambigua + Ceriodaphnia)/2, 'Cop.perL' = Cyclops.scutifer/2) %>%
  select(pond,date,Clad.perL,Cop.perL,density) %>%
  rename('zd' = density)
zoops2$date <- as.Date(zoops2$date, format = '%d.%m.%y')
zoo <- bind_rows(zoops_tot,zoops2)
rm(zoops2, zoops_tot)

zoo$date <- as.numeric(format(zoo$date, '%j'))
zoo %<>% filter(date %in% date.range) %>% arrange(date,pond)
zoo$week <- rep(0:7,each=96)
zoo %<>% select(pond, week, Clad.perL:zd)

zoops_com <- filter(zoops_com, date == '27.07.17', pond %!in% to.rm) 

## binding and adding treatments

phyto$pH <- treat$pH.local[match(phyto$pond,treat$pond.ID)]
phyto$disp <- treat$dispersal[match(phyto$pond,treat$pond.ID)]
phyto$str <- treat$pH.var[match(phyto$pond,treat$pond.ID)]
phyto$pH <- factor(phyto$pH)
phyto$pond <- as.factor(phyto$pond)


#### Optional filters ####

#getting rid of pre-treatment samples and Phase II data points

#getting rid of two ponds which were accidently acidified to a lower pH than intended

#getting rid of homogeneous metacommunities for now for simplicty (fix later)

#### Question 1: does dispersal & and local pH 
#zoo

zoo <- zoo %>% filter(week < 7, str == 'heterogeneous')
zoo$stime <- zoo$week

# m8 <- bam(zd_l ~ s(stime, k=7) + s(stime, pond, bs='re',m=1,k=5), data=zoo, method='ML')
# m9 <- bam(zd_l ~ pH + s(stime, by = pH, k=7) + s(stime, pond, bs='re',m=1,k=5), data=zoo, method='ML')
# m10 <- bam(zd_l ~ s(stime, by = disp, k=7) + s(stime, pond, bs='re',m=1,k=5), data=zoo, method='ML')
# m11 <- bam(zd_l ~ s(stime, by = str, k=7) + s(stime, pond, bs='re',m=1,k=5), data=zoo, method='ML')
# m12 <- bam(zd_l ~ pH*disp + s(stime, by = pH, k=7) + s(stime, by = disp, k=7) + s(stime, by = interaction(pH,disp), k=7) + s(stime, pond, bs='re',m=1,k=5), data=zoo, method='ML')
# m13 <- bam(zd_l ~ s(stime, by = pH, k=7) + s(stime, by = str, k=7) + s(stime, by = interaction(pH,str), k=7) + s(stime, pond, bs='re',m=1,k=5), data=zoo, method='ML')
# m14 <- bam(zd_l ~ s(stime, by = pH, k=7) + s(stime, by = disp, k=7) + s(stime, by = str, k=7) + s(stime, by = interaction(pH,disp,str), k=7) + s(stime, pond, bs='re',m=1,k=5), data=zoo, method='ML')
#
#AIC(m8,m9,m10,m11,m12,m13,m14)
#AIC(m8,m9,m10,m11,m12,m13,m14)$AIC - min(AIC(m8,m9,m10,m11,m12,m13,m14)$AIC)


#### ####

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

#### -- ####

zoops_com <- zoops_com %>% filter(date == '27.07.17', pond %!in% to.rm) %>% 
  filter(abundance != 0)
com <- zoops_com %>% select(Eubosmina.longispina:Ceriodaphnia) %>% as.matrix
row.names(com) <- zoops_com$pond
com.trt <- zoops_com %>% select(pond) %>% as.data.frame
com.trt <- left_join(com.trt, treat, by = c('pond'='pond.ID')) %>% select(pond, MC.ID:pH.var )

dm <- vegdist(log1p(com),method = 'bray')
adonis(dm~com.trt$pH.local+com.trt$dispersal:com.trt$pH.var+com.trt$dispersal:com.trt$pH.var:com.trt$pH.local)

