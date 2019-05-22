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

#### Load and format plankton data ####

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
  filter(pond %!in% to.rm) %>% group_by(date, pond) %>% summarise_if(is.numeric, mean)
phyto$date <- as.Date(phyto$date, format = '%d/%m/%Y')
phyto <- arrange(phyto, by = date)

phyto <- phyto %>% filter(nchar(pond) == 2)
phyto$week <- rep(c(0:10,12,16), each=96)

phyto$pH <- treat$pH.local[match(phyto$pond,treat$pond.ID)]
phyto$disp <- treat$dispersal[match(phyto$pond,treat$pond.ID)]
phyto$str <- treat$pH.var[match(phyto$pond,treat$pond.ID)]
phyto$pH <- factor(phyto$pH)
phyto$pond <- as.factor(phyto$pond)

#fixing the 1 zero value, which does not fit on a log plot - error on instrument is 0.1 ug, I add 0.01 ug
phyto$total[which.min(phyto$total)] <- 0.01
phyto$logchla <- log(phyto$total)

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
zoops_tot$date <- format(zoops_tot$date, '%j')
zoops_tot$date <- as.numeric(zoops_tot$date) - 143

zoops_com <- read.csv('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2017/data/zoops_community.csv', stringsAsFactors = F) 

zoops2 <- zoops_com %>% filter(pond %!in% to.rm) %>%
  filter(pond != 'NA') %>%
  mutate('Clad.perL' = (Eubosmina.longispina + Diaphanosoma + Sida.crystallina + Chydorus.sphaericus + Daphnia.pulex + Simocephalus + Daphnia.ambigua + Ceriodaphnia)/2, 'Cop.perL' = Cyclops.scutifer/2) %>%
  select(pond,date,Clad.perL,Cop.perL,density) %>%
  rename('zd' = density)
zoops2$date <- as.Date(zoops2$date, format = '%d.%m.%y')
zoops2$date <- format(zoops2$date, '%j')
zoops2$date <- as.numeric(zoops2$date) - 143

zoo <- bind_rows(zoops_tot,zoops2) %>% arrange(date)
rm(zoops2, zoops_tot)  

weeks <- data.frame('day' = c(16,23,30,37,44,51,58,65,72,79,86,100,126),
                    'week' = c(0:10,12,16))

zoo$week <- weeks$week[match(zoo$date,weeks$day)]
rm(weeks)
zoo$zd_l <- log1p(zoo$zd)
zoo$log_clad <- log1p(zoo$Clad.perL)
zoo$log_cop <- log1p(zoo$Cop.perL)

zoo$pH <- treat$pH.local[match(zoo$pond,treat$pond.ID)]
zoo$disp <- treat$dispersal[match(zoo$pond,treat$pond.ID)]
zoo$str <- treat$pH.var[match(zoo$pond,treat$pond.ID)]

zoo$pH <- factor(zoo$pH)
zoo$pond <- as.factor(zoo$pond)

#### Optional filters ####

#getting rid of time zero and Phase II data points

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

