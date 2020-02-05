#### Vincent Fugere, 2018 - 2020
#### LEAP 2017 experiment
#### Code to examine effect of dispersal and pH heterogeneity on zooplankton metacommunities

rm(list=ls())

library(tidyverse)
library(readxl)
library(writexl)
library(magrittr)
library(scales)
library(chron)

library(brms)
library(vegan)
library(mgcv)
library(itsadug)
library(ape)

source('/Users/vincentfugere/Google Drive/Recherche/PhD/R/functions/utils.R')

## get data

source('thibodeau_zoo.R')
source('LEAP_data_format.R')

#### 1. impact of acidification on diversity and CC, without any effect of dispersal.

pdf('~/Desktop/effect_of_pH.pdf',width=8.5 ,height = 11, pointsize = 12, onefile = T)

#thibodeau

par(mfrow=c(3,2),cex=0.7,mar=c(4,4,0.5,0.5))
boxplot(density~treatment,thibodeau,ylab='total density')
boxplot(richness~treatment,thibodeau,ylab='richness')
boxplot(hill~treatment,thibodeau,ylab=expression(e^Shannon))
boxplot(Nauplii~treatment,thibodeau,ylab='% nauplii')
boxplot(bosm_cerio~treatment,thibodeau,ylab='% bosm+cerio')
boxplot(chyd_daph_diaphan~treatment,thibodeau,ylab='% chyr+daph+diaph')
par(mfrow=c(1,1))

#LEAP, end of Phase 1 only, no dispersal only

to.keep <- filter(treat, dispersal == 'N') %>% select(pond.ID,pH.local) %>% as.data.frame
leap <- zoops_com %>% filter(pond %in% to.keep$pond.ID) %>% left_join(to.keep, c('pond' = "pond.ID"))
empty.columns <- names(leap[2:15])[colSums(leap[,2:15]) == 0]
leap <- leap %>% select(-empty.columns)
leap.com <- leap %>% select(`Bosmina longirostris`:copepodids)

leap$density <- apply(leap.com, 1, 'sum')/2
leap$richness <- specnumber(leap.com)
leap$hill <- exp(diversity(leap.com, index='shannon'))

par(mfrow=c(5,3),cex=0.7,mar=c(4,4,0.5,0.5))
boxplot(density~pH.local,leap,ylab='total density',xlab='pH')
boxplot(richness~pH.local,leap,ylab='richness',xlab='pH')
boxplot(hill~pH.local,leap,ylab=expression(e^Shannon),xlab='pH')
for(i in 2:11){
  sub <- leap[,c(i,12)]
  species <- names(sub)[1]
  boxplot(log1p(sub[,1])~sub[,2],ylab=species,xlab='pH')
}
par(mfrow=c(1,1))

dev.off()

#### Exploring pH, disp, and effects on local communities ####

#data <- filter(data, week < 7)

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
  plot(y~week,tempdata,type='n',yaxt='n',xaxt='n',xlim=c(0,7),ylim=range(tempdata$y),ann=F,bty='l')
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
  plot(y~week,tempdata,type='n',yaxt='n',xaxt='n',xlim=c(0,7),ylim=range(tempdata$y),ann=F,bty='l')
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

library(lme4)
m1 <- lmer(zd_log ~ pH.trt_f + disp_f + str_f + pH.trt_f:disp_f + disp_f:str_f + pH.trt_f:disp_f:str_f + (1|MC_f/pond_f), data)
summary(m1)
anova(m1)
plot(m1)

dd <- hete
dd$var <- dd$total_log
m1 <- update(m1, newdata = dd)
#m1 <- brm(var ~ pH.trt_f * disp_f + (1|MC_f/pond_f), dd, cores=4)
summary(m1)
plot(m1)
pp_check(m1)
marginal_effects(m1)

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
  filter(zd > 0) %>% 
  filter(pH.var == 'heterogeneous')
trt.sub <- select(com, pond, MC.ID:pH.var)
com <- select(com, `Bosmina longirostris`:copepodids)
row.names(com) <- trt.sub$pond

dm <- vegdist(com,method = 'jaccard')
adonis(dm~trt.sub$pH.local*trt.sub$dispersal*trt.sub$pH.var+trt.sub$upstream)
dm <- vegdist(com,method = 'bray')
adonis(dm~trt.sub$pH.local*trt.sub$dispersal*trt.sub$pH.var)
dm <- vegdist(log1p(com),method = 'bray')
adonis(dm~trt.sub$pH.local*trt.sub$dispersal*trt.sub$pH.var)

ordi <- metaMDS(com, distance = 'bray', k = 2, autotransform = FALSE, trymax = 500)
g<-ordi$points[,1:2]

trt.sub <- cbind(trt.sub,g)
boxplot(MDS1~dispersal*pH.local,trt.sub)
boxplot(MDS2~dispersal*pH.local,trt.sub)
with(trt.sub, table(pH.local,dispersal))

par(mfrow=c(1,2))
zoo.last <- filter(zoo, week == 6)
boxplot(log_cop~disp*pH,zoo.last,main='copepod')
boxplot(log_clad~disp*pH,zoo.last,main='cladoceran')

pdf('~/Desktop/nmds_plot.pdf',width = 5,height = 5,pointsize = 12)
plot(g[,2] ~ g[,1], type = "n",yaxt='n',xaxt='n',ann=F)
title(xlab='NMDS dimension 1',cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='NMDS dimension 2',cex.lab=1,line = 2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
pch.vec <- as.numeric(as.factor(trt.sub$dispersal))-1
pch.vec[trt.sub$pH.var == 'heterogeneous'] <- pch.vec[trt.sub$pH.var == 'heterogeneous'] + 15
points(g[,2] ~ g[,1],pch=pch.vec,col=cols[as.numeric(as.factor(trt.sub$pH.local))])
points(g[,2] ~ g[,1],pch=c(15,16,17)[as.numeric(as.factor(trt.sub$dispersal))],col=cols[as.numeric(as.factor(trt.sub$pH.local))])
labels <- ordi$species[,]
names <- rownames(labels) %>% str_replace('\\.',' ') %>%  make.italic
text(labels, names, cex = 0.7, col = 1)
legend('topright',bty='n',legend=bquote('Stress ='~.(round(ordi$stress,2))))
dev.off()

##### gamms ####

pdf('~/Desktop/hete.pdf',width=9,height=9,pointsize = 8)

par(mfrow=c(4,3),cex=1)
col2 <- c('#1b9e77','#d95f02','#7570b3')

zoo <- hete
zoo$stime <- scale(zoo$week)

gam.zoo.disp.tot <- bam(zd_log ~ pH.trt_f*disp_f + s(stime, by = pH.trt_f, k=7) + s(stime, by = disp_f, k=7) + s(stime, by = interaction(pH.trt_f,disp_f), k=7) + s(stime, pond_f, bs='re',m=1,k=5), data=zoo, method='fREML')
gam.zoo.disp.clad <- bam(Clad.perL_log ~ pH.trt_f*disp_f + s(stime, by = pH.trt_f, k=7) + s(stime, by = disp_f, k=7) + s(stime, by = interaction(pH.trt_f,disp_f), k=7) + s(stime, pond_f, bs='re',m=1,k=5), data=zoo, method='fREML')
gam.zoo.disp.cop <- bam(Cop.perL_log ~ pH.trt_f*disp_f + s(stime, by = pH.trt_f, k=7) + s(stime, by = disp_f, k=7) + s(stime, by = interaction(pH.trt_f,disp_f), k=7) + s(stime, pond_f, bs='re',m=1,k=5), data=zoo, method='fREML')

#pH 8.5
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='8.5',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(crustacean density)'),main='pH 8.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='8.5',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='8.5',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)
legend('bottomleft',bty='n',legend = c('no dispersal','stepping stone','global dispersal'),pch=16,col=col2,y.intersp=1)

plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='8.5',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(cladoceran density)'),main='pH 8.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='8.5',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='8.5',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='8.5',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(copepod density)'),main='pH 8.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='8.5',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='8.5',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)


#pH 7
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='7',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(crustacean density)'),main='pH 7',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='7',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='7',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='7',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(cladoceran density)'),main='pH 7',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='7',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='7',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='7',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(copepod density)'),main='pH 7',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='7',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='7',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

#pH 5.5
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='5.5',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(crustacean density)'),main='pH 5.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='5.5',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='5.5',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='5.5',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(cladoceran density)'),main='pH 5.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='5.5',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='5.5',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='5.5',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(copepod density)'),main='pH 5.5',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='5.5',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='5.5',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

#pH 4
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='4',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(crustacean density)'),main='pH 4',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='4',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.tot, view="stime", cond=list(pH.trt_f='4',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='4',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(cladoceran density)'),main='pH 4',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='4',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.clad, view="stime", cond=list(pH.trt_f='4',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='4',disp_f='N'), col=col2[1], rm.ranef=T, hide.label = T,xlab='weeks',ylab=expression(log[e]~'1+(copepod density)'),main='pH 4',rug=F,bty='l',legend_plot_all = F, h0=NA,ylim=c(-1,5))
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='4',disp_f='L'), col=col2[2], rm.ranef=T, add=T, rug=F)
plot_smooth(gam.zoo.disp.cop, view="stime", cond=list(pH.trt_f='4',disp_f='G'), col=col2[3], rm.ranef=T, add=T, rug=F)

dev.off()
