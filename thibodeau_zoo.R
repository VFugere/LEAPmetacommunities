library(tidyverse)
library(readxl)
library(magrittr)

zoops <- read_excel('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/Thibodeau data/zoops.xls')
zoops %<>% group_by(Taxon, Treatment, Date) %>%
  summarize(dens = mean(Abundance.nb.L)) %>%
  mutate(logdens = log1p(dens))
zoops$d <- as.Date(zoops$Date, format='%yyyy-%mm-%dd')

#zoops <- filter(zoops, Treatment != 'Pulse')

ggplot(data = zoops, aes(x = d, y = dens)) + geom_line(aes(color = Taxon, lty = Treatment))
ggplot(data = zoops, aes(x = d, y = logdens)) + geom_line(aes(color = Taxon, lty = Treatment))

#### ordination ####

zoops_com <- read_excel('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/Thibodeau data/zoops.xls') %>%
  filter(day == 35) %>% select(-Date) %>% filter(Taxon != 'Total')

zoops_com$site <- paste0('bag',1:24)
colnames(zoops_com) <- tolower(colnames(zoops_com))

treat <- select(zoops_com, site, treatment)
com <- zoops_com %>% select(taxon,dens) %>% spread(key=taxon, value=dens, fill=0, drop=F)


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




#### relationship between chl. a and biovolume ####

algae <- read_excel('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/Thibodeau data/algae.xls')
algae %>% select(Ank.con:Epi.sp) %>% rowSums(na.rm = T) -> algae$totalbv
algae$totalbv <- algae$totalbv*10^-9

chla <- read_excel('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/Thibodeau data/chlorophyll_clean.xlsx') %>%
  rename(chla = 'Chl.a.(ug/L)')

phyto <- algae %>% select(Mesocosm,date,totalbv) %>% inner_join(chla, by = c('Mesocosm','date'))
plot(totalbv ~ chla, phyto, pch = 16, xlab='chlorophyll (ug/L)', ylab='biovolume (mm3/L)')
