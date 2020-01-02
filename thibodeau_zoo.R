library(tidyverse)
library(readxl)
library(magrittr)
library(vegan)

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

zoops_com$site <- rep(paste0('bag',1:9),7)
colnames(zoops_com) <- tolower(colnames(zoops_com))

treat <- select(zoops_com[1:9,], site, treatment)
com <- zoops_com %>% select(-day,-biomass.ug.l, -treatment) %>%
  spread(key=taxon, value=abundance.nb.l)

com <- select(com, Bosmina:Nauplii) %>% as.matrix
row.names(com) <- treat$site

dm <- vegdist(com,method = 'jaccard')
dm <- vegdist(com,method = 'bray')
dm <- vegdist(log1p(com),method = 'bray')
adonis(dm~treat$treatment)

ordi <- metaMDS(com, distance = 'jaccard', binary=T, k = 2, autotransform = T, trymax = 500)
g<-ordi$points[,1:2]

treat <- cbind(treat,g)
boxplot(MDS1~treatment,treat)
boxplot(MDS2~treatment,treat)

cols<-c('firebrick2','gold2','forestgreen','darkblue')
cols <- cols[c(4,3,2)]
#pdf('~/Desktop/nmds_plot.pdf',width = 5,height = 5,pointsize = 12)
plot(g[,2] ~ g[,1], type = "n",yaxt='n',xaxt='n',ann=F)
title(xlab='NMDS dimension 1',cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='NMDS dimension 2',cex.lab=1,line = 2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
points(g[,2] ~ g[,1],pch=16,col=cols[as.numeric(as.factor(treat$treatment))])
labels <- ordi$species[,]
names <- rownames(labels)
text(labels, names, cex = 0.7, col = 1)
legend('topright',bty='n',legend=bquote('Stress ='~.(round(ordi$stress,2))))
#dev.off()


#### relationship between chl. a and biovolume ####

algae <- read_excel('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/Thibodeau data/algae.xls')
algae %>% select(Ank.con:Epi.sp) %>% rowSums(na.rm = T) -> algae$totalbv
algae$totalbv <- algae$totalbv*10^-9

chla <- read_excel('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/Thibodeau data/chlorophyll_clean.xlsx') %>%
  rename(chla = 'Chl.a.(ug/L)')

phyto <- algae %>% select(Mesocosm,date,totalbv) %>% inner_join(chla, by = c('Mesocosm','date'))
plot(totalbv ~ chla, phyto, pch = 16, xlab='chlorophyll (ug/L)', ylab='biovolume (mm3/L)')
