library(tidyverse)
library(readxl)
library(magrittr)

zoops <- read_excel('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/Thibodeau data/zoops.xls')
zoops %<>% group_by(Taxon, Treatment, Date) %>%
  summarize(dens = mean(Abundance.nb.L)) %>%
  mutate(logdens = log1p(dens))
zoops$d <- as.Date(zoops$Date, format='%yyyy-%mm-%dd')

zoops <- filter(zoops, Treatment != 'Pulse')

ggplot(data = zoops, aes(x = d, y = dens)) + geom_line(aes(color = Taxon, lty = Treatment)) + theme(plot.subtitle = element_text(vjust = 1), 
    axis.line = element_line(size = 0.5,  linetype = "solid"),
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    panel.background = element_rect(fill = NA))

ggplot(data = zoops, aes(x = d, y = logdens)) + geom_line(aes(color = Taxon, lty = Treatment))

algae <- read_excel('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/Thibodeau data/algae.xls')
algae %>% select(Ank.con:Epi.sp) %>% rowSums(na.rm = T) -> algae$totalbv
algae$totalbv <- algae$totalbv*10^-9

chla <- read_excel('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/Thibodeau data/chlorophyll_clean.xlsx') %>%
  rename(chla = 'Chl.a.(ug/L)')

#relationship between chl. a and biovolume

phyto <- algae %>% select(Mesocosm,date,totalbv) %>% inner_join(chla, by = c('Mesocosm','date'))

plot(totalbv ~ chla, phyto, pch = 16, xlab='chlorophyll (ug/L)', ylab='biovolume (mm3/L)')
