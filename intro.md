### Introduction  ###

Phase 1 of the LEAP 2017 experiment is an experiment in itself to answer interesting ecological questions. The main question is: do dispersal, habitat heterogeneity

 whether dispersal in metacommunities of 4 ponds can buffer local communities against local environmental degradation (acidification). Previous tests of the spatial insurance hypothesis with plankton communities in mesocosms have shown that dispersal from the regional species pool can mitigate impacts of stressors, whereby communities with immigration were less impacted than isolated communities without immigration. Examples of such experiments include [Thompson et al. (2012) _J. Anim. Ecol._](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2656.2011.01908.x) and [Limberger et al. (2019) _Ecol. Lett._](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13365) In these experiments, the dispersal treatment included an inoculum that pooled organisms from many local sites, i.e. the experimenters sampled many lakes in the study region, mixed all samples in a container, and then used this mixture to inoculate "ponds with dispersal" intermittently. Control ponds received no such inoculum. Basically, filling mesocosms with water originating from many lakes produces communities that are more resistant to disturbances than mesososms filled with water from a single lake. In a sense, this could be considered an effect of both dispersal (removing all barriers to dispersal within a landscape) and diversity (increasing the size of the species pool and thus the mean local diversity of communities). I understand that this is part of the spatial insurance hypothesis, but perhaps landscapes with a higher mean local diversity are more resilient irrespective of the presence of dispersal among the local communities that they contain. Andy, am I missing something?

In our experiment, we filled all ponds/metacommunities with the same species pool. We assume that the initial diversity of all local communities is equivalent. Then, we only manipulate dispersal within metacommunities, such that the only difference between no vs. local vs. global dispersal is the movement of individuals. The fact that we manipulate the spatial structure of dispersal (local vs global), not only its presence, is also novel. Our objective is to determine whether a 1% weekly dispersal event can mitigate some of the negative effects of acidification on zooplankton communities. In this report, I analyze our recently-obtained zooplankton composition data for LEAP 2017, and combine this with previously-published data from a lake mesocosm experiment done by a student from Bea's lab. To me the two experiments  are linked: in Bea's experiment, we see the effect of acidification on Lac Hertel zoo communities under conditions that are closest to the natural lake conditions. We notice negative effects of acidification. We then ask: can dispersal prevent some of these effects? To answer this question, we need a metacommunity-scale experiment with many more mesocosms (LEAP), which inevitably trades-off some realism for replication.

### Datasets ###

- Zooplankton community composition data from [Thibodeau et al. (2016) _Proc B._](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2015.1215) This is an acification experiment done on the floating dock on Lac Hertel. The phytoplankton response has been analyzed thoroughly but the zooplankton composition data has not been analyzed (or at least no analysis was published). We have data from 9 in-lake mesocosms: 3 served as controls, 3 were temporarily acidified to a pH of 5 (pulse perturbation) with a single acidification event, and 3 were maintained at a pH of 5 over several weeks with regular titrations (press perturbation). I will refer to this dataset as the 'Thibodeau' dataset throughout.

- LEAP 2017, Phase 1, before 'lethal acidification' to a pH of 3 in Phase 2. We have ~ 8 weeks of data during which weekly dispersal and pH treatments were applied. 96 ponds were arrayed into 4-pond metacommunities varying in their mean pH, stress heterogeneity (gradient of 4 pH vs. all 4 ponds maintained at the same pH), presence vs. absence of dispersal, and dispersal mode (local vs. global). Local dispersal in this context would represent a landscape in which downstream polluted lakes receive immigrants from upstream pristine lakes. Global dispersal provides equal dispersal probability among all lakes within a region. Homogeneous metacommunities in which all ponds have the same pH serve as controls to test for the effect of dispersal _per se_, in a homogeneous landscape. You can access a drawing of the design [**here**](https://www.dropbox.com/s/eknyap3z2fxlat7/LEAP%202017%20-%20design.pdf?dl=0). The response variables are:

  1. Zoo abundance: total abundance of copepods and cladocerans. Weekly measurements.
  2. Zoo composition: species-level abundance at the end of the experiment. 1 measurement per pond.
  3. Fluoroprobe: biomass of browns & greens, and total chlorophyll. Cyanos & cryptos very rare so excluded. Weekly measurements (8 measurements per pond).
  4. Some physico-chemical parameters and estimates of ecosystem metabolism. 3-8 measurements per pond depending on variable.

Dispersal is predicted to have different effects in homogeneous vs. heterogeneous metacommunities. Morever, these two types of metacommunities differ in the number of replicates we have for each dispersal treatment. Because we only included a single homogeneous metacommunity of 4 ponds at each pH level, we only have one replicate metacommunity per pH*dispersal treatment combination for homogeneous metacommunities. In contrast, for heterogeneous metacommunities, we have 4 replicate metacommunities for each dispersal treatment, i.e. 16 ponds per dispersal treatment. For these reasons, I often analyze homogeneous and heteregeneous metacommunities separately. Our predictions are mostly for heterogeneous metacommunities, such that homogeneous metacommunities can be seen as a form of control treatment. This distinction is perhaps confusing and unnecessary: perhaps we could simply report the results from the 48 ponds from heterogeneous MCs, although it would be a little sad not to use the full dataset. Let me know what you think about this.

<!-- ### Predictions & Questions ### -->

<!-- Does acidification impact zooplankton communities? -->

<!-- 1. Acifidication to a low pH (< 6) will affect zooplankton composition and have negative effects on zooplankton density and diversity (alpha). This will be visible in both the Thibodeau dataset (especially for the more extreme 'press' treatment) and LEAP dataset. At LEAP, we push the stress gradient farther (pH of 4), which might lead to complete extinction (we know this is the case but we could also have guessed it from the acidification literature). -->

<!-- 2. Beta diversity in metacommunities should also diminish with pH stress, as harsher environments are expected to filter out a significant share of the regional pool, thus reducing the importance of ecological drift in community assembly (e.g. [Chase 2007](https://www.pnas.org/content/104/44/17430.short)). This effect will be visible in homogeneous metacommunities containing 4 ponds of a single pH. If that pH is lower, then so should beta diversity within the MC. All heteregenous MCs have the same pH gradient (4, 5.5, 7, 8.5) so beta diversity would not be affected by pH in these MCs. -->

<!-- Can dispersal prevent some of these negative effects? (LEAP only) -->

<!-- 3. In hetereogeneous metacommunities, local diversity and abundance of ponds with a pH of 4 or 5.5 should be higher when they are connected to pH 7-8.5 ponds with dispersal than when they are isolated. We will test this by measuring final density and diversity at the end of Phase 1, and stability over Phase 1: all should be higher in MCs with dispersal. -->

<!-- 4. In homogeneous metacommunities, at any given stress level (pH), dispersal should lower beta diversity. I am not sure what would be the prediction for local diversity. It could go either way: a positive effect of dispersal on alpha diversity but only in stressful environments, i.e. if resistant species are rare and have a low probability of colonizing mesocosms at the onset of the experiment. Or, a negative effect of dispersal on local diversity, assuming that dispersal will spread competitively dominant species. Its hard to guess given that our "1% dispersal rate" is a weekly transfer, not a continuous flow of dispersers as in mathematical models based on which predictions can be made. -->

<!-- 5. Global vs. local: not sure about this one. In low pH ponds, local dispersal should have stronger effects on local diversity and density than global dispersal since a greater volume of water is transferred from a less degraded environment. For example, a pH 5.5 pond would receive ~ 10L from a pH 7 pond under local dispersal, but only 2.5L under global dispersal (+2.5 L from a pH8 pond, 2.5L from itself, but also 2.5L from a pH4 pond which might lack plankton altogether).  -->

### Effects of pH on zooplankton ###

Lets first look at this in the Thibodeau dataset, quantifying for each bag: 1) total zooplankton density, 2) species richness, 3) alpha diversity (Hill number, or exponent of Shannon index), and 4) community composition. I could not fit an ordination to this data because there are only 9 mesocosms with 7 species, and treatment effects are very strong (such that many points sit on top of each other in multivariate space). By simply looking at the community matrix, I noticed that: 1) Nauplii disappear completely in the press treatment but increase in abundance to > 95% relative abundance in the pulse treatment; 2) the cladocerans *Bosmina* and *Ceriodaphnia* disappear completely in both acid treatments; 3) in the press treatment, each bag is dominated (>80% relative abundance) by a different cladoceran, one of: _Chydorus_, _Daphnia_, or _Diaphanosoma_. Density and diversity are also affected by the treatments. These results are shown in Figure 1.

```{r}
source(here::here('./thibodeau_zoo.R'))
```

```{r, include=TRUE, fig.width=5,fig.height=6,fig.cap="Impacts of acidification on zooplankton communities in in-lake mesocosms"}
par(mfrow=c(3,2),cex=0.7,mar=c(4,4,0.5,0.5))
boxplot(density~treatment,thibodeau,ylab='total density')
boxplot(richness~treatment,thibodeau,ylab='richness')
boxplot(hill~treatment,thibodeau,ylab=expression(e^Shannon))
boxplot(Nauplii~treatment,thibodeau,ylab='% nauplii')
boxplot(bosm_cerio~treatment,thibodeau,ylab='% bosm+cerio')
boxplot(chyd_daph_diaphan~treatment,thibodeau,ylab='% chyr+daph+diaph')
par(mfrow=c(1,1))
```

```{r}
source(here::here('./LEAP_data_format.R'))
```

Lets now do a similar analysis with the LEAP data, only including isolated ponds from metacommunities without dispersal. This shows the effect of pH in the absence of dispersal. Again, we look at density, diversity, and composition.
