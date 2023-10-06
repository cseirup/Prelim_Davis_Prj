

# Allometric equations ----------------------------------------------------

#Complete tree biomass: Young's complete tree, Aaron Teets code; DBH=diameter at breast height, includes bark, to biomass in kilograms, comp = complete
compmass=function(Species, DBH){
  if(Species=='PIRU'){comp.b0=1.10651;comp.b1=2.298388}
  else if(Species=='ABBA'){comp.b0=0.8161535;comp.b1=2.413978}
  else if(Species=='PIGL'){comp.b0=1.10651;comp.b1=2.298388}#using PIRU values
  else if(Species=='PIST'){comp.b0=0.5724865;comp.b1=2.467798}
  else if(Species=='TSCA'){comp.b0=0.8645225;comp.b1=2.38591}
  else if(Species=='THOC'){comp.b0=1.32942;comp.b1=1.919051}
  else if(Species=='BEAL'){comp.b0=1.345053;comp.b1=2.335473}
  else if(Species=='BEPA'){comp.b0=0.74343;comp.b1=2.639837}
  else if(Species=='BECO'){comp.b0=0.74343;comp.b1=2.639837}#using BEPA values
  else if(Species=='ACRU'){comp.b0=1.187596;comp.b1=2.37025}
  else if(Species=='AMELANCHIER'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='ACPE'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='ACSP'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='ACER'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='SODE'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='PRPE'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='UNKCON'){comp.b0=1.10651;comp.b1=2.298388} #using PIRU
  else {comp.b0=0;comp.b1=0}
  comp=exp(comp.b0+comp.b1*(log(DBH/2.54)))*0.4536
  return(comp)}

#calculating carbon mass from carbon mass using % from Lamlom and Savidge 2003
compcarb=function(Species, biomass){
  if(Species=='PIRU'){per.carb = 0.5039} #using picea glauca 0.5039
  else if(Species=='PIGL'){per.carb = 0.5039}
  else if(Species=='ABBA'){per.carb = 0.5008}
  else if(Species=='PIST'){per.carb = 0.4974}
  else if(Species=='TSCA'){per.carb = 0.5033}
  else if(Species=='THOC'){per.carb = 0.5172}
  else if(Species=='BEAL'){per.carb = 0.4627}
  else if(Species=='BEPA'){per.carb = 0.4837}
  else if(Species=='BECO'){per.carb = 0.4837}#using BEPA
  else if(Species=='ACRU'){per.carb = 0.4864}
  else if(Species=='ACER'){per.carb = 0.4864} #using ACRU
  else if(Species=='AMELANCHIER'){per.carb = 0.4864} #using ACRU
  else if(Species=='ACPE'){per.carb = 0.4864} #using ACRU
  else if(Species=='ACSP'){per.carb = 0.4864} #using ACRU
  else if(Species=='SODE'){per.carb = 0.4864} #using ACRU
  else if(Species=='PRPE'){per.carb = 0.4953} #using Prunus serotina
  else if(Species=='UNKCON'){per.carb = 0.5039} #using Picea glauca
  else {per.carb = 0.0}
  carbon=per.carb*biomass
  return(carbon)}


#Function to sub sample the 10x200m plot as 5 10x10m plots (same as Davis). Just trees.
sample_fun_trees <- function(df, site){
  df2 <- filter(df, Site == site)
  plots <- data.frame(SubPlotID = unique(df2[,"SubPlotID"]))
  colnames(plots) <- "SubPlotID" # bug handling for purrr::map
  n <- 5

  samp <- data.frame(SubPlotID = sample(plots$SubPlotID, n, replace = FALSE)) %>%
    dplyr::arrange(SubPlotID)

  # Set up unique naming column, so plots selected more than once have a unique ID.
  samp$case <- as.factor(stringr::str_pad(rownames(samp), nchar(n), side ="left", pad = 0))
  
  #Joining list of selected subplots back to their data
  df_samp <- dplyr::left_join(samp, df2, by = "SubPlotID") %>%
    dplyr::arrange(case)

  #calculating total density, basal area, biomass for the 5 subplots selected
  df_calc <- df_samp %>%  summarise(Site = first(Site),
                                    sum_stems = sum(num_stem),
                                    sum_BA_cm2 = sum(BA_cm2),
                                    sum_carb = sum(carbon*num_stem),
                                    sum_comp = sum(comp*num_stem)) 
  
  #converting from plot to ha and kg to megagram. Total area for the 5 subplots = 500m
  df_calc2 <- df_calc %>%  mutate(density_ha = sum_stems * (10000/500), 
                                  BA_m2ha = sum_BA_cm2/500,
                                  carbonmass_Mgha = (sum_carb * (10000/500))/1000, 
                                  biomass_Mgha = (sum_comp * (10000/500))/1000) %>% 
                           mutate(across(where(is.numeric), round, 0)) %>% 
                           select(Site, density_ha, BA_m2ha, biomass_Mgha, carbonmass_Mgha) %>% data.frame()
}


#Function to calculate confidence intervals on bootstrap
bootCI <- function(df){
  boot_CI <- df %>% select(-Site)
  boot_CIs <- data.frame(t(apply(boot_CI, 2,
                                 quantile, probs = c(0.025, 0.5, 0.975), na.rm = T)))
  num_boot_fail <- num_reps - max(boot_CI$boot)
  
  boot_CIs$metric <- rownames(boot_CIs)
  boot_CIs$num_boots <- max(boot_CI$boot)
  boot_CIs$Site <- first(df$Site)
  colnames(boot_CIs) <- c("lower95", "median", "upper95", "metric", "num_boots", "Site")
  print(boot_CIs)
}


#Function for having unlabeled tickmarks on graph axis
every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}


#Species composition plotting functions

#For making the tree plot
Comptreeplot <- function(df, SiteCode, sp_group){ # df = data frame, Site = 2 letter site code, sp_group = value containing selected species for the site
  df2 <- df %>% filter(Site == SiteCode & Species %in% sp_group) #%>% 
               # mutate(Species = fct_reorder(Species, desc(BA_m2ha)))
  ggplot(aes(x = Species, y = BA_m2ha, fill = SampleEventNum), data = df2)+
    geom_col(position = position_dodge2(preserve = "single"))+
    #xlab('Species')+
    ylab(SiteCode)+ 
    theme(axis.text = element_text(size = 9), # change axis label size
          strip.text = element_text(size = 10), # change facet text size
          axis.title.x = element_blank(), # change axis title size
          axis.text.x = element_text(size = 10, face = "italic"),
          #axis.title.y = element_blank(),
          legend.position = "none")+
    scale_x_discrete(labels = sp_names, guide = guide_axis(n.dodge=3))+
    scale_fill_manual(name = "Year", labels = c("1" = '1959', "2" = '2020s'), 
                      values = c("1" = '#a1d99b', "2" = '#31a354'))+
    theme_FHM() 
}

#For making the sapling plot
Compsapplot <- function(df, SiteCode, sp_group){ # df = data frame, Site = 2 letter site code, sp_group = value containing selected species for the site
  df2 <- df %>% filter(Site == SiteCode & Species %in% sp_group) #%>% 
   # mutate(Species = fct_reorder(Species, desc(BA_m2ha)))
  ggplot(aes(x = Species, y = BA_m2ha, fill = SampleEventNum), data = df2)+
    geom_col(position = position_dodge2(preserve = "single"))+
    #xlab('Species')+
   # ylab(bquote('Basal area ('~m^2*'/ha)'))+ 
    theme(axis.text = element_text(size = 9), # change axis label size
          strip.text = element_text(size = 10), # change facet text size
          axis.title.x = element_blank(), # change axis title size
          axis.text.x = element_text(size = 10, face = "italic"),
          axis.title.y = element_blank(),
          legend.position = "none")+
    scale_x_discrete(labels = sp_names, guide = guide_axis(n.dodge=3))+
    scale_fill_manual(name = "Year", labels = c("1" = '1959', "2" = '2020s'), 
                      values = c("1" = '#a1d99b', "2" = '#31a354'))+
    theme_FHM() 
}

#For making the seedling plot
Compseedplot <- function(df, SiteCode, sp_group){ # df = data frame, SiteCode = 2 letter site code, sp_group = value containing selected species for the site
  df2 <- df %>% filter(Site == SiteCode & Species %in% sp_group)# %>% 
   # mutate(Species = fct_reorder(Species, desc(den_m2)))
  ggplot(aes(x = Species, y = den_m2, fill = SampleEventNum), data = df2)+
    geom_col(position = position_dodge2(preserve = "single"))+
    theme(axis.text = element_text(size = 9), # change axis label size
          strip.text = element_text(size = 10), # change facet text size
          axis.title.x = element_blank(), # change axis title size
          axis.text.x = element_text(size = 10, face = "italic"),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10))+
    scale_x_discrete(labels = sp_names, guide = guide_axis(n.dodge=3))+
    scale_fill_manual(name = "Year", labels = c("1" = '1959', "2" = '2020s'), 
                      values = c("1" = '#a1d99b', "2" = '#31a354'))+
    theme_FHM() 
}

#For ranking species composition by basal area
rankBA <- function(df, SiteCode){ # df = data frame, SiteCode = 2 letter site code,
  df %>% filter(Site == SiteCode & BA_m2ha>0 ) %>% 
  group_by(Species, SampleEventNum) %>% 
  summarize(BA_m2ha = sum(BA_m2ha)) %>% 
  arrange(desc(BA_m2ha))
}

#For ranking species composition by basal area
rankDEN <- function(df, SiteCode){ # df = data frame, SiteCode = 2 letter site code,
  df %>% filter(Site == SiteCode & den_m2>0 ) %>% 
    group_by(Species, SampleEventNum) %>% 
    summarize(den_m2 = sum(den_m2)) %>% 
    arrange(desc(den_m2))
}

