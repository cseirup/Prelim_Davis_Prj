#code for some simple graphs of expected results

library(tidyverse)
library(readxl)

df <- read_excel("../data/Exp_results_biomass.xlsx")
df_long <- df %>% pivot_longer(-Site, names_to = "Sample_Year", values_to = "biomass")

biom_plot <- ggplot(df_long, aes(x = Sample_Year, y = biomass, group = Site, colour = Site))+
              geom_point(size = 1.5) +
              geom_line(size = 1) +
              theme_bw()+
              ylab("Hypothetical Above Ground Biomass")+
              xlab("Sample Year")+
              expand_limits(y=0)+
              expand_limits(x=1)+
              scale_x_discrete(expand = expansion(mult = c(0.3, .3)))

biom_plot
ggsave("./figures/exp_biomass.jpg", biom_plot, width = 4, height = 4)


moss <- read_excel("../data/Exp_results_bryo.xlsx")
moss_long <- moss %>% pivot_longer(-Site, names_to = "Sample_Year", values_to = "moss")

moss_plot <- ggplot(moss_long, aes(x = Sample_Year, y = moss, group = Site, colour = Site))+
            geom_point(size = 1.5) +
             geom_line(size = 1) +
              theme_bw()+
              expand_limits(y=0)+
              ylab("Hypothetical % Bryophyte Cover")+
              xlab("Sample Year")+
              expand_limits(x=1)+
              scale_x_discrete(expand = expansion(mult = c(0.3, .3)))
  

moss_plot
ggsave("./figures/exp_moss.jpg", moss_plot, width = 4, height = 4)


netn <- read_excel("../data/Exp_results_netn.xlsx")

netn_plot <- ggplot(netn, aes(x = X, y = Y, colour = Plot_type))+
  geom_point(size = 1.5)+
  theme_bw()+
  expand_limits(y=c(-3,9))+
  expand_limits(x=c(-3,9))
 

netn_plot
ggsave("./figures/exp_netn.jpg", netn_plot)


ord <- read_excel("../data/Exp_results_ord.xlsx")
ord_plot <- ggplot(ord, aes(x = X, y = Y, colour = Event, shape = Site))+
  geom_point(size = 6)+
  scale_color_discrete(name = "Sampling Event", labels = c("2020", "1959"))+
  scale_shape_discrete(labels = c("Blackwoods", "Beech Mtn", "Otter Pt", "Pemetic Mtn"))+
  labs(x = "Axis 1", y = "Axis 2")+
  theme_bw(base_size = 18)+
  expand_limits(y=c(5,10))+
  expand_limits(x=c(5,10))


ord_plot
ggsave("./figures/exp_ord2.jpg", ord_plot)

