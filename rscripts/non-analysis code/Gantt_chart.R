#Attempting a Gantt chart to describe the project management

install.packages("remotes")
remotes::install_github("giocomai/ganttrify")
library(ganttrify)
library(tidyverse)

gc_data <- read.csv("../data/Gantt_chart.csv")
gc_data <- rename(gc_data, wp = ï..wp )# I have no idea why it's coming it weird

ganttrify(project = gc_data,
          project_start_date = "2020-01",
          font_family = "Roboto Condensed")

ggsave