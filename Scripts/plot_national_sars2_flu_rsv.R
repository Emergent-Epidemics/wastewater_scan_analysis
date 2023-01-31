#Plot US data on SARS-CoV-2, influenza, RSV from wastewater scan
#SV Scarpino
#Jan. 2023

###########
#libraries#
###########
library(zoo)
library(ggplot2)
library(wesanderson)
library(dplyr)

######
#Data#
######
sars2_dat <- read.csv("~/desktop/WWSCAN_SARS-CoV-2_N_Gene_-_all_variants_Multiple_Locations_20230130.csv")

sars2_dat$collection_date <- as.POSIXct(strptime(sars2_dat$collection_date, format = "%Y-%m-%d"))

rsv_dat <- read.csv("~/desktop/WWSCAN_Respiratory_syncytial_virus_(RSV)_RSV_Multiple_Locations_20230130.csv")

rsv_dat$collection_date <- as.POSIXct(strptime(rsv_dat$collection_date, format = "%Y-%m-%d"))

flu_dat <- read.csv("~/desktop/WWSCAN_Influenza_Influenza_A_Multiple_Locations_20230130.csv")

flu_dat$collection_date <- as.POSIXct(strptime(flu_dat$collection_date, format = "%Y-%m-%d"))

##########
#Analysis#
##########
sars2_dat$pop_norm <- sars2_dat$N_Gene_gc_g_dry_weight_pmmov/sars2_dat$sewershed_pop

sars2 <- by(sars2_dat$pop_norm, sars2_dat$collection_date, median, na.rm = TRUE)

rsv_dat$pop_norm <- as.numeric(rsv_dat$RSV_gc_g_dry_weight_pmmov)/rsv_dat$sewershed_pop

rsv <- by(rsv_dat$pop_norm, rsv_dat$collection_date, median, na.rm = TRUE)

flu_dat$pop_norm <- as.numeric(flu_dat$Influenza_A_gc_g_dry_weight_pmmov)/rsv_dat$sewershed_pop

flu <- by(flu_dat$pop_norm, flu_dat$collection_date, median, na.rm = TRUE)

##########
#Plotting#
##########
disease <- c(rep("SARS-CoV-2", length(sars2)), rep("RSV", length(rsv)), rep("influenza", length(flu)))
dates <- c(names(sars2), names(rsv), names(flu))
prev <- c(rollmean(as.numeric(sars2), k = 7, na.pad = TRUE), rollmean(as.numeric(rsv), k = 7, na.pad = TRUE), rollmean(as.numeric(flu), k = 7, na.pad = TRUE))

cols <- wes_palette(name = "FantasticFox1", n = 3)
dat.plot <- data.frame(dates, disease, prev)
use <- which(as.Date(dates) > as.Date("2022-09-01"))
dat.plot.use <- dat.plot[use,]
dat.plot.use <- dat.plot.use %>%
  group_by(disease) %>%
  mutate(prev = prev/max(prev, na.rm = TRUE))

ggplot(dat.plot.use, aes(x = as.Date(dates), y = prev, color = disease, group = disease)) + geom_line(linewidth = 1) + scale_color_manual(values = cols, name = "Pathogen") + xlab("Date (2022-23)") + ylab("Wastewater signal (normalized viral concentration)") + theme(legend.position = c(0.3, 0.8), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffff75", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 15), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) 