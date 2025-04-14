##########
# By: Ethan Gyllenhaal
# Updated 10 April 2025
#
# R script used for making species richness map for antbirds
# originally tried to adapt someone else's code, but it was outdated, and I remade everything
## with a lot of help from https://damariszurell.github.io/EEC-MGC/a4_RangeMaps.html
# NOTE: this script takes a lot of time and memory, up to 8G if I remember correctly
########

# not sure which of these are used in final script, sorry
library(raster)
library(sp)
library(rangeBuilder)
library(spdep)
library(ncf)
library(dismo)
library(letsR)
library(rworldmap)
library(rangeMapper)
library(sf)
library(tidyverse)
library(terra)
library(letsR) # main package unique to this

setwd("path/to/dir/")

#read in range polygons in shapefile format
# obtained from BirdLife: https://datazone.birdlife.org/contact-us/request-our-data
full_ranges <- st_read("species/BOTW_2024_2.gpkg") # ALL SPECIES, big file
# should try this function in the future: lets.presab.birds

# read in species list of Thamnophilidae
# copied from birdlife taxonomy to avoid mismatches
species <- read.csv("species.txt", header=F)

# use a join to limit all species dataframe to just antbirds
antbird_ranges <- left_join(species, full_ranges, by = join_by(V1==sci_name))
# rename column
names(antbird_ranges)[names(antbird_ranges) == 'V1'] <- 'sci_name'
# convert to sf
antbird_sf <- st_as_sf(antbird_ranges)

# letsR makes a presence-absence matric with given geography limits
ant_rich <- lets.presab(antbird_sf, resol=0.1, xmn = -92, xmx = -33, ymn = -33, ymx = 18)
# makes a pallete for plotting
ant_colors <- colorRampPalette(c("white", "slateblue4"))
plot(ant_rich)


# once I figured this out I did it for two of my study genera
# figured I'd include it as an example for anyone interested!
# (obviously not used here, but Pachycephala diversity is highest on continents on Aplonis on the Solomon Islands)
# ((Australia has a depauperate core bc it's arid, but NG's has a very species rich mountain range))

aplonis_ranges <- filter(full_ranges, grepl("Aplonis", sci_name))
names(aplonis_ranges)[names(aplonis_ranges) == 'V1'] <- 'sci_name'
aplonis_sf <- st_as_sf(aplonis_ranges)
aplonis_rich <- lets.presab(aplonis_sf, resol=0.5, xmn = 85, xmx = 180, ymn = -35, ymx = 30)
plot(aplonis_rich)

pachy_ranges <- filter(full_ranges, grepl("Pachycephala", sci_name))
names(pachy_ranges)[names(pachy_ranges) == 'V1'] <- 'sci_name'
pachy_sf <- st_as_sf(pachy_ranges)
pachy_rich <- lets.presab(pachy_sf, resol=0.5, xmn = 85, xmx = 180, ymn = -45, ymx = 25)
whistler_colors <- colorRampPalette(c("yellow1", "black"))
plot(pachy_rich, col_rich=whistler_colors)
