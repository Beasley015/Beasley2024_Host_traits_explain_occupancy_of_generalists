###################################################
# Building multi-scale MSOM to test links between #
# ectoparasite life history and relative          #
# importance of covariates, host specificity      #
# Ch. 2 of dissertation                           #
# Spring 2021                                     #
###################################################

# Setup --------------------------------------------
# Load packages
library(tidyverse)

# Load data
mamm.raw <- read.csv("MammRawData.csv")
ecto.raw <- read.csv("EctosAll.csv")
veg.raw <- read.csv("VegRawData.csv")

# Set seed
set.seed(15)

# Clean data ------------------------
# Mammals
mamm.2020 <- mamm.raw %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
  filter(Date > as.Date("2020-01-01")) %>%
  select(Site:Day, Tag, Abbrev:Mass, Ecto) %>%
  filter(Tag != "") %>%
  filter(!Abbrev %in% c("", "PE??")) %>%
  rename("SampleNo." = Ecto)

# Veg
veg.2020 <- veg.raw %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
  filter(Date > as.Date("2020-01-01"))

# Merge mammal and parasite data ----------------------------
mamm.ecto <- mamm.2020 %>%
  left_join(., ecto.raw, by = "SampleNo.") %>%
  filter(Order != "")

# Mammal summary stats and early tables/figures ------------------
# Number of mammal species
length(unique(mamm.ecto$Abbrev))

# Individuals per mammal species
mamm.abund <- mamm.2020 %>%
  distinct_at(.vars = vars(Site, Abbrev, Tag), .keep_all = T) %>%
  group_by(Abbrev) %>%
  summarise(Count = n()) 
# Create rank-abundance figure
mamm.abund$Abbrev <- reorder(mamm.abund$Abbrev, -mamm.abund$Count)

ggplot(data = mamm.abund, aes(x = Abbrev, y = Count))+
  geom_col(fill = "light gray", color = "black")+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0, max(mamm.abund$Count+10)))+
  labs(x = "Species", y = "Abundance")+
  theme_bw()+
  theme(panel.grid = element_blank())

# Parasite summary and tables/figures -----------------------
# Number of ectoparasite orders
length(unique(mamm.ecto$Order))

# Number of ecto families
length(unique(mamm.ecto$Family))

# Number of ecto species
nrow(unique(mamm.ecto[c("Order", "Genus", "Species")]))-2
# Need to subtract 2 because I can't figure out how to get 
# a couple of non-ID'd Ixodidae out of there

# Prevalence by host species

# Avg. load per host species

# Raster of parasite specs present on host specs

