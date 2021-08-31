###################################################
# Building multi-scale MSOM to test links between #
# ectoparasite life history and relative          #
# importance of covariates, host specificity      #
# Ch. 3 of dissertation                           #
# Spring 2021                                     #
###################################################

# Setup --------------------------------------------
# Load packages
library(tidyverse)
library(viridis)
library(abind)

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
  filter(!Abbrev %in% c("", "PE??", "SOCI")) %>%
  rename("SampleNo." = Ecto)

# Ectos 
ecto.clean <- ecto.raw %>%
  filter(Other != "Springtail?") %>%
  .[!(.$Order == "Ixodida" & .$Species == ""),]

# Veg
veg.2020 <- veg.raw %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
  filter(Date > as.Date("2020-01-01"))

# Merge datasets and adjust day by sampling period --------------------
breaks <- c(as.Date("2020-05-25"), as.Date("2020-06-17"),
            as.Date("2020-07-12"), as.Date("2020-08-08"))

mamm.ecto <- mamm.2020 %>%
  left_join(., ecto.raw, by = "SampleNo.") %>%
  mutate(trap.sesh = cut(as.Date(Date, format = "%Y-%m-%d"), breaks, 
                         labels = as.character(1:3))) %>%
  mutate(Day = case_when(trap.sesh == "2" ~ as.character(Day+3),
                         trap.sesh == "3" ~ as.character(Day+6),
                         TRUE ~ as.character(Day))) %>%
  mutate(Day = as.numeric(Day))

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
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("MammAbund.jpeg")

# Parasite summary and tables/figures -----------------------
# Number of ectoparasite orders
length(unique(ecto.clean$Order))

# Number of ecto families
length(unique(ecto.clean$Family))

# Number of ecto species
nrow(unique(ecto.clean[c("Order", "Genus", "Species")]))

# Abundance per order
mamm.ecto %>%
  filter(is.na(Order) == F & Order != "") %>%
  group_by(Order) %>%
  summarise(Count = n()) 
# So many fleas

# table <- mamm.ecto %>%
#   filter(is.na(Order) == F & Order != "") %>%
#   group_by(Order, Genus, Species) %>%
#   summarise(Count = n())

# Prevalence by host species
mamm.prev <- mamm.ecto %>%
  group_by(Abbrev) %>%
  summarise(Prevalence = 
              length(Order[Order != "" & is.na(Order) == F])/
              length(Order)) 

mamm.prev$Abbrev <- reorder(mamm.prev$Abbrev, 
                            -mamm.prev$Prevalence)

ggplot(data = mamm.prev, aes(x = Abbrev, y = Prevalence))+
  geom_col(color = "black", fill = "lightgray")+
  labs(x = "Species")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "MammPrev.jpeg")

# Avg. load per capture event per species
mamm.load <- mamm.ecto %>%
  group_by(Abbrev, SampleNo.) %>%
  filter(SampleNo. != "") %>%
  summarise(ecto = n())

ggplot(data = mamm.load, aes(x = ecto))+
  geom_bar(color = "black", fill = "lightgray")+
  labs(x = "Parasite Load", y = "Number of Samples")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,10,1))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 200))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("EctoLoad.jpeg")

# Raster of parasite specs present on host specs
mamm.raster <- mamm.ecto %>%
  filter(SampleNo. != "") %>%
  select(Abbrev, Order, Family, Genus, Species) %>%
  group_by(Abbrev) %>%
  distinct() %>% 
  filter(Order != "" & is.na(Order)==F) %>%
  .[!(.$Order == "Ixodida" & .$Species == ""),] %>%
  unite(ecto, Genus, Species) %>%
  mutate(Occ = 1)

ggplot(mamm.raster, aes(x = Abbrev, y = ecto))+
  geom_raster(aes(fill = Order))+
  scale_fill_viridis_d(limits = c("Siphonaptera", "Ixodida",
                                  "Mesostigmata", "Diptera"))+
  labs(x = "Host Species", y = "Parasite Species")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("ectosraster.jpeg", height = 4, width = 8.5,
#        units = "in")

# Hosts per ecto species
n.host <- mamm.ecto %>%
  filter(SampleNo. != "") %>%
  select(Abbrev, Order, Genus, Species) %>%
  distinct() %>% 
  filter(Order != "" & is.na(Order)==F) %>%
  .[!(.$Order == "Ixodida" & .$Species == ""),] %>%
  unite(ecto, Genus, Species) %>%
  group_by(Order, ecto) %>%
  summarise(hosts = n())

ggplot(data = n.host, aes(x = hosts, fill = Order))+
  geom_bar(color = "black")+
  scale_fill_manual(values = gray.colors(n = 4, end = 0.85,
                                         start = 0.05))+
  labs(x = "Number of Hosts", y = "Number of Parasite Species")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("NumOfHosts.jpeg")

# Prepping data for MSOM ----------------------------
# Clean data
mecto.clean <- mamm.ecto %>%
  filter(!(SampleNo. != "" & is.na(Species) == T))

# Select columns and add occupancy column
mecto.smol <- mecto.clean %>%
  select(Tag, Site, Day, Genus, Species) %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
  mutate(Ecto = case_when(Ecto == "_" ~ NA_character_,
                          Ecto == "" ~ NA_character_,
                          TRUE ~ Ecto)) %>%
  mutate(Occ = case_when(Ecto == NA ~ NA_real_, TRUE ~ 1))



