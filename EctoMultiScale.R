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
  left_join(., ecto.clean, by = "SampleNo.") %>%
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

# Prepping data for MSOM ----------------------
# Clean data
mecto.clean <- mamm.ecto %>%
  filter(!(SampleNo. != "" & is.na(Species) == T)) %>%
  filter(Tag != "RECAP")

# Select columns and add occupancy column
mecto.smol <- mecto.clean %>%
  select(Tag, Site, Day, Genus, Species) %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
  mutate(Ecto = case_when(Ecto == "_" ~ NA_character_,
                          Ecto == "" ~ NA_character_,
                          TRUE ~ Ecto)) %>%
  mutate(Occ = case_when(is.na(Ecto) == T ~ 0, TRUE ~ 1))

# Fill in missing host/ecto pairs and host/day pairs
combos <- distinct(expand_grid(mecto.smol$Tag, mecto.smol$Ecto))

mecto.joined <- mecto.smol %>%
  full_join(combos, by = c("Tag" = "mecto.smol$Tag", 
                           "Ecto" ="mecto.smol$Ecto")) %>%
  filter(is.na(Ecto) == F)
  
combos2 <- mecto.joined %>%
  group_by(Tag) %>%
  expand(Tag, Day, Ecto)
  
mecto.full <- mecto.joined %>%
  right_join(combos2) %>%
  mutate(Occ = case_when(is.na(Occ) == T ~ 0, TRUE ~ Occ)) %>%
  arrange(Tag) %>%
  fill(Site, Day) %>%
  distinct()

# Expand trap days
mecto.days <- mecto.full %>%
  group_by(Tag) %>%
  pivot_wider(names_from = Day, values_from = Occ, values_fn = sum) %>%
  select(Tag, Site, Ecto, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`) %>%
  mutate_at(.vars = '4':'12', .funs = function(x) case_when(x > 1 ~ 1,
                                                            x == 1 ~ 1,
                                                            x == 0 ~ 0)) %>%
  ungroup()

# Shift non-0 values left to convert to capture no.
caplist <- list()
for(i in 1:nrow(mecto.days)){
  x <- mecto.days[i,4:12]
  y <- c(x[!is.na(x)], x[is.na(x)])
  y <- data.frame(as.list(y))
  colnames(y) <- 1:9

  z <- cbind(mecto.days[i,1:3], y)
  caplist[[i]] <- z
}

mecto.caps <- do.call("rbind", caplist)
mecto.caps <- mecto.caps[,-c(9:12)]

# Convert to list
ecto.list <- split(mecto.caps, mecto.caps$Ecto)

# Create vector for site and drop from list
site.vec <- as.numeric(factor(ecto.list[[1]]$Site))

ecto.list <- lapply(ecto.list, 
                    function(x) x <- select(x, -c(Site, Ecto, Tag)))

# Coerce to array
ecto.ar <- array(as.numeric(unlist(ecto.list)), 
                 dim=c(nrow(ecto.list[[1]]), ncol(ecto.list[[1]]), 
                       length(ecto.list)))

# Write model -------------------------
cat("
    model{
    # Define hyperpriors
    
    # Priors
    for(i in 1:necto){
    
      a0[i]
      
      b0[i]
      
      c0[i]
    
      # Occupancy model: Site
      for(j in 1:nsite){
      
        logit(psi[j,i]) <- a0[i]
        Z[j,i] ~ dbern(psi[j,i])
    
        # Occupancy model: Host
        for(k in 1:nind){
        
          logit(theta[k[site.vec==j], i]) <- b0[i]
          mu.theta[k[site.vec==j], i] <- theta[k[site.vec==j],i]*Z[j,i]
          Y[k[site.vec==j],i] ~ dbern(mu.theta[k[site.vec==j], i])
    
          # Detection model
          for(l in 1:ncap){
          
          logit(p[k[site.vec==j],l,i]) <- c0[i]
          mu.p[k[site.vec==j],l,i] <- p[k[site.vec==j],l,i]*
                                      theta[k[site.vec==j],i]
          obs[k[site.vec==j],l,i] ~ dbern(mu.p[k[site.vec==j],l,i])
    
          }
        }
      }
    }
    }
    ", file="ectomod.txt")