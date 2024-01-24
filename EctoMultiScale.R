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
library(vegan)
library(R2jags)
library(multiApply)
library(agricolae)
library(gridExtra)
library(pscl)
library(spData)
library(grid)
library(rstatix)
library(patchwork)
library(sf)
library(terra)

# Load data
mamm.raw <- read.csv("MammRawData.csv")
ecto.raw <- read.csv("EctosAll.csv")
veg.raw <- read.csv("VegRawData.csv")

# Set seed
set.seed(15)

# Clean data ------------------------
# Mammals
for(i in 1:length(mamm.raw$Date)){
  if(nchar(mamm.raw$Date[i]) < 10){
    mamm.raw$Date[i] <- paste("0", mamm.raw$Date[i], sep = "")
  }
}

mamm.2020 <- mamm.raw %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
  filter(Date > as.Date("2020-01-01")) %>%
  select(Site:Day, Tag:Mass, Ecto, Desc.) %>%
  filter(Tag != "") %>%
  filter(Desc. == "") %>%
  filter(!Abbrev %in% c("", "PE??", "SOCI")) %>%
  rename("SampleNo." = Ecto) %>%
  select(!Desc.) %>%
  mutate(Genus = paste(substr(Genus, start = 1, stop = 1), ".")) %>%
  unite(col = HostName, Genus, Species, sep = " ") 

# Ectos 
ecto.clean <- ecto.raw %>%
  filter(Other != "Springtail?") %>%
  .[!(.$Order == "Ixodida" & .$Species == ""),] #%>%
  # rename(Ecto.genus = Genus)

# Veg
veg.2020 <- veg.raw %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
  filter(Date > as.Date("2020-01-01"))

# Merge datasets and adjust day by sampling period ----------------
breaks <- c(as.Date("2020-05-25"), as.Date("2020-06-17"),
            as.Date("2020-07-12"), as.Date("2020-08-08"))

mamm.ecto <- mamm.2020 %>%
  left_join(., ecto.clean, by = "SampleNo.") %>%
  mutate(trap.sesh = cut(as.Date(Date, format = "%Y-%m-%d"), breaks, 
                         labels = as.character(1:3))) %>%
  mutate(Day = case_when(trap.sesh == "2" ~ as.character(Day+3),
                         trap.sesh == "3" ~ as.character(Day+6),
                         TRUE ~ as.character(Day))) %>%
  mutate(Day = as.numeric(Day)) %>%
  # Run this line to keep Mesostigmata aggregated:
  mutate(Genus = case_when(Order == "Psocodea" ~ "Unknown Psocodea",
                           Order == "Mesostigmata" ~
                             "Unknown Mesostigmata",
                           TRUE ~ Genus))
  #Run this line to disaggregate Mesostigmata
  # mutate(Genus = case_when(Order == "Psocodea" ~ "Unknown Psocodea",
  #                          Order == "Mesostigmata" & Genus == "" ~
  #                            "Unknown Mesostigmata",
  #                          TRUE ~ Genus))

# Mammal summary stats and early tables/figures ------------------
# Number of mammal species
length(unique(mamm.ecto$Abbrev))

# Individuals per mammal species
mamm.abund <- mamm.2020 %>%
  distinct_at(.vars = vars(Site, Abbrev, Tag), .keep_all = T) %>%
  group_by(HostName) %>%
  summarise(Count = n()) 

# Create rank-abundance figure
mamm.abund$SpecName <- reorder(mamm.abund$HostName, -mamm.abund$Count)

ggplot(data = mamm.abund, aes(x = HostName, y = Count))+
  geom_col(fill = "light gray", color = "black")+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0, max(mamm.abund$Count+10)))+
  labs(x = "Species", y = "Abundance")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1,
                                   hjust = 1))

# ggsave("MammAbund.jpeg", width = 4, height = 3, units = "in")

# Parasite summary and tables/figures -----------------------
# Number of ectoparasite orders
length(unique(ecto.clean$Order))

# Number of ecto families
length(unique(ecto.clean$Family))

# Number of ecto species
nrow(unique(ecto.clean[c("Order", "Genus", "Species")]))

# Abundance per order
ecto.abund <- mamm.ecto %>%
  filter(is.na(Order) == F & Order != "") %>%
  group_by(Order) %>%
  summarise(Count = n()) 
# So many fleas

ggplot(data = ecto.abund, aes(x = Order, y = Count))+
  geom_col()+
  labs(y = "Abundance")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

# ggsave(filename = "ectoabund.jpg", width = 7, height = 4.5,
#        units = "in")

# table <- mamm.ecto %>%
#   filter(is.na(Order) == F & Order != "") %>%
#   group_by(Order, Genus, Species) %>%
#   summarise(Count = n())

# Prevalence by host species
mamm.prev <- mamm.ecto %>%
  group_by(HostName) %>%
  summarise(Prevalence = 
              length(Order[Order != "" & is.na(Order) == F])/
              length(Order)) 

mamm.prev$HostName <- reorder(mamm.prev$HostName, 
                            -mamm.prev$Prevalence)

ggplot(data = mamm.prev, aes(x = HostName, y = Prevalence))+
  geom_col(color = "black", fill = "lightgray")+
  labs(x = "Species")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1,
                                   hjust = 1))

# ggsave(filename = "MammPrev.jpeg", width = 4, height = 3, 
#        units = "in")

# Avg. load per capture event per species
mamm.load <- mamm.ecto %>%
  group_by(HostName, SampleNo.) %>%
  filter(SampleNo. != "") %>%
  summarise(ecto = n())

ggplot(data = mamm.load, aes(x = ecto))+
  geom_bar(color = "black", fill = "lightgray")+
  labs(x = "Parasite Load", y = "Number of Samples")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,10,1))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 200))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("EctoLoad.jpeg", width = 4, height = 3, units = "in")

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
  # Comment the following to aggregate I.scapularis:
  # mutate(Species = case_when(Genus == "Ixodes" &
  #                              Species == "scapularis" ~
  #                              paste(Species, Other, sep = "_"),
  #                            TRUE ~ Species)) %>%
  select(Abbrev, Order, Family, Genus, Species) %>%
  distinct() %>% 
  filter(Order != "" & is.na(Order)==F) %>%
  unite(ecto, Genus, Species) %>%
  group_by(Order, Family, ecto) %>%
  summarise(hosts = n())

ggplot(data = n.host, aes(x = hosts, fill = Order))+
  geom_bar(color = "black")+
  scale_fill_manual(values = gray.colors(n = 5, end = 0.85,
                                         start = 0.05))+
  labs(x = "Number of Hosts", y = "Number of Parasite Species")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("NumOfHosts.jpeg")

# by host subfamily:
n.host.sub <- mamm.ecto %>%
  filter(SampleNo. != "") %>%
  mutate(Species = case_when(Genus == "Ixodes" & 
                               Species == "scapularis" ~
                               paste(Species, Other, sep = "_"),
                             TRUE ~ Species)) %>%
  select(Abbrev, Order, Family, Genus, Species) %>%
  mutate(HostSub = case_when(Abbrev %in% c("ZAHU", "NAIN")
                             ~ "Zapodinae",
                             Abbrev %in% c("PELE", "PEMA")
                             ~ "Neotominae",
                             Abbrev %in% c("MIPE", "MYGA")
                             ~ "Arvicolinae",
                             Abbrev == "TAST" ~ "Xerinae",
                             Abbrev == "BLBR" ~ "Soricinae")) %>%
  select(-Abbrev) %>%
  distinct() %>%
  filter(Order != "" & is.na(Order)==F) %>%
  unite(ecto, Genus, Species) %>%
  group_by(Order, Family, ecto) %>%
  summarise(hosts = n())

ggplot(data = n.host.sub, aes(x = hosts, fill = Order))+
  geom_bar(color = "black")+
  scale_fill_manual(values = gray.colors(n = 5, end = 0.85,
                                         start = 0.05))+
  labs(x = "Number of Hosts", y = "Number of Parasite Species")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("hostsubfam.jpeg")

# by host family:
n.host.fam <- mamm.ecto %>%
  filter(SampleNo. != "") %>%
  mutate(Species = case_when(Genus == "Ixodes" & 
                               Species == "scapularis" ~
                               paste(Species, Other, sep = "_"),
                             TRUE ~ Species)) %>%
  select(Abbrev, Order, Family, Genus, Species) %>%
  mutate(HostSub = case_when(Abbrev %in% c("ZAHU", "NAIN")
                             ~ "Dipodidae",
                             Abbrev %in% c("PELE", "PEMA", 
                                           "MIPE", "MYGA")
                             ~ "Cricetidae",
                             Abbrev == "TAST" ~ "Sciuridae",
                             Abbrev == "BLBR" ~ "Soricidae")) %>%
  select(-Abbrev) %>%
  distinct() %>%
  filter(Order != "" & is.na(Order)==F) %>%
  unite(ecto, Genus, Species) %>%
  group_by(Order, Family, ecto) %>%
  summarise(hosts = n())

ggplot(data = n.host.fam, aes(x = hosts, fill = Order))+
  geom_bar(color = "black")+
  scale_fill_manual(values = gray.colors(n = 5, end = 0.85,
                                         start = 0.05))+
  labs(x = "Number of Hosts", y = "Number of Parasite Species")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("nhostfam.jpeg")

# Test host richness/parasite richness relationship
rich <- mamm.ecto %>%
  filter(SampleNo. != "") %>%
  select(Site, Abbrev, Order, Family, Genus, Species) %>%
  group_by(Abbrev) %>%
  distinct() %>% 
  filter(Order != "" & is.na(Order)==F) %>%
  unite(ecto, Genus, Species) %>%
  mutate(Occ = 1) %>%
  group_by(Site) %>%
  summarise(Mamm.rich = length(unique(Abbrev)), 
            ecto.rich = length(unique(ecto)))

summary(lm(data = rich, ecto.rich~Mamm.rich))

ggplot(data = rich, aes(x = Mamm.rich, y = ecto.rich))+
  geom_point(size =2)+
  # geom_smooth(method = 'lm', se = F, color = "black")+
  labs(x = "Host Richness", y = "Ectoparasite Richness")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("hostectorich.jpeg", width = 5, height = 3, units = "in")

# Prepping data for MSOM ----------------------
# Clean data
mecto.clean <- mamm.ecto %>%
  filter(!(SampleNo. != "" & is.na(Species) == T)) %>%
  # Include this line to filter "rare" hosts (<10 individuals):
  filter(!(Abbrev %in% c("NAIN", "MYGA"))) %>%
  filter(Tag != "RECAP" & Tag != "DESC")

# Select columns and add occupancy column
mecto.smol <- mecto.clean %>%
  select(Tag, Site, Day, Genus, Species, Other) %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
  #Comment this line to keep Iscap aggregated:
  mutate(Ecto = case_when(Ecto == "Ixodes_scapularis" ~
                            paste(Ecto, Other, sep = "_"),
                          TRUE ~ Ecto)) %>%
  select(-Other) %>%
  mutate(Ecto = case_when(Ecto == "_" ~ NA_character_,
                          Ecto == "" ~ NA_character_,
                          TRUE ~ Ecto)) %>%
  mutate(Occ = case_when(is.na(Ecto) == T ~ 0, TRUE ~ 1))

# Fill in missing host/ecto pairs and host/day pairs
combos <- distinct(expand_grid(mecto.smol$Tag, mecto.smol$Ecto))

mecto.joined <- mecto.smol %>%
  full_join(combos, by = c("Tag" = "mecto.smol$Tag", 
                           "Ecto" ="mecto.smol$Ecto")) 

combos2 <- mecto.joined %>%
  expand(nesting(Tag, Day), Ecto)
  
mecto.full <- mecto.joined %>%
  right_join(combos2) %>%
  mutate(Occ = case_when(is.na(Occ) == T ~ 0, TRUE ~ Occ)) %>%
  arrange(Tag) %>%
  fill(Site, Day) %>%
  distinct()

# Expand trap days
mecto.days <- mecto.full %>%
  group_by(Tag) %>%
  pivot_wider(names_from = Day, values_from = Occ, 
              values_fn = sum) %>%
  select(Tag, Site, Ecto, `1`, `2`, `3`, `4`, `5`, `6`, 
         `7`, `8`, `9`) %>%
  mutate_at(.vars = '4':'12', 
            .funs = function(x) case_when(x >= 1 ~ 1,
                                          x == 0 ~ 0)) %>%
  ungroup() %>%
  filter(!is.na(Ecto))

# Shift non-0 values left to convert to capture no.
cap.convert <- function(dat, cols, cols2){
  caplist <- list()
  for(i in 1:nrow(dat)){
    x <- dat[i,cols]
    y <- c(x[!is.na(x)], x[is.na(x)])
    y <- data.frame(as.list(y))
    colnames(y) <- 1:9

    z <- cbind(dat[i,cols2], y)
    caplist[[i]] <- z
  }
  return(caplist)
}

caplist <- cap.convert(mecto.days, 4:12, 1:3)

mecto.caps <- do.call("rbind", caplist)
mecto.caps <- mecto.caps[,-12]

# Convert to list
ecto.list <- split(mecto.caps, mecto.caps$Ecto)

# Create vector for site and drop from list
site.vec <- as.integer(factor(ecto.list[[1]]$Site))
tag.id <- 1:length(site.vec)

# Create matrix of id's per site
taglist <- list()
for(i in 1:length(unique(site.vec))){
  taglist[[i]] <- print(tag.id[site.vec==i])
}

library(plyr)
tagmat <- t(ldply(taglist, rbind))
detach('package:plyr', unload = T)

ecto.list2 <- lapply(ecto.list, 
                    function(x) x <- select(x, -c(Site, Ecto, Tag)))

# Coerce to array
ecto.ar <- array(as.numeric(unlist(ecto.list2)), 
                 dim=c(nrow(ecto.list2[[1]]), 8, #8 = max caps
                       length(ecto.list2)))

# Site covariate: Veg cover type -----------------------
veg.site <- veg.2020 %>%
  group_by(Site, Habitat) %>%
  summarise(across(.cols = c(Canopy:X.DeadVeg), mean, na.rm = T))

comp.pc <- prcomp(veg.site[,10:15])

comp.cov1 <- comp.pc$x[,1]
comp.cov1 <- as.vector(scale(comp.cov1))

comp.cov2 <- comp.pc$x[,2]
comp.cov2 <- as.vector(scale(comp.cov2))

comp.frame <- cbind(veg.site, pc1 = comp.cov1, 
                    pc2 = comp.cov2)

p <- ggplot(data = comp.frame, aes(x = pc1, y = pc2, color = Habitat))+
  geom_point(size = 2)+
  scale_color_viridis_d(end = 0.95)+
  labs(x = "PC1", y = "PC2")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())+
  annotation_custom(grob = textGrob(label = "Leaf Litter"), 
                    ymin = -2.2, ymax = -2.2, xmin = 1.4, 
                    xmax = 1.4)+
  annotation_custom(grob = textGrob(label = "Grass/Forb"),
                    ymin = -2.2, ymax = -2.2, xmin = -0.9,
                    xmax = -0.9)+
  annotation_custom(grob = textGrob(label = "Grass",
                                    rot = 90),
                    ymin = -1.5, ymax = -1.6, xmin = -1.6,
                    xmax = -1.5)+
  annotation_custom(grob = textGrob(label = "Forb", rot = 90),
                    ymin = 1.5, ymax = 1.5, xmin = -1.6,
                    xmax = -1.6)

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

# ggsave(gt, filename = "pc1.jpeg", width = 4.5, height = 3.5, units = "in")

# Site covariate: Host Abundance ----------------
abund.frame <- mecto.clean %>%
  select(Site, Tag) %>%
  distinct() %>%
  group_by(Site) %>%
  summarise(HostAbund = n())

abund.vec <- as.vector(scale(abund.frame$HostAbund))

# Site covariate: Shannon diversity
host.div <- mecto.clean %>%
  select(Site, Abbrev, Tag) %>%
  distinct() %>%
  group_by(Site, Abbrev) %>%
  summarise(abund = n()) %>%
  summarise(shannon = diversity(abund, index = "shannon"))

div.vec <- as.vector(scale(host.div$shannon))

# Host covariate: Species -----------------------
# Get species abbrev of each host
host.spec <- mecto.clean %>%
  select(Tag, Abbrev) %>%
  distinct() %>%
  arrange(Tag)

# Convert to factor
hostspec <- factor(host.spec$Abbrev, 
                   levels = unique(host.spec$Abbrev))

# Host covariate: Adj. body mass -----------------------
host.mass <- mecto.clean %>%
  select(Tag, Abbrev, Mass) %>%
  arrange(Tag) %>%
  group_by(Tag, Abbrev) %>%
  summarise(Mass = mean(Mass, na.rm = T)) %>%
  ungroup()

host.z <- host.mass %>%
  group_by(Abbrev) %>%
  mutate(mass.scaled = as.vector(scale(Mass)))

hostmass <- host.z$mass.scaled

# Host covariate: Sex -------------------------
host.sex <- mecto.clean %>%
  select(Tag, Abbrev, Sex) %>%
  arrange(Tag) %>%
  group_by (Tag, Abbrev) %>%
  distinct()

hostsex <- host.sex$Sex
hostsex[hostsex == ""] <- NA
hostsex <- as.numeric(factor(hostsex))-1

# Detection covariate: Julian date -----------------------
dates.frame <- mecto.clean %>%
  mutate(julian = format(Date, "%j")) %>%
  mutate(julian = as.vector(scale(as.numeric(julian)))) %>%
  dplyr::select(Tag, Day, julian) %>%
  distinct() %>%
  pivot_wider(names_from = Day, values_from = julian) %>%
  arrange(Tag) %>%
  relocate(Tag, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`) %>%
  filter(Tag %in% mecto.caps$Tag)

dateslist <- cap.convert(dates.frame, 2:10, 1)

dates <- do.call("rbind", dateslist)
dates <- dates[,c(-1, -10)]
  
# Write model -------------------------
cat("
    model{
    # Define hyperpriors
    
    mean.a0 ~ dunif(0,1)
    a0.mean <- log(mean.a0)-log(1-mean.a0)
    tau.a0 ~ dgamma(0.1, 0.1)
    
    mean.a1 ~ dunif(0,1)
    a1.mean <- log(mean.a1)-log(1-mean.a1)
    tau.a1 ~ dgamma(0.1, 0.1)
    
    mean.a2 ~ dunif(0,1)
    a2.mean <- log(mean.a2)-log(1-mean.a2)
    tau.a2 ~ dgamma(0.1, 0.1)
    
    mean.a3 ~ dunif(0,1)
    a3.mean <- log(mean.a3)-log(1-mean.a3)
    tau.a3 ~ dgamma(0.1, 0.1)
    
    mean.a4 ~ dunif(0,1)
    a4.mean <- log(mean.a4)-log(1-mean.a4)
    tau.a4 ~ dgamma(0.1, 0.1)
    
    mean.b0 ~ dunif(0,1)
    b0.mean <- log(mean.b0)-log(1-mean.b0)
    tau.b0 ~ dgamma(0.1, 0.1)
    
    mean.b2 ~ dunif(0,1)
    b2.mean <- log(mean.b2)-log(1-mean.b2)
    tau.b2 ~ dgamma(0.1, 0.1)
    
    mean.b3 ~ dunif(0,1)
    b3.mean <- log(mean.b3)-log(1-mean.b3)
    tau.b3 ~ dgamma(0.1, 0.1)
    
    mean.c0 ~ dunif(0,1)
    c0.mean <- log(mean.c0)-log(1-mean.c0)
    tau.c0 ~ dgamma(0.1, 0.1)
    
    mean.c1 ~ dunif(0,1)
    c1.mean <- log(mean.c1)-log(1-mean.c1)
    tau.c1 ~ dgamma(0.1, 0.1)
    
    mean.c2 ~ dunif(0,1)
    c2.mean <- log(mean.c2)-log(1-mean.c2)
    tau.c2 ~ dgamma(0.1, 0.1)
    
    probs ~ dunif(0,1)
    
    # Priors
    for(i in 1:necto){
    
      a0[i] ~ dnorm(a0.mean, tau.a0)
      a1[i] ~ dnorm(a1.mean, tau.a1)
      a2[i] ~ dnorm(a2.mean, tau.a2)
      a3[i] ~ dnorm(a3.mean, tau.a3)
      a4[i] ~ dnorm(a4.mean, tau.a4)
      
      b0[i] ~ dnorm(b0.mean, tau.b0)
      b1[i,1] <- 0
      for(s in 2:K){
        b1[i,s] ~ dnorm(0, 1)
      }
      b2[i] ~ dnorm(b2.mean, tau.b2)
      b3[i] ~ dnorm(b3.mean, tau.b3)
      b4[i,1] <- 0
      for(s in 2:K){
        b4[i,s] ~ dnorm(0,1)
      }
      
      c0[i] ~ dnorm(c0.mean, tau.c0)
      c1[i] ~ dnorm(c1.mean, tau.c1)
      c2[i] ~ dnorm(c2.mean, tau.c2)
    
      # Occupancy model: Site
      for(j in 1:nsite){
      
        logit(psi[j,i]) <- a0[i] + a1[i]*deadveg[j] + 
                            a2[i]*grassforb[j] + a3[i]*abund[j]
                            + a4[i]*div[j]
        Z[j,i] ~ dbern(psi[j,i])

        # Occupancy model: Host
        for(k in tagmat[1:hostvec[j],j]){
          # May need to create matrix of host sex
          sex[k,i] ~ dbern(probs)
        
          logit(theta[k,i]) <- b0[i] + b1[i,host[k]] + b2[i]*mass[k]
                                      + b3[i]*sex[k,i] 
                                      + b4[i,host[k]]*deadveg[j]
          Y[k,i] ~ dbern(theta[k,i]*Z[j,i])
          
          # Detection model
          for(l in 1:ncap[k]){
          
            logit(p[k,l,i]) <- c0[i] + c1[i]*l + c2[i]*julian[k,l]
            obs[k,l,i] ~ dbern(p[k,l,i]*Y[k,i])
    
          }
        }
      }
    }
    }
    ", file="ectomod.txt")

# Set up model --------------------------
# Define no. of ectos, sites, hosts, etc.
necto <- length(unique(mecto.caps$Ecto))

nsite <- length(unique(site.vec))

ncap <- apply(ecto.ar[,,1], 1, function(x) length(na.omit(x)))

site.occ <- mecto.caps %>%
  dplyr::select(-Tag) %>%
  rowwise() %>%
  mutate(Occ = max(c_across(`1`:`8`), na.rm = T)) %>%
  dplyr::select(Site, Ecto, Occ) %>%
  distinct() %>%
  group_by(Site, Ecto) %>%
  summarise(Occ = max(Occ)) %>%
  arrange(Site) %>%
  pivot_wider(names_from = Ecto, values_from = Occ)

sitemax <- as.matrix(site.occ[,-c(1)])

hostmax <- apply(ecto.ar, c(1,3), max, na.rm = T)#[,-c(1,4)]

hostvec <- apply(tagmat, 2, function(x) length(na.omit(x)))

# Define data to send and params to keep
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4', 'Y', 'Z')

# Init values
inits <- function(){
  inits <- list(
    Z = sitemax,
    Y = hostmax
  )
}

# Send model to JAGS
# model <- jags(model.file = 'ectomod_V1V2D.txt', data = datalist,
#               n.chains = 3, parameters.to.save = params,
#               inits = inits, n.burnin = 15000, n.iter = 30000,
#               n.thin = 20)

# Save model
# saveRDS(model, file = "ectomod.rds")
# saveRDS(model, file = "mesos_agg.rds")
# saveRDS(model, file = "iscap_agg.rds")
# saveRDS(model, file = "sansrarehosts.rds")

# Load model
# model <- readRDS("ectomod.rds")
# model <- readRDS("sansabund.rds")
# model <- readRDS("mesos_agg.rds")
# model <- readRDS("iscap_agg.rds")
model <- readRDS("sansrarehosts.rds")

# Model selection: site ------------------------
# Full model
model <- jags(model.file = 'ectomod.txt', data = datalist,
              n.chains = 3, parameters.to.save = params,
              inits = inits, n.burnin = 9000, n.iter = 20000,
              n.thin = 10)

samples.site.full <- jags.samples(model$model, 
                             variable.names = c("WAIC"),
                             n.iter = 5000, n.burnin = 1000, 
                             n.thin = 1, type = "mean")

# Intercept model
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates)

params <- c('b1', 'b2', 'b3', 'b4', 'c1', 'c2', 'psi','theta', 
            'Z', 'Y')

model.site.intercept <- jags(model.file = 'ectomod_SiteIntercept.txt',
                             data = datalist, n.chains = 3, 
                             parameters.to.save = params,
                             inits = inits, n.burnin = 9000, 
                             n.iter = 20000, n.thin = 10)

samples.site.intercept <- jags.samples(model.site.intercept$model, 
                                       variable.names = c("WAIC"),
                                       n.iter = 5000, n.burnin = 1000,
                                       n.thin = 1, type = "mean")

# Veg 1 only
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates)

params <- c('a1', 'b1', 'b2', 'b3', 'b4', 'c1', 'c2', 'psi',
            'theta', 'Z', 'Y')

model.site.deadveg <- jags(model.file = 'ectomod_Veg1.txt',
                             data = datalist, n.chains = 3, 
                             parameters.to.save = params,
                             inits = inits, n.burnin = 9000, 
                             n.iter = 20000, n.thin = 10)

samples.site.deadveg <- jags.samples(model.site.deadveg$model, 
                                       variable.names = c("WAIC"),
                                       n.iter = 5000, n.burnin = 1000,
                                       n.thin = 1, type = "mean")

# Veg 2 only
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, grassforb = comp.cov2,
                 deadveg = comp.cov1,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates)

params <- c('a2', 'b1', 'b2', 'b3', 'b4', 'c1', 'c2', 'psi',
            'theta', 'Z', 'Y')

model.site.grassforb <- jags(model.file = 'ectomod_Veg2.txt',
                           data = datalist, n.chains = 3, 
                           parameters.to.save = params,
                           inits = inits, n.burnin = 9000, 
                           n.iter = 20000, n.thin = 10)

samples.site.grassforb <- jags.samples(model.site.grassforb$model, 
                                     variable.names = c("WAIC"),
                                     n.iter = 5000, n.burnin = 1000,
                                     n.thin = 1, type = "mean")

# Host Abundance
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec)

params <- c('a3', 'b1', 'b2', 'b3', 'b4', 'c1', 'c2',
            'psi', 'theta')

model.site.hostabund <- jags(model.file = 'ectomod_Hostabund.txt',
                             data = datalist, n.chains = 3, 
                             parameters.to.save = params,
                             inits = inits, n.burnin = 9000, 
                             n.iter = 20000, n.thin = 10)

samples.site.hostabund <- jags.samples(model.site.hostabund$model, 
                                       variable.names = c("WAIC"),
                                       n.iter = 5000, n.burnin = 1000,
                                       n.thin = 1, type = "mean")

# Host Diversity
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, div = div.vec)

params <- c('a4', 'b1', 'b2', 'b3', 'b4', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.site.hostdiv <- jags(model.file = 'ectomod_Hostdiv.txt',
                             data = datalist, n.chains = 3, 
                             parameters.to.save = params,
                             inits = inits, n.burnin = 9000, 
                             n.iter = 20000, n.thin = 10)

samples.site.hostdiv <- jags.samples(model.site.hostdiv$model, 
                                       variable.names = c("WAIC"),
                                       n.iter = 5000, n.burnin = 1000,
                                       n.thin = 1, type = "mean")

# Veg1 + Veg2
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates)

params <- c('a1', 'a2', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.site.veggies <- jags(model.file = 'ectomod_Veg.txt',
                           data = datalist, n.chains = 3, 
                           parameters.to.save = params,
                           inits = inits, n.burnin = 9000, 
                           n.iter = 20000, n.thin = 10)

samples.site.veggies <- jags.samples(model.site.veggies$model, 
                                     variable.names = c("WAIC"),
                                     n.iter = 5000, n.burnin = 1000,
                                     n.thin = 1, type = "mean")

# Veg1 + Abundance
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec)

params <- c('a1', 'a3', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.site.veg1abund <- jags(model.file = 'ectomod_Veg1Abund.txt',
                           data = datalist, n.chains = 3, 
                           parameters.to.save = params,
                           inits = inits, n.burnin = 9000, 
                           n.iter = 20000, n.thin = 10)

samples.site.veg1abund <- jags.samples(model.site.veg1abund$model, 
                                     variable.names = c("WAIC"),
                                     n.iter = 5000, n.burnin = 1000,
                                     n.thin = 1, type = "mean")

# Veg1 + Diversity
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, div = div.vec)

params <- c('a1', 'a4', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.site.veg1div <- jags(model.file = 'ectomod_Veg1Div.txt',
                             data = datalist, n.chains = 3, 
                             parameters.to.save = params,
                             inits = inits, n.burnin = 9000, 
                             n.iter = 20000, n.thin = 10)

samples.site.veg1div <- jags.samples(model.site.veg1div$model, 
                                       variable.names = c("WAIC"),
                                       n.iter = 5000, n.burnin = 1000,
                                       n.thin = 1, type = "mean")

# Veg2 + Abundance
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec)

params <- c('a2', 'a3', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.site.veg2abund <- jags(model.file = 'ectomod_Veg2abund.txt',
                           data = datalist, n.chains = 3, 
                           parameters.to.save = params,
                           inits = inits, n.burnin = 9000, 
                           n.iter = 20000, n.thin = 10)

samples.site.veg2abund <- jags.samples(model.site.veg2abund$model, 
                                     variable.names = c("WAIC"),
                                     n.iter = 5000, n.burnin = 1000,
                                     n.thin = 1, type = "mean")

# Veg2 + Diversity
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, div = div.vec)

params <- c('a2', 'a4', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.site.veg2div <- jags(model.file = 'ectomod_Veg2div.txt',
                             data = datalist, n.chains = 3, 
                             parameters.to.save = params,
                             inits = inits, n.burnin = 9000, 
                             n.iter = 20000, n.thin = 10)

samples.site.veg2div <- jags.samples(model.site.veg2div$model, 
                                       variable.names = c("WAIC"),
                                       n.iter = 5000, n.burnin = 1000,
                                       n.thin = 1, type = "mean")

# Abundance + Diversity
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a3', 'a4', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.site.abunddiv <- jags(model.file = 'ectomod_AbundDiv.txt',
                           data = datalist, n.chains = 3, 
                           parameters.to.save = params,
                           inits = inits, n.burnin = 9000, 
                           n.iter = 20000, n.thin = 10)

samples.site.abunddiv <- jags.samples(model.site.abunddiv$model, 
                                     variable.names = c("WAIC"),
                                     n.iter = 5000, n.burnin = 1000,
                                     n.thin = 1, type = "mean")

# Veg1 + Veg2 + Abundance
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec)

params <- c('a1', 'a2', 'a3', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.site.v1v2a <- jags(model.file = 'ectomod_V1V2A.txt',
                            data = datalist, n.chains = 3, 
                            parameters.to.save = params,
                            inits = inits, n.burnin = 9000, 
                            n.iter = 20000, n.thin = 10)

samples.site.v1v2a <- jags.samples(model.site.v1v2a$model, 
                                      variable.names = c("WAIC"),
                                      n.iter = 5000, n.burnin = 1000,
                                      n.thin = 1, type = "mean")

# Veg1 + Veg2 + Diversity
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, div = div.vec)

params <- c('a1', 'a2', 'a4', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.site.v1v2d <- jags(model.file = 'ectomod_V1V2D.txt',
                         data = datalist, n.chains = 3, 
                         parameters.to.save = params,
                         inits = inits, n.burnin = 15000, 
                         n.iter = 30000, n.thin = 20)

# saveRDS(model.site.v1v2d, file="sansabund.rds")

samples.site.v1v2d <- jags.samples(model.site.v1v2d$model, 
                                   variable.names = c("WAIC"),
                                   n.iter = 5000, n.burnin = 1000,
                                   n.thin = 1, type = "mean")

# Veg1 + Abundance + Diversity
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a3', 'a4', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.site.v1ad <- jags(model.file = 'ectomod_v1aD.txt',
                         data = datalist, n.chains = 3, 
                         parameters.to.save = params,
                         inits = inits, n.burnin = 9000, 
                         n.iter = 20000, n.thin = 10)

samples.site.v1ad <- jags.samples(model.site.v1ad$model, 
                                   variable.names = c("WAIC"),
                                   n.iter = 5000, n.burnin = 1000,
                                   n.thin = 1, type = "mean")

# Veg2 + Abundance + Diversity
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a2', 'a3', 'a4', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.site.v2ad <- jags(model.file = 'ectomod_v2aD.txt',
                        data = datalist, n.chains = 3, 
                        parameters.to.save = params,
                        inits = inits, n.burnin = 9000, 
                        n.iter = 20000, n.thin = 10)

samples.site.v2ad <- jags.samples(model.site.v2ad$model, 
                                  variable.names = c("WAIC"),
                                  n.iter = 5000, n.burnin = 1000,
                                  n.thin = 1, type = "mean")

# Writing WAIC table
WAIC.list <- list("GrassForb + HostAbundance + HostDiversity" = 
                    samples.site.v1ad[[1]])

# Compare
waic.frame.site <- data.frame(model = names(WAIC.list), 
                         WAIC = unlist(lapply(WAIC.list, sum)))

# load existing WAIC table
waic.table <- read.csv(file = "waic_site.csv")

# combine with other WAICs and write table
waic.frame.site <- rbind(waic.frame.site, waic.table)

waic.frame.site <- arrange(waic.frame.site, WAIC)

print(waic.frame.site)

# write.csv(waic.frame.site, row.names = F,
#           file = "waic_site.csv")

# Model selection: host ---------------------
# Full model
model <- jags(model.file = 'ectomod.txt', data = datalist,
              n.chains = 3, parameters.to.save = params,
              inits = inits, n.burnin = 9000, n.iter = 20000,
              n.thin = 10)

samples.host.full <- jags.samples(model$model, 
                                  variable.names = c("WAIC"),
                                  n.iter = 5000, n.burnin = 1000, 
                                  n.thin = 1, type = "mean")

# Intercept model
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, K = length(unique(hostspec)),
                 deadveg = comp.cov1, grassforb = comp.cov2,
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2','a3', 'a4', 'c1', 'c2', 'psi', 'theta', 
            'Z', 'Y')

model.host.inter <- jags(model.file = 'ectomod_HostIntercept.txt', 
              data = datalist, n.chains = 3, 
              parameters.to.save = params, inits = inits, 
              n.burnin = 9000, n.iter = 20000, n.thin = 10)

samples.host.inter <- jags.samples(model.host.inter$model, 
                                  variable.names = c("WAIC"),
                                  n.iter = 5000, n.burnin = 1000, 
                                  n.thin = 1, type = "mean")

# Host spec
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)),
                 deadveg = comp.cov1, grassforb = comp.cov2,
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b1', 'c1', 'c2',
            'psi', 'theta')

model.host.spec <- jags(model.file = 'ectomod_Spec.txt', 
                         data = datalist, n.chains = 3, 
                         parameters.to.save = params, inits = inits, 
                         n.burnin = 9000, n.iter = 20000, n.thin = 10)

samples.host.spec <- jags.samples(model.host.spec$model, 
                                   variable.names = c("WAIC"),
                                   n.iter = 5000, n.burnin = 1000, 
                                   n.thin = 1, type = "mean")

# Host mass
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b2', 'c1', 'c2',
            'psi', 'theta')

model.host.mass <- jags(model.file = 'ectomod_Mass.txt', 
                        data = datalist, n.chains = 3, 
                        parameters.to.save = params, inits = inits, 
                        n.burnin = 9000, n.iter = 20000, n.thin = 10)

samples.host.mass <- jags.samples(model.host.mass$model, 
                                  variable.names = c("WAIC"),
                                  n.iter = 5000, n.burnin = 1000, 
                                  n.thin = 1, type = "mean")

# Host sex
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b3', 'c1', 'c2',
            'psi', 'theta')

model.host.sex <- jags(model.file = 'ectomod_Sex.txt', 
                        data = datalist, n.chains = 3, 
                        parameters.to.save = params, inits = inits, 
                        n.burnin = 9000, n.iter = 20000, n.thin = 10)

samples.host.sex <- jags.samples(model.host.sex$model, 
                                  variable.names = c("WAIC"),
                                  n.iter = 5000, n.burnin = 1000, 
                                  n.thin = 1, type = "mean")

# Interaction term
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)),
                 deadveg = comp.cov1, grassforb = comp.cov2,
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.host.inter <- jags(model.file = 'ectomod_Interaction.txt', 
                       data = datalist, n.chains = 3, 
                       parameters.to.save = params, inits = inits, 
                       n.burnin = 9000, n.iter = 20000, n.thin = 10)

samples.host.inter <- jags.samples(model.host.inter$model, 
                                 variable.names = c("WAIC"),
                                 n.iter = 5000, n.burnin = 1000, 
                                 n.thin = 1, type = "mean")

# Spec + mass
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'c1', 'c2',
            'psi', 'theta')

model.host.specmass <- jags(model.file = 'ectomod_SpecMass.txt', 
                       data = datalist, n.chains = 3, 
                       parameters.to.save = params, inits = inits, 
                       n.burnin = 9000, n.iter = 20000, n.thin = 10)

samples.host.specmass <- jags.samples(model.host.specmass$model, 
                                 variable.names = c("WAIC"),
                                 n.iter = 5000, n.burnin = 1000, 
                                 n.thin = 1, type = "mean")

# Spec + sex
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)),
                 deadveg = comp.cov1, grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b1', 'b3', 'c1', 'c2',
            'psi', 'theta')

model.host.specsex <- jags(model.file = 'ectomod_SpecSex.txt', 
                            data = datalist, n.chains = 3, 
                            parameters.to.save = params, 
                            inits = inits, n.burnin = 9000, 
                            n.iter = 20000, n.thin = 10)

samples.host.specsex <- jags.samples(model.host.specsex$model, 
                                      variable.names = c("WAIC"),
                                      n.iter = 5000, n.burnin = 1000, 
                                      n.thin = 1, type = "mean")

# Mass + sex
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta')

model.host.sexmass <- jags(model.file = 'ectomod_SexMass.txt', 
                           data = datalist, n.chains = 3, 
                           parameters.to.save = params, 
                           inits = inits, n.burnin = 9000, 
                           n.iter = 20000, n.thin = 10)

samples.host.sexmass <- jags.samples(model.host.sexmass$model, 
                                     variable.names = c("WAIC"),
                                     n.iter = 5000, n.burnin = 1000, 
                                     n.thin = 1, type = "mean")

# Mass + Interaction
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b2', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.host.massinter <- jags(model.file = 'ectomod_MassInter.txt', 
                           data = datalist, n.chains = 3, 
                           parameters.to.save = params, 
                           inits = inits, n.burnin = 9000, 
                           n.iter = 20000, n.thin = 10)

samples.host.massinter <- jags.samples(model.host.massinter$model, 
                                     variable.names = c("WAIC"),
                                     n.iter = 5000, n.burnin = 1000, 
                                     n.thin = 1, type = "mean")

# Sex + Interaction
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)),
                 deadveg = comp.cov1, grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.host.sexinter <- jags(model.file = 'ectomod_SexInter.txt', 
                             data = datalist, n.chains = 3, 
                             parameters.to.save = params, 
                             inits = inits, n.burnin = 9000, 
                             n.iter = 20000, n.thin = 10)

samples.host.sexinter <- jags.samples(model.host.sexinter$model, 
                                       variable.names = c("WAIC"),
                                       n.iter = 5000, 
                                       n.burnin = 1000, n.thin = 1,
                                      type = "mean")

# Spec + Interaction
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)),
                 deadveg = comp.cov1, grassforb = comp.cov2,
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b1', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.host.specinter <- jags(model.file = 'ectomod_SpecInter.txt', 
                            data = datalist, n.chains = 3, 
                            parameters.to.save = params, 
                            inits = inits, n.burnin = 9000, 
                            n.iter = 20000, n.thin = 10)

samples.host.specinter <- jags.samples(model.host.specinter$model, 
                                      variable.names = c("WAIC"),
                                      n.iter = 5000, 
                                      n.burnin = 1000, n.thin = 1,
                                      type = "mean")

# Mass + Sex + Interaction
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.host.msi <- jags(model.file = 'ectomod_MSI.txt', 
                             data = datalist, n.chains = 3, 
                             parameters.to.save = params, 
                             inits = inits, n.burnin = 9000, 
                             n.iter = 20000, n.thin = 10)

samples.host.msi <- jags.samples(model.host.msi$model, 
                                       variable.names = c("WAIC"),
                                       n.iter = 5000, 
                                       n.burnin = 1000, n.thin = 1,
                                       type = "mean")

# Spec + Mass + Interaction
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.host.SpMI <- jags(model.file = 'ectomod_SpMI.txt', 
                       data = datalist, n.chains = 3, 
                       parameters.to.save = params, 
                       inits = inits, n.burnin = 9000, 
                       n.iter = 20000, n.thin = 10)

samples.host.SpMI <- jags.samples(model.host.SpMI$model, 
                                 variable.names = c("WAIC"),
                                 n.iter = 5000, 
                                 n.burnin = 1000, n.thin = 1,
                                 type = "mean")

# Spec + Sex + Interaction
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)),
                 deadveg = comp.cov1, grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b1', 'b3', 'c1', 'c2',
            'psi', 'theta', 'b4')

model.host.SpSI <- jags(model.file = 'ectomod_SpSI.txt', 
                        data = datalist, n.chains = 3, 
                        parameters.to.save = params, 
                        inits = inits, n.burnin = 9000, 
                        n.iter = 20000, n.thin = 10)

samples.host.SpSI <- jags.samples(model.host.SpSI$model, 
                                  variable.names = c("WAIC"),
                                  n.iter = 5000, 
                                  n.burnin = 1000, n.thin = 1,
                                  type = "mean")

# Spec + Mass + Sex
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates, abund = abund.vec,
                 div = div.vec)

params <- c('a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'b3', 'c1', 'c2',
            'psi', 'theta')

model.host.SpMS <- jags(model.file = 'ectomod_SpMS.txt', 
                        data = datalist, n.chains = 3, 
                        parameters.to.save = params, 
                        inits = inits, n.burnin = 9000, 
                        n.iter = 20000, n.thin = 10)

samples.host.SpMS <- jags.samples(model.host.SpMS$model, 
                                  variable.names = c("WAIC"),
                                  n.iter = 5000, 
                                  n.burnin = 1000, n.thin = 1,
                                  type = "mean")

# Set up WAIC table
WAIC.list <- list("Species + Mass + Sex" = 
                    samples.host.SpMS[[1]])

# Compare
waic.frame.host <- data.frame(model = names(WAIC.list), 
                              WAIC = unlist(lapply(WAIC.list, sum)))

# load existing WAIC table
waic.table <- read.csv(file = "waic_host.csv")

# combine with other WAICs and write table
waic.frame.host <- rbind(waic.frame.host, waic.table)

waic.frame.host <- arrange(waic.frame.host, WAIC)

print(waic.frame.host)

write.csv(waic.frame.host, row.names = F,
          file = "waic_host.csv")

# Covariate results ----------------------
# Veg composition: forest or not
a1 <- model$BUGSoutput$sims.list$a1

a1s <- data.frame(mean = apply(a1, 2, mean),
           lo = apply(a1, 2, quantile, 0.125),
           hi = apply(a1, 2, quantile, 0.875)) %>%
  mutate(sig = NA) %>%
  mutate(sig = case_when(lo > 0 | hi < 0 ~ "Yes",
                         TRUE ~ ""))
rownames(a1s) <- colnames(site.occ)[-c(1)]

veg1 <- ggplot(data = a1s, aes(x = rownames(a1s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_point(aes(x = rownames(a1s), y = 11, color = sig), shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  geom_hline(yintercept = 0)+
  labs(x = "Ectoparasite Species", y = "Leaf Litter Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        panel.grid = element_blank(), legend.position = "None",
        plot.margin = unit(c(5.5,5.5,5.5,60), "pt"))

# ggsave(filename = "vegforest_75_sansrares.jpeg", width = 10,
#        height = 4,units = "in")

# Veg comp: grass or forb
a2 <- model$BUGSoutput$sims.list$a2

a2s <- data.frame(mean = apply(a2, 2, mean),
           lo = apply(a2, 2, quantile, 0.125),
           hi = apply(a2, 2, quantile, 0.875)) %>%
      mutate(sig = case_when(lo > 0 | hi < 0 ~ "Yes",
                         TRUE ~ ""))
rownames(a2s) <- colnames(site.occ)[-c(1)]

veg2 <- ggplot(data = a2s, aes(x = rownames(a2s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = rownames(a1s), y = 4.5, color = sig), shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Grass/Forb Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None",
        plot.margin = unit(c(5.5,5.5,5.5,60), "pt"))

# ggsave(filename = "vegforb_75_sansrares.jpeg", width = 10,
#        height = 4, units = "in")

(veg1/veg2) + 
  plot_annotation(tag_levels = "a")

# ggsave(filename = "vegplts_mesosdisagg.jpeg", width = 10, height = 6,
#        units = "in", dpi = 600)

# Host abundance
a3 <- model$BUGSoutput$sims.list$a3

a3s <- data.frame(mean = apply(a3, 2, mean),
                  lo = apply(a3, 2, quantile, 0.125),
                  hi = apply(a3, 2, quantile, 0.875)) %>%
  mutate(sig = case_when(lo > 0 | hi < 0 ~ "Yes",
                         TRUE ~ ""))
rownames(a3s) <- colnames(site.occ)[-c(1)]

ggplot(data = a3s, aes(x = rownames(a3s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = rownames(a1s), y = 1.5, color = sig), shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Host Abundance")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None")

# ggsave(filename = "hostabund_75_interac.jpeg", width = 9, 
#        height = 4, units = "in")

# Host diversity
a4 <- model$BUGSoutput$sims.list$a4

a4s <- data.frame(mean = apply(a4, 2, mean),
                  lo = apply(a4, 2, quantile, 0.125),
                  hi = apply(a4, 2, quantile, 0.875)) %>%
  mutate(sig = case_when(lo > 0 | hi < 0 ~ "Yes",
                         TRUE ~ ""))
rownames(a4s) <- colnames(site.occ)[-c(1)]

ggplot(data = a4s, aes(x = rownames(a4s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = rownames(a1s), y = 10.5, color = sig), shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Host Diversity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None",
        plot.margin = unit(c(5.5,5.5,5.5,60), "pt"))

# ggsave(filename = "diversity_75_sansrares.jpeg", width = 9,
#        height = 4, units = "in")

# Effect of host spec
b1 <- model$BUGSoutput$sims.list$b1

b1lo <- as.data.frame(apply(b1, c(2,3), quantile, 0.125))
colnames(b1lo) <- levels(hostspec)
b1lo$Ecto <- sort(unique(mecto.caps$Ecto))

b1lo <- pivot_longer(b1lo, cols = -Ecto, 
                     names_to = "host", values_to = "lo")

b1hi <- as.data.frame(apply(b1, c(2,3), quantile, 0.875))
colnames(b1hi) <- levels(hostspec)
b1hi$Ecto <- sort(unique(mecto.caps$Ecto))

b1hi <- pivot_longer(b1hi, cols = -Ecto, 
                     names_to = "host", values_to = "hi")

b1mean <- as.data.frame(apply(b1, c(2,3), mean))
colnames(b1mean) <- levels(hostspec)
b1mean$Ecto <- sort(unique(mecto.caps$Ecto))

b1mean <- pivot_longer(b1mean, cols = -Ecto,
                       names_to = "host", values_to = "mean")

hilo <- inner_join(b1hi, b1lo, by = c("host", "Ecto"))
hilo$mean <- b1mean$mean

hilo.host <- hilo[which(hilo$lo > 0 | hilo$hi < 0),]

# Make interval plot for each ectoparasite species:
ecto.specs <- unique(hilo.host$Ecto)

plotlist <- list()
for(i in 1:length(ecto.specs)){
  spec <- hilo[hilo$Ecto == ecto.specs[i],]
  spec <- spec[spec$host != "PELE",]
  
  plotlist[[i]] <- ggplot(data = spec, aes(x = mean, y = host))+
    geom_point()+
    geom_errorbar(aes(xmin = lo, xmax = hi))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    labs(x = "Coefficient", y = "Host Species", title = ecto.specs[i])+
    scale_x_continuous(limits = c(-5, 5))+
    theme_bw(base_size = 12)+
    theme(panel.grid = element_blank())
}

# for(i in 1:length(plotlist)){
#   ggsave(plotlist[[i]], filename = paste("host", i, ".jpeg",
#                                          sep = ""),
#          width = 5, height = 3, units = "in")
# }

# host mass
b2 <- model$BUGSoutput$sims.list$b2

b2s <- data.frame(mean = colMeans(b2),
           lo = apply(b2, 2, quantile, 0.125),
           hi = apply(b2, 2, quantile, 0.875)) %>%
  mutate(sig = NA) %>%
  mutate(sig = case_when(lo > 0 | hi < 0 ~ "Yes",
                         TRUE ~ ""))
rownames(b2s) <- colnames(site.occ)[-c(1)]

ggplot(data = b2s, aes(x = rownames(b2s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = rownames(b2s), y = -0.8, color = sig), 
             shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Host Mass Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None",
        plot.margin = unit(c(5.5,5.5,5.5,60), "pt"))

# ggsave(filename = "hostmass_75_sansrares.jpeg", width = 9,
#        height = 4, units = "in")

# host sex
b3 <- model$BUGSoutput$sims.list$b3

b3s <- data.frame(mean = apply(b3, 2, mean),
           lo = apply(b3, 2, quantile, 0.125),
           hi = apply(b3, 2, quantile, 0.875)) %>%
  mutate(sig = NA) %>%
  mutate(sig = case_when(lo > 0 | hi < 0 ~ "Yes",
                         TRUE ~ ""))
rownames(b3s) <- colnames(site.occ)[-c(1)]

ggplot(data = b3s, aes(x = rownames(b3s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = rownames(b3s), y = 1.5, color = sig), 
             shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Host Sex Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None",
        plot.margin = unit(c(5.5,5.5,5.5,60), "pt"))

# ggsave(filename = "hostsex_75_sansrares.jpeg", width = 10,
#        height = 4, units = "in")

# Host/habitat interaction term
b4 <- model$BUGSoutput$sims.list$b4

b4lo <- as.data.frame(apply(b4, c(2,3), quantile, 0.125))
colnames(b4lo) <- levels(hostspec)
b4lo$Ecto <- sort(unique(mecto.caps$Ecto))

b4lo <- pivot_longer(b4lo, cols = -Ecto, 
                     names_to = "host", values_to = "lo")

b4hi <- as.data.frame(apply(b4, c(2,3), quantile, 0.875))
colnames(b4hi) <- levels(hostspec)
b4hi$Ecto <- sort(unique(mecto.caps$Ecto))

b4hi <- pivot_longer(b4hi, cols = -Ecto, 
                     names_to = "host", values_to = "hi")

b4mean <- as.data.frame(apply(b4, c(2,3), mean))
colnames(b4mean) <- levels(hostspec)
b4mean$Ecto <- sort(unique(mecto.caps$Ecto))

b4mean <- pivot_longer(b4mean, cols = -Ecto,
                       names_to = "host", values_to = "mean")

hilo <- inner_join(b4hi, b4lo, by = c("host", "Ecto"))
hilo$mean <- b4mean$mean

hilo.smol <- hilo[which(hilo$lo > 0 | hilo$hi < 0),]

# Make plots for sig. ectoparasite species:
ecto.specs <- unique(hilo.smol$Ecto)

# Get host-level occupancy:
theta <- model$BUGSoutput$sims.list$theta
theta.mean <- apply(theta, c(2,3), mean)

colnames(theta.mean) <- sort(unique(mecto.caps$Ecto))
theta.frame <- as.data.frame(theta.mean) %>%
  mutate(Tag = host.spec$Tag) %>%
  pivot_longer(cols = !Tag, names_to = "Ecto", values_to = "Occ")

# Set up data:
interac.frame <- mecto.caps %>%
  filter(Ecto %in% hilo.smol$Ecto) %>%
  select(Tag, Site, Ecto) %>%
  left_join(comp.frame[,c("Site", "pc1")], by = "Site") %>%
  left_join(host.spec, by = "Tag") %>%
  left_join(theta.frame, by = c("Tag", "Ecto"))

plotlist <- list()
for(i in 1:length(ecto.specs)){
  plotlist[[i]] <- ggplot(data = interac.frame[interac.frame$Ecto==ecto.specs[i],], 
                     aes(x = pc1, y = Occ, color = Abbrev, 
                         fill = Abbrev))+
    geom_point()+
    geom_smooth(method = 'lm', se = T)+
    scale_color_viridis_d(end = 0.95, name = "Host")+
    scale_fill_viridis_d(end = 0.95, name = "Host")+
    labs(x = "PC1", y = "Occupancy Probability (Host)",
         title = ecto.specs[i])+
    theme_bw(base_size = 12)+
    theme(panel.grid = element_blank())
}

(plotlist[[1]]|plotlist[[2]])/
  (plotlist[[3]]|plotlist[[4]])+
  plot_layout(guides = "collect")&
  plot_annotation(tag_levels = "a")

# ggsave(filename = "interaction_plts_sansrares.jpeg", width = 8,
#        height = 6, units = "in", dpi = 600)

# Effect of capture no
c1 <- model$BUGSoutput$sims.list$c1

quantile(c1, c(0.125, 0.875))

c1s <- data.frame(mean = apply(c1, 2, mean), 
           lo = apply(c1, 2, quantile, 0.125), 
           hi = apply(c1, 2, quantile, 0.875)) %>%
  mutate(sig = NA) %>%
  mutate(sig = case_when(lo > 0 | hi < 0 ~ "Yes",
                         TRUE ~ ""))
rownames(c1s) <- colnames(site.occ)[-c(1)]

c1_plt <- ggplot(data = c1s, aes(x = rownames(c1s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = rownames(c1s), y = 0.6, color = sig), shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Capture no. Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        panel.grid = element_blank(), legend.position = "None",
        plot.margin = unit(c(5.5,5.5,5.5,60), "pt"))

# ggsave(filename = "capno_75_sansrares.jpeg", width = 10,
#        height = 4, units = "in")

# Julian date
c2 <- model$BUGSoutput$sims.list$c2
quantile(c2, c(0.125, 0.875))

c2s <- data.frame(mean = apply(c2, 2, mean), 
           lo = apply(c2, 2, quantile, 0.125), 
           hi = apply(c2, 2, quantile, 0.875)) %>%
  mutate(sig = case_when(lo > 0 | hi < 0 ~ "Yes",
                         TRUE ~ ""))
rownames(c2s) <- colnames(site.occ)[-c(1)]

c2_plot <- ggplot(data = c2s, aes(x = rownames(c2s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = rownames(c2s), y = 1.2, color = sig), shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Julian Date Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None",
        plot.margin = unit(c(5.5,5.5,5.5,60), "pt"))

# ggsave(filename = "datecov_75_sansrares.jpeg", width = 10,
#        height = 4, units = "in")

(c1_plt / c2_plot)+
  plot_annotation(tag_levels = "a")

ggsave(filename = "detcovs_iscapagg.jpeg", width = 10, height = 8, 
       units = "in", dpi = 600)

# Bayesian r-squared ---------------------
# Calculation for model variance
var.yhat <- function(y,n=necto){
  yvar <- matrix(NA, ncol = n)
  for(i in 1:necto){
    yvar[,i] <- var(y[,i])
  }
  return(yvar)
}

var.res <- function(x,n=necto){
  resvar <- matrix(NA, ncol = n)
  for(i in 1:necto){
    resvar[,i] <- (1/necto)*sum(x[,i]*(1-x[,i]))
  }
  return(resvar)
}

# Calculation for Bayesian r-squared
bayes.r2 <- function(x,y){
  r2 <- y/(y+x)
  return(r2)
}

# Site-level
psi <- model$BUGSoutput$sims.list$psi

site.yhat <- apply(psi, 1, var.yhat)
site.res <- apply(psi, 1, var.res)
site.r2 <- Apply(data = list(site.res, site.yhat), margins = 2, 
                 fun = bayes.r2)[[1]]

# Host-level
theta <- model$BUGSoutput$sims.list$theta

host.yhat <- apply(theta, 1, var.yhat)
host.res <- apply(theta, 1, var.res)
host.r2 <- Apply(data = list(host.res, host.yhat), margins = 2, 
                 fun = bayes.r2)[[1]]

# R2 Figures ------------------------
# Vector of species names
ecto.specs <- names(ecto.list2)

# Get estimated # hosts
# occ.out <- model$BUGSoutput$sims.list$Y
# occ.mat <- apply(occ.out, c(2,3), mean)
# colnames(occ.mat) <- ecto.specs
# colnames(occ.mat)[15] <- "UnknownMesostigmata"
# 
# host.frame <- as.data.frame(cbind(occ.mat, hostspec)) %>%
#   mutate(hostID = 1:nrow(.)) %>%
#   pivot_longer(Ctenophthalmus_pseudagyrtes:UnknownMesostigmata, 
#                values_to = "OccProb", names_to = "Ecto") %>%
#   filter(OccProb > 0.75) %>%
#   group_by(Ecto) %>%
#   summarise(nhost = length(unique(hostspec)))
# Infested species did not change- keep old code

# Site-level r-squared
mean.site.r2 <- apply(site.r2, 1, mean, na.rm = T)
median.site.r2 <- apply(site.r2, 1, median)
site.df <- data.frame(Ecto = ecto.specs, site.r2 = mean.site.r2,
                      median.r2 = median.site.r2)

# Get family and order names for each species
site.df <- site.df %>%
  left_join(n.host, by = c("Ecto" = "ecto")) %>%
  mutate(Order = case_when(startsWith(Ecto, "Ixodes") ~ "Ixodida",
                           TRUE ~ Order)) %>%
  mutate(hosts = as.numeric(hosts)) %>%
  mutate(hosts = case_when(Ecto == "Ixodes_scapularis_larva" ~ 5,
                           Ecto == "Ixodes_scapularis_nymph" ~ 5,
                           TRUE ~ hosts))

# If counting host by subfamily:
# site.df <- site.df %>%
#   left_join(n.host.sub, by = c("Ecto" = "ecto"))

# If counting host by family:
# site.df <- site.df %>%
#   left_join(n.host.fam, by = c("Ecto" = "ecto"))

# full fig
site.r2.fig <- ggplot(data = site.df, aes(x = Order, 
                                        y = site.r2))+
  geom_boxplot(fill = "lightgray", outlier.shape = NA)+
  # geom_jitter(aes(color = Order))+
  scale_color_viridis_d(end = 0.95)+
  labs(y = bquote("Site"~R^2))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), legend.position = "none")

# ggsave(filename = "siter2_order_full.jpeg", height = 4, width = 6,
#        units = "in", dpi = 600)

# Host-level r-squared
mean.host.r2 <- apply(host.r2, 1, mean, na.rm = T)
host.df <- data.frame(Ecto = ecto.specs, host.r2 = mean.host.r2)

# Host species:
host.df <- host.df %>%
  left_join(n.host, by = c("Ecto" = "ecto"))

# Host subfamilies:
# host.df <- host.df %>%
#   left_join(n.host.sub, by = c("Ecto" = "ecto"))

# Host families:
# host.df <- host.df %>%
#   left_join(n.host.fam, by = c("Ecto" = "ecto"))

# Merge with fleas
host.df <- host.df %>%
  mutate(Order = case_when(startsWith(Ecto, "Ixodes") ~ "Ixodida",
                           TRUE ~ Order)) %>%
  mutate(hosts = as.numeric(hosts)) %>%
  mutate(hosts = case_when(Ecto == "Ixodes_scapularis_larva" ~ 5,
                           Ecto == "Ixodes_scapularis_nymph" ~ 5,
                           TRUE ~ hosts))

# Final fig
host.r2.fig <- ggplot(data = host.df, aes(x = Order, y = host.r2))+
  geom_boxplot(fill = 'lightgray', outlier.shape = NA)+
  # geom_jitter(aes(color = Order))+
  scale_color_viridis_d(end = 0.95)+
  labs(y = bquote("Host"~R^2))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "hostr2_order_full.jpeg", height = 4, width = 6,
       # units = 'in', dpi = 600)

r2fig <- site.r2.fig + host.r2.fig +
  plot_annotation(tag_levels = "a")

# ggsave(filename = "r2fig_rare.jpeg", height = 3, width = 10,
#        units = 'in',dpi = 600)

# R-squared vs. host specificity -------------------
summary(lm(data = site.df, formula = hosts~site.r2))
summary(lm(data = host.df, formula = hosts~host.r2))
# This is wild

# look at site r2 sans o. leucopus:
summary(lm(data = site.df[site.df$Ecto != "Orchopeas_leucopus",],
           formula = hosts~site.r2))
qplot(data = site.df[site.df$Ecto != "Orchopeas_leucopus",],
      x = site.r2, y = hosts)
# no pattern when o. leucopus is dropped

r2site <- ggplot(data = site.df, aes(x = site.r2, y = hosts))+
  geom_point(aes(color = Order), size = 2)+
  # geom_smooth(se = F, method = 'lm', color = "black")+
  scale_color_viridis_d(end = 0.95)+
  labs(x = bquote("Site"~R^2), y = "Number of Hosts")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

r2host <- ggplot(data = host.df, aes(x = host.r2, y = hosts))+
  geom_point(aes(color = Order), size = 2)+
  geom_smooth(se = F, method = 'lm', color = "black")+
  scale_color_viridis_d(end = 0.95)+
  labs(x = bquote("Host"~R^2), y = "Number of Hosts")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

r2site + r2host +
  plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "a")

# ggsave(filename = "r2nhost_sansrares.jpeg", width = 10,
#        height = 3, units = "in", dpi = 600)

# Ecto abundance vs specificity -------------
abund.hosts <- mecto.clean %>%
  select(SampleNo., Genus, Species, Other) %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
  #Comment this line to keep Iscap aggregated:
  mutate(Ecto = case_when(Ecto == "Ixodes_scapularis" ~
                            paste(Ecto, Other, sep = "_"),
                          TRUE ~ Ecto)) %>%
  filter(Ecto != "") %>%
  group_by(Ecto) %>%
  summarise(Abund = n()) %>%
  left_join(n.host, by = c("Ecto" = "ecto"))

ggplot(data = abund.hosts, aes(x = hosts, y = log(Abund)))+
  geom_point(aes(color = Order), size = 2)+
  geom_smooth(method = 'lm', color = "black")+
  scale_color_viridis_d(end = 0.9)+
  labs(x = "Hosts", y = "log(Abundance)")+
  theme_bw()+
  theme(panel.grid = element_blank())

summary(lm(data = abund.hosts, log(Abund)~hosts))

# ggsave(filename="specificityabund.jpeg", width = 5, height = 3,
#        units = "in")

# P. leucopus captured ectos --------
# Get parasite counts
pele.necto <- mecto.clean %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
  #Comment this line to keep Iscap aggregated:
  # mutate(Ecto = case_when(Ecto == "Ixodes_scapularis" ~
  #                           paste(Ecto, Other, sep = "_"),
  #                         TRUE ~ Ecto)) %>%
  select(Site, Day, Tag, HostName, Ecto) %>%
  filter(HostName == "P . leucopus") %>%
  left_join(comp.frame[,c("pc1", "Site","Habitat")], by = "Site") %>%
  mutate(Ecto = case_when(Ecto == "" ~ NA,
                          TRUE ~ Ecto)) %>%
  group_by(Tag, pc1) %>%
  summarise(necto = length(na.omit(Ecto))) 

# Check for overdispersion: mean & var should be close
mean(pele.necto$necto)
var(pele.necto$necto)
# pretty overdispersed: zero-infl neg binomial it is

pele.mod <- zeroinfl(necto ~ pc1 | pc1, data = pele.necto, 
                 dist = 'negbin')
summary(pele.mod)

# Dispersion Statistic
E2 <- resid(pele.mod, type = "pearson")
N  <- nrow(pele.necto)
p  <- length(coef(pele.mod)) + 1 # '+1' is due to theta
sum(E2^2) / (N - p)
# close to 1- looks good

# Plot it
pele.plt <- ggplot(data = pele.necto, aes(x = pc1, y = necto))+
  geom_jitter()+
  geom_smooth(aes(y = predict(pele.mod, type = "count")), 
              color = "black")+
  labs(x = "Vegetation Coefficient (PC1)", 
       y = "Ectoparasites Collected")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "pele_infest_iscapagg.jpeg", width = 5, height = 3,
#        units = "in", dpi = 600)
   
# T. striatus captured ectos --------------
# Get parasite counts
tast.necto <- mecto.clean %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
  #Comment this line to keep Iscap aggregated:
  # mutate(Ecto = case_when(Ecto == "Ixodes_scapularis" ~
  #                           paste(Ecto, Other, sep = "_"),
  #                         TRUE ~ Ecto)) %>%
  dplyr::select(Site, Day, Tag, HostName, Ecto) %>%
  filter(HostName == "T . striatus")  %>%
  left_join(comp.frame[,c("pc1", "Site")], by = "Site") %>%
  mutate(Ecto = case_when(Ecto == "" ~ NA,
                          TRUE ~ Ecto)) %>%
  group_by(Tag, pc1) %>%
  summarise(necto = length(na.omit(Ecto))) 

# Check for overdispersion: mean & var should be close
mean(tast.necto$necto)
var(tast.necto$necto)
# pretty overdispersed: zero-infl neg binomial it is

tast.mod <- zeroinfl(necto ~ pc1 | pc1, data = tast.necto, 
                     dist = 'negbin')
summary(tast.mod)

# Dispersion Statistic
E2 <- resid(tast.mod, type = "pearson")
N  <- nrow(tast.necto)
p  <- length(coef(tast.mod)) + 1 # '+1' is due to theta
sum(E2^2) / (N - p)
# not as close to 1 but still ok

# Plot it
tast.plt <- ggplot(data = tast.necto, aes(x = pc1, y = necto))+
  geom_jitter()+
  geom_smooth(aes(y = predict(tast.mod, type = "count")), 
              color = "black")+
  # add geom_smooth if I can get it to work
  labs(x = "Vegetation Coefficient (PC1)", 
       y = "Ectoparasites Collected")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "tast_infest.jpeg", width = 5, height = 3,
#        units = "in", dpi = 600)

ecto.habs <- mecto.clean %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
  select(Site, Day, Tag, HostName, Ecto) %>%
  filter(HostName == "T . striatus") %>%
  left_join(comp.frame[,c("pc1", "Site","Habitat")], by = "Site") %>%
  select(Ecto,Habitat) %>%
  filter(Ecto != "") %>%
  group_by(Ecto, Habitat) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = Habitat, values_from = count, 
              values_fill = 0)

props <- apply(ecto.habs[,2:4], 2, function(x) x/sum(x))
rownames(props) <- ecto.habs$Ecto
# Not enough parasites in fields/farms for this to be informative

# both species
(pele.plt|tast.plt) +
  plot_annotation(tag_levels = "a")

# ggsave(filename = "infest_both.jpeg", width = 7, height = 3,
#        units = "in", dpi = 600)

# B brevicauda captured ectos ------------
# probably could write a function for this
# Get parasite counts
blbr.necto <- mecto.clean %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
  #Comment this line to keep Iscap aggregated:
  mutate(Ecto = case_when(Ecto == "Ixodes_scapularis" ~
                            paste(Ecto, Other, sep = "_"),
                          TRUE ~ Ecto)) %>%
  dplyr::select(Site, Day, Tag, HostName, Ecto) %>%
  filter(HostName == "B . brevicauda")  %>%
  left_join(comp.frame[,c("pc1", "Site")], by = "Site") %>%
  mutate(Ecto = case_when(Ecto == "" ~ NA,
                          TRUE ~ Ecto)) %>%
  group_by(Tag, pc1) %>%
  summarise(necto = length(na.omit(Ecto))) 

# Only 2 parasites from BLBR; this will not work

# Check for spatial autocorrelation/plot sites---------------
traplines <- st_read("UpdatedTracks.kml")

traplines <- st_crop(traplines, xmin = -73.4, xmax = -72.75, 
                   ymin = 44.1, ymax = 44.6)

# Convert to midpoints
st_line_midpoints <- function(sf_lines = NULL){
  
  g <- st_geometry(sf_lines)
  
  g_mids <- lapply(g, function(x){
    
    coords <- as.matrix(x)
    
    # this is just a copypaste of View(maptools:::getMidpoints):
    get_mids <- function (coords){
      dist <- sqrt((diff(coords[, 1])^2 + (diff(coords[, 2]))^2))
      dist_mid <- sum(dist)/2
      dist_cum <- c(0, cumsum(dist))
      end_index <- which(dist_cum > dist_mid)[1]
      start_index <- end_index - 1
      start <- coords[start_index, ]
      end <- coords[end_index, ]
      dist_remaining <- dist_mid - dist_cum[start_index]
      mid <- start + (end - start) * (dist_remaining/dist[start_index])
      return(mid)
    }
    
    mids <- st_point(get_mids(coords))
  })
}
trap.list <- st_line_midpoints(traplines)

# Build the data frame
trap.df <- as.data.frame(do.call(rbind, trap.list))
names(trap.df) <- c("X", "Y")
trap.df$Name <- traplines$Name

# Add habitat
trap.df <- trap.df %>%
  left_join(y=select(veg.site, Site, Habitat), 
            by = c("Name" = "Site"))

# Plot it
alldat <- map_data("county")

chittenden <- alldat[alldat$region=="vermont" &
                       alldat$subregion=="chittenden",]

ggplot(data = chittenden, aes(x = long, y = lat, fill = Habitat))+
  geom_polygon(color = "black", fill = "lightgray")+
  geom_point(data = trap.df, aes(x = X, y = Y), 
              size = 3, alpha = 0.6,
              pch = 21, color = "black")+
  scale_fill_viridis_d()+
  # North arrow isn't working
  # north()+
  theme_classic(base_size = 14)+
  theme(axis.title = element_blank(), axis.text = element_blank())

# ggsave(file = "sitemap.jpeg", dpi = 600)

# Use MMRR to check for spatial autocorrelation
# If autocorrelated, geographic distance will be sig predictor

vegdist <- as.matrix(dist(comp.cov1,upper = T))

trap.xy <- trap.df[-11, 1:2]
eucdist <- as.matrix(dist(trap.xy, upper = T))

mamm.counts <- mamm.2020 %>%
  distinct_at(.vars = vars(Site, Abbrev, Tag), .keep_all = T) %>%
  group_by(HostName, Site) %>%
  summarise(Count = n()) %>%
  pivot_wider(names_from=HostName, values_from=Count,
              values_fill = 0)

host.dist <- as.matrix(dist(mamm.counts[,-1]))

ecto.dist <- as.matrix(dist(sitemax, upper = T))

MMRR<-function(Y,X,nperm=999){
  #compute regression coefficients and test statistics
  nrowsY<-nrow(Y)
  y<-unfold(Y)
  if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
  Xmats<-sapply(X,unfold)
  fit<-lm(y~Xmats)
  coeffs<-fit$coefficients
  summ<-summary(fit)
  r.squared<-summ$r.squared
  tstat<-summ$coefficients[,"t value"]
  Fstat<-summ$fstatistic[1]
  tprob<-rep(1,length(tstat))
  Fprob<-1
  
  #perform permutations
  for(i in 1:nperm){
    rand<-sample(1:nrowsY)
    Yperm<-Y[rand,rand]
    yperm<-unfold(Yperm)
    fit<-lm(yperm~Xmats)
    summ<-summary(fit)
    Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
    tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
  }
  
  #return values
  tp<-tprob/(nperm+1)
  Fp<-Fprob/(nperm+1)
  names(r.squared)<-"r.squared"
  names(coeffs)<-c("Intercept",names(X))
  names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
  names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
  names(Fstat)<-"F-statistic"
  names(Fp)<-"F p-value"
  return(list(r.squared=r.squared,
              coefficients=coeffs,
              tstatistic=tstat,
              tpvalue=tp,
              Fstatistic=Fstat,
              Fpvalue=Fp))
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR

unfold<-function(X){
  x<-vector()
  for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
  return(x)
}

# Check for independence between x values
MMRR(vegdist, list(eucdist)) 
MMRR(host.dist, list(eucdist))
MMRR(host.dist, list(vegdist))
# We good

# Test for autocorrelation
MMRR(ecto.dist, list(vegdist, eucdist, host.dist))
# No spatial autocorrelation

# Counts by habitat ------------------------
parasite.frame <- mecto.clean %>%
  filter(is.na(Order) == F) %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
  #Comment this line to keep Iscap aggregated:
  mutate(Ecto = case_when(Ecto == "Ixodes_scapularis" ~
                            paste(Ecto, Other, sep = "_"),
                          TRUE ~ Ecto)) %>%
  select(Site, Abbrev, Ecto) %>%
  left_join(select(veg.site, Site, Habitat), by = "Site")

host.props <- mecto.clean %>%
  select(Site, Tag, Abbrev) %>%
  distinct() %>%
  left_join(select(veg.site, Site, Habitat), by = "Site") %>%
  group_by(Habitat, Abbrev) %>%
  summarise(count = n()) %>%
  mutate(prop = count/sum(count)) %>%
  ungroup()

host.all <- expand(host.props,Habitat,Abbrev)
host.all <- full_join(host.props, host.all) 

ggplot(host.all, aes(x = Abbrev, y = count, fill = Habitat))+
  geom_col(position = "dodge")+
  labs(x = "Host Species", y = "Captured Individuals")+
  scale_fill_viridis_d(end = 0.9)+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "hostcounts.jpeg", width = 10,
#        height = 3, units = "in", dpi = 600)

# O. leucopus
orle.abund <- parasite.frame %>%
  filter(Ecto == "Orchopeas_leucopus") %>%
  group_by(Habitat, Abbrev) %>%
  summarise(count = n()) %>%
  mutate(ecto.prop = count/sum(count)) %>%
  ungroup() %>%
  left_join(dplyr::select(host.props,-count), 
            by = c("Abbrev", "Habitat")) %>%
  mutate(adj_count = ecto.prop/prop) %>%
  ungroup() %>%
  mutate(Abbrev = factor(Abbrev, levels = c("PELE", "PEMA", "MYGA",
                                            "ZAHU", "MIPE"))) %>%
  mutate(Habitat = factor(Habitat, levels = c("Farm","Field",
                                              "Forest")))

orle.all <- expand(orle.abund,Habitat,Abbrev)
orle.all <- full_join(orle.abund, orle.all)  

orle.fig <- ggplot(data = orle.all, aes(x = Abbrev, y = adj_count, 
                            fill = Habitat))+
  geom_col(position = "dodge")+
  labs(x = "Host Species", y = "Adjusted Count")+
  scale_fill_viridis_d(end = 0.9)+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# Larval I. scap
iscl.abund <- parasite.frame %>%
  filter(Ecto == "Ixodes_scapularis_larva") %>%
  group_by(Habitat, Abbrev) %>%
  summarise(count = n()) %>%
  mutate(ecto.prop = count/sum(count)) %>%
  ungroup() %>%
  left_join(dplyr::select(host.props,-count), 
            by = c("Abbrev", "Habitat")) %>%
  mutate(adj_count = ecto.prop/prop) %>%
  ungroup() %>%
  mutate(Abbrev = factor(Abbrev, levels = c("PELE", "MIPE", "PEMA",
                                            "ZAHU"))) %>%
  mutate(Habitat = factor(Habitat, levels = c("Farm","Field", 
                                              "Forest")))

iscl.all <- expand(iscl.abund,Habitat,Abbrev)
iscl.all <- full_join(iscl.abund, iscl.all)  

iscl.fig <- ggplot(data = iscl.all, aes(x = Abbrev, y = adj_count,
                                        fill = Habitat))+
  geom_col(position = "dodge")+
  labs(x = "Host Species", y = "Adjusted Count")+
  scale_fill_viridis_d(end = 0.9)+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# Nymph iscap
iscn.abund <- parasite.frame %>%
  filter(Ecto == "Ixodes_scapularis_nymph") %>%
  group_by(Habitat, Abbrev) %>%
  summarise(count = n()) %>%  
  mutate(ecto.prop = count/sum(count)) %>%
  ungroup() %>%
  left_join(dplyr::select(host.props,-count), 
            by = c("Abbrev", "Habitat")) %>%
  mutate(adj_count = ecto.prop/prop) %>%
  ungroup() %>%
  mutate(Abbrev = factor(Abbrev, levels = c("MIPE", "TAST", "ZAHU",
                                            "PELE", "PEMA"))) %>%
  mutate(Habitat = factor(Habitat, levels = c("Farm", "Field", 
                                              "Forest"))) %>%
  ungroup()

iscn.all <- expand(iscn.abund,Habitat,Abbrev)
iscn.all <- full_join(iscn.abund, iscn.all)  

iscn.fig <- ggplot(data = iscn.all, 
                   aes(x = Abbrev, y = adj_count,
                       fill = Habitat))+
  geom_col(position = "dodge")+
  labs(x = "Host Species", y = "Adjusted Count")+
  scale_fill_viridis_d(end = 0.9)+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# All iscap
isc.abund <- parasite.frame %>%
  filter(Ecto == "Ixodes_scapularis") %>%
  group_by(Habitat, Abbrev) %>%
  summarise(count = n()) %>%  
  mutate(ecto.prop = count/sum(count)) %>%
  ungroup() %>%
  left_join(dplyr::select(host.props,-count), 
            by = c("Abbrev", "Habitat")) %>%
  mutate(adj_count = ecto.prop/prop) %>%
  ungroup() %>%
  mutate(Abbrev = factor(Abbrev, levels = c("MIPE", "TAST", "ZAHU",
                                            "PELE", "PEMA"))) %>%
  mutate(Habitat = factor(Habitat, levels = c("Farm", "Field", 
                                              "Forest"))) %>%
  ungroup()

isc.all <- expand(isc.abund,Habitat,Abbrev)
isc.all <- full_join(isc.abund, iscn.all)  

isc.fig <- ggplot(data = isc.all, 
                   aes(x = Abbrev, y = adj_count,
                       fill = Habitat))+
  geom_col(position = "dodge")+
  labs(x = "Host Species", y = "Adjusted Count")+
  scale_fill_viridis_d(end = 0.9)+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# Megabothris quirini
mequ.abund <- parasite.frame %>%
  filter(Ecto == "Megabothris_quirini") %>%
  group_by(Habitat, Abbrev) %>%
  summarise(count = n()) %>%
  mutate(ecto.prop = count/sum(count)) %>%
  ungroup() %>%
  left_join(dplyr::select(host.props,-count), 
              by = c("Abbrev", "Habitat")) %>%
  mutate(adj_count = ecto.prop/prop) %>%
  ungroup() %>%
  mutate(Abbrev = factor(Abbrev, levels = c("MIPE", "PELE", "ZAHU",
                                              #"MYGA", "NAIN",
                                            "PEMA"))) %>%
  mutate(Habitat = factor(Habitat, levels = c("Farm", "Field", 
                                                "Forest"))) %>%
  ungroup()

mequ.all <- expand(mequ.abund,Habitat,Abbrev)
mequ.all <- full_join(mequ.abund, mequ.all)  

mequ.fig <- ggplot(data = mequ.all, aes(x = Abbrev, y = adj_count, 
                                        fill = Habitat))+
  geom_col(position = "dodge")+
  labs(x = "Host Species", y = "Adjusted Count")+
  scale_fill_viridis_d(end = 0.9)+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(plot = mequ.fig, filename = "mequ_counts_adj.jpeg", 
#        width = 4, height = 3, units = "in", dpi = 600)

(iscl.fig + iscn.fig)/
  (orle.fig + mequ.fig) +
  plot_annotation(tag_levels = "a")+
  plot_layout(guides = "collect")

ggsave(filename = "generalist_counts_adj.jpeg", width = 10,
       height = 6, units = "in", dpi = 600)

# ggsave(plot = isc.fig, filename = "isc_counts_adj.jpeg", width = 4,
#        height = 3, units = "in", dpi = 600)
