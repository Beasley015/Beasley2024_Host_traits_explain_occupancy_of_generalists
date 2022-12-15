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
library(R2jags)
library(multiApply)
library(agricolae)
library(gridExtra)
library(grid)
library(rstatix)
library(patchwork)

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
  select(Site:Day, Tag, Abbrev:Mass, Ecto, Desc.) %>%
  filter(Tag != "") %>%
  filter(Desc. == "") %>%
  filter(!Abbrev %in% c("", "PE??", "SOCI")) %>%
  rename("SampleNo." = Ecto) %>%
  select(!Desc.)

# Ectos 
ecto.clean <- ecto.raw %>%
  filter(Other != "Springtail?") %>%
  .[!(.$Order == "Ixodida" & .$Species == ""),]

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
  mutate(Genus = case_when(Order == "Psocodea" ~ "Unknown Psocodea",
                           Order == "Mesostigmata" & Genus == "" ~
                             "Unknown Mesostigmata",
                           TRUE ~ Genus))

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
                                  "Mesostigmata", "Diptera",
                                  "Psocodea"))+
  labs(x = "Host Species", y = "Parasite Species")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("ectosraster.jpeg", height = 4, width = 8.5,
#        units = "in")

# Hosts per ecto species
n.host <- mamm.ecto %>%
  filter(SampleNo. != "") %>%
  mutate(Species = case_when(Genus == "Ixodes" & 
                               Species == "scapularis" ~
                               paste(Species, Other, sep = "_"),
                             TRUE ~ Species)) %>%
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
  geom_point()+
  geom_smooth(method = 'lm', se = F, color = "black")+
  labs(x = "Host Richness", y = "Ectoparasite Richness")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# Prepping data for MSOM ----------------------
# Clean data
mecto.clean <- mamm.ecto %>%
  filter(!(SampleNo. != "" & is.na(Species) == T)) %>%
  filter(Tag != "RECAP" & Tag != "DESC")

# Select columns and add occupancy column
mecto.smol <- mecto.clean %>%
  select(Tag, Site, Day, Genus, Species, Other) %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
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

mecto.caps <- mecto.caps %>%
  group_by(Tag) %>%
  filter(is.na(`2`) == F | sum(`1`) >= 1) %>%
  ungroup()

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
                 dim=c(nrow(ecto.list2[[1]]), 8, 
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
  theme_bw(base_size = 14)+
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

# ggsave(filename = "pc1.jpeg", width = 4, height = 4, units = "in",)

# Site covariate: Veg complexity -----------------
height.pc <- prcomp(veg.site[,4:9])

height.cov <- height.pc$x[,1]

height.cov <- as.vector(scale(height.cov))

height.frame <- cbind(veg.site, pc1 = height.cov, 
                      pc2 = height.pc$x[,2])

p <- ggplot(data = height.frame, aes(x = pc1, y = pc2, 
                                     color = Habitat))+
  geom_point(size = 2)+
  scale_color_viridis_d(end = 0.95)+
  labs(x = "PC1 (Decreasing Complexity)", y = "PC2")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())+
  annotation_custom(grob = textGrob(label = "10 cm tall",
                                    rot = 90),
                    ymin = -0.5, ymax = -0.5, xmin = -2.35,
                    xmax = -2.35)+
  annotation_custom(grob = textGrob(label = "30 cm tall", rot = 90),
                    ymin = 0.4, ymax = 0.4, xmin = -2.35,
                    xmax = -2.35)

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

# ggsave(filename = "pc2.jpeg", width = 4, height = 4, 
#        units = "in")

# Both pc's ---------------------
pcs <- data.frame(Site = height.frame$Site, 
                  Habitat = height.frame$Habitat,
                  CompPC = comp.frame$pc1, 
                  HeightPC = height.frame$pc1)

ggplot(data = pcs, aes(x = CompPC, y = HeightPC, color = Habitat))+
  geom_point(size = 2)+
  labs(x = "Composition PC", y = "Vertical Structure PC")+
  scale_color_viridis_d(end = 0.95)+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

cor.test(x = pcs$HeightPC, y = pcs$CompPC)
# Veg cover and complexity highly correlated- use cover only

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
      
      b0[i] ~ dnorm(b0.mean, tau.b0)
      b1[i,1] <- 0
      for(s in 2:K){
        b1[i,s] ~ dnorm(0, 1)
      }
      b2[i] ~ dnorm(b2.mean, tau.b2)
      b3[i] ~ dnorm(b3.mean, tau.b3)
      
      c0[i] ~ dnorm(c0.mean, tau.c0)
      c1[i] ~ dnorm(c1.mean, tau.c1)
      c2[i] ~ dnorm(c2.mean, tau.c2)
    
      # Occupancy model: Site
      for(j in 1:nsite){
      
        logit(psi[j,i]) <- a0[i] + a1[i]*deadveg[j] + 
                            a2[i]*grassforb[j] 
        Z[j,i] ~ dbern(psi[j,i])

        # Occupancy model: Host
        for(k in tagmat[1:hostvec[j],j]){
          # May need to create matrix of host sex
          sex[k,i] ~ dbern(probs)
        
          logit(theta[k,i]) <- b0[i] + b1[i,host[k]] + b2[i]*mass[k]
                                      + b3[i]*sex[k,i]
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
                 julian = dates)

params <- c('a1', 'a2', 'b1', 'b2', 'b3', 'c1', 'c2', 'psi',
            'theta')

# Init values
inits <- function(){
  inits <- list(
    Z = sitemax,
    Y = hostmax
  )
}

# Send model to JAGS
# model <- jags(model.file = 'ectomod.txt', data = datalist,
#               n.chains = 3, parameters.to.save = params,
#               inits = inits, n.burnin = 5000, n.iter = 12000,
#               n.thin = 12)

# Save model
# saveRDS(model, file = "ectomod.rds")
# saveRDS(model, file = "sans1cap.rds")
# saveRDS(model, file = "sans0caps0ecto.rds")

# Load model
model <- readRDS("ectomod.rds")
# model <- readRDS("sans1cap.rds")
# model <- readRDS("sans0caps0ecto.rds")

# Model selection: site ------------------------
# Full model
model <- jags(model.file = 'ectomod.txt', data = datalist,
              n.chains = 3, parameters.to.save = params,
              inits = inits, n.burnin = 4000, n.iter = 12000,
              n.thin = 15)

samples.site.full <- jags.samples(model$model, 
                             variable.names = c("WAIC"),
                             n.iter = 5000, n.burnin = 1000, 
                             n.thin = 1, type = "mean")

# Intercept model
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, 
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates)

params <- c('b1', 'b2', 'b3', 'c1', 'c2', 'psi','theta', 'Z', 'Y')

model.site.intercept <- jags(model.file = 'ectomod_SiteIntercept.txt',
                             data = datalist, n.chains = 3, 
                             parameters.to.save = params,
                             inits = inits, n.burnin = 4000, 
                             n.iter = 12000, n.thin = 15)

samples.site.intercept <- jags.samples(model.site.intercept$model, 
                                       variable.names = c("WAIC"),
                                       n.iter = 5000, n.burnin = 1000,
                                       n.thin = 1, type = "mean")

# Cov 1 only
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates)

params <- c('a1', 'b1', 'b2', 'b3', 'c1', 'c2', 'psi',
            'theta', 'Z', 'Y')

model.site.deadveg <- jags(model.file = 'ectomodSite1.txt',
                             data = datalist, n.chains = 3, 
                             parameters.to.save = params,
                             inits = inits, n.burnin = 4000, 
                             n.iter = 12000, n.thin = 10)

samples.site.deadveg <- jags.samples(model.site.deadveg$model, 
                                       variable.names = c("WAIC"),
                                       n.iter = 5000, n.burnin = 1000,
                                       n.thin = 1, type = "mean")

# Cov 2 only
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)), 
                 mass = hostmass, grassforb = comp.cov2,
                 sex = matrix(rep(hostsex, necto), ncol = necto),
                 julian = dates)

params <- c('a2', 'b1', 'b2', 'b3', 'c1', 'c2', 'psi',
            'theta', 'Z', 'Y')

model.site.grassforb <- jags(model.file = 'ectomodSite2.txt',
                           data = datalist, n.chains = 3, 
                           parameters.to.save = params,
                           inits = inits, n.burnin = 4000, 
                           n.iter = 12000, n.thin = 10)

samples.site.grassforb <- jags.samples(model.site.grassforb$model, 
                                     variable.names = c("WAIC"),
                                     n.iter = 5000, n.burnin = 1000,
                                     n.thin = 1, type = "mean")

# Writing WAIC table
WAIC.list <- list("Grass/Forb" = samples.site.grassforb[[1]])

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
              inits = inits, n.burnin = 4000, n.iter = 12000,
              n.thin = 10)

samples.host.full <- jags.samples(model$model, 
                                  variable.names = c("WAIC"),
                                  n.iter = 5000, n.burnin = 1000, 
                                  n.thin = 1, type = "mean")

# Intercept model
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, deadveg = comp.cov1, 
                 grassforb = comp.cov2, julian = dates)

params <- c('a1', 'a2', 'c1', 'c2', 'psi', 'theta', 'Z', 'Y')

model.host.inter <- jags(model.file = 'ectomod_HostIntercept.txt', 
              data = datalist, n.chains = 3, 
              parameters.to.save = params, inits = inits, 
              n.burnin = 4000, n.iter = 12000, n.thin = 10)

samples.host.inter <- jags.samples(model.host.inter$model, 
                                  variable.names = c("WAIC"),
                                  n.iter = 5000, n.burnin = 1000, 
                                  n.thin = 1, type = "mean")

# Host spec
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, host = hostspec, 
                 K = length(unique(hostspec)), deadveg = comp.cov1, 
                 grassforb = comp.cov2, julian = dates)

params <- c('a1', 'a2', 'b1', 'c1', 'c2', 'psi', 'theta', 'Z', 'Y')

model.host.spec <- jags(model.file = 'ectomodSpec.txt', 
                         data = datalist, n.chains = 3, 
                         parameters.to.save = params, inits = inits, 
                         n.burnin = 4000, n.iter = 7500, n.thin = 10)

samples.host.spec <- jags.samples(model.host.spec$model, 
                                   variable.names = c("WAIC"),
                                   n.iter = 5000, n.burnin = 1000, 
                                   n.thin = 1, type = "mean")

# Host mass
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2, julian = dates)

params <- c('a1', 'a2', 'b2', 'c1', 'c2', 'psi', 'theta', 'Z', 'Y')

model.host.mass <- jags(model.file = 'ectomodMass.txt', 
                        data = datalist, n.chains = 3, 
                        parameters.to.save = params, inits = inits, 
                        n.burnin = 4000, n.iter = 12000, n.thin = 10)

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
                 julian = dates)

params <- c('a1', 'a2', 'b3', 'c1', 'c2', 'psi', 'theta', 'Z', 'Y')

model.host.sex <- jags(model.file = 'ectomodSex.txt', 
                        data = datalist, n.chains = 3, 
                        parameters.to.save = params, inits = inits, 
                        n.burnin = 4000, n.iter = 12000, n.thin = 10)

samples.host.sex <- jags.samples(model.host.sex$model, 
                                  variable.names = c("WAIC"),
                                  n.iter = 5000, n.burnin = 1000, 
                                  n.thin = 1, type = "mean")

# Spec + mass
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat,
                 hostvec = hostvec, host = hostspec, 
                 K = length(unique(hostspec)), 
                 mass = hostmass, deadveg = comp.cov1, 
                 grassforb = comp.cov2, julian = dates)

params <- c('a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'psi',
            'theta', 'Z', 'Y')

model.host.specmass <- jags(model.file = 'ectomodSpecMass.txt', 
                       data = datalist, n.chains = 3, 
                       parameters.to.save = params, inits = inits, 
                       n.burnin = 4000, n.iter = 12000, n.thin = 10)

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
                 julian = dates)

params <- c('a1', 'a2', 'b1', 'b3', 'c1', 'c2', 'psi',
            'theta', 'Z', 'Y')

model.host.specsex <- jags(model.file = 'ectomodSpecSex.txt', 
                            data = datalist, n.chains = 3, 
                            parameters.to.save = params, 
                            inits = inits, n.burnin = 4000, 
                            n.iter = 12000, n.thin = 10)

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
                 julian = dates)

params <- c('a1', 'a2', 'b2', 'b3', 'c1', 'c2', 'psi',
            'theta', 'Z', 'Y')

model.host.sexmass <- jags(model.file = 'ectomodSexMass.txt', 
                           data = datalist, n.chains = 3, 
                           parameters.to.save = params, 
                           inits = inits, n.burnin = 4000, 
                           n.iter = 12000, n.thin = 10)

samples.host.sexmass <- jags.samples(model.host.sexmass$model, 
                                     variable.names = c("WAIC"),
                                     n.iter = 5000, n.burnin = 1000, 
                                     n.thin = 1, type = "mean")

# Set up WAIC table
WAIC.list <- list("SexMass" = samples.host.sexmass[[1]])

# Compare
waic.frame.host <- data.frame(model = names(WAIC.list), 
                              WAIC = unlist(lapply(WAIC.list, sum)))

# load existing WAIC table
waic.table <- read.csv(file = "waic_host.csv")

# combine with other WAICs and write table
waic.frame.host <- rbind(waic.frame.host, waic.table)

waic.frame.host <- arrange(waic.frame.host, WAIC)

print(waic.frame.host)

# write.csv(waic.frame.host, row.names = F,
#           file = "waic_host.csv")

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

ggplot(data = a1s, aes(x = rownames(a1s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_point(aes(x = rownames(a1s), y = 7, color = sig), shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  geom_hline(yintercept = 0)+
  labs(x = "Ectoparasite Species", y = "Veg Composition Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None")

# ggsave(filename = "vegforest.jpeg", width = 9,
#        height = 4,units = "in")

# Veg comp: grass or forb
a2 <- model$BUGSoutput$sims.list$a2

a2s <- data.frame(mean = apply(a2, 2, mean),
           lo = apply(a2, 2, quantile, 0.125),
           hi = apply(a2, 2, quantile, 0.875)) %>%
      mutate(sig = case_when(lo > 0 | hi < 0 ~ "Yes",
                         TRUE ~ ""))
rownames(a2s) <- colnames(site.occ)[-c(1)]

ggplot(data = a2s, aes(x = rownames(a2s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = rownames(a1s), y = 2, color = sig), shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Vertical Structure Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None")

# ggsave(filename = "vegforb.jpeg", width = 9, height = 4,
#        units = "in")

# Effect of host spec
b1 <- model$BUGSoutput$sims.list$b1

b1lo <- as.data.frame(apply(b1, c(2,3), quantile, 0.125))
colnames(b1lo) <- levels(hostspec)
b1lo$Ecto <- sort(unique(mecto.caps$Ecto))#[-c(1,4)]

b1lo <- pivot_longer(b1lo, cols = -Ecto, 
                     names_to = "host", values_to = "lo")

b1hi <- as.data.frame(apply(b1, c(2,3), quantile, 0.875))
colnames(b1hi) <- levels(hostspec)
b1hi$Ecto <- sort(unique(mecto.caps$Ecto))#[-c(1,4)]

b1hi <- pivot_longer(b1hi, cols = -Ecto, 
                     names_to = "host", values_to = "hi")

b1mean <- as.data.frame(apply(b1, c(2,3), mean))
colnames(b1mean) <- levels(hostspec)
b1mean$Ecto <- sort(unique(mecto.caps$Ecto))#[-c(1,4)]

b1mean <- pivot_longer(b1mean, cols = -Ecto,
                       names_to = "host", values_to = "mean")

hilo <- inner_join(b1hi, b1lo, by = c("host", "Ecto"))
hilo$mean <- b1mean$mean

hilo[which(hilo$lo > 0 | hilo$hi < 0),]

# Make interval plot for each ectoparasite species:
ecto.specs <- unique(hilo$Ecto)

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
  geom_point(aes(x = rownames(b2s), y = -1, color = sig), shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Host Mass Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None")

# ggsave(filename = "hostmass.jpeg", width = 9, height = 4,
#        units = "in")

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
  geom_point(aes(x = rownames(b3s), y = -1.1, color = sig), 
             shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Host Sex Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None")

# ggsave(filename = "hostsex.jpeg", width = 9, height = 4,
#        units = "in")

# Effect of capture no
c1 <- model$BUGSoutput$sims.list$c1

c1s <- data.frame(mean = apply(c1, 2, mean), 
           lo = apply(c1, 2, quantile, 0.125), 
           hi = apply(c1, 2, quantile, 0.875)) %>%
  mutate(sig = NA) %>%
  mutate(sig = case_when(lo > 0 | hi < 0 ~ "Yes",
                         TRUE ~ ""))
rownames(c1s) <- colnames(site.occ)[-c(1)]

ggplot(data = c1s, aes(x = rownames(c1s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = rownames(c1s), y = 0.8, color = sig), shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Capture no. Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None")

# ggsave(filename = "capno.jpeg", width = 9, height = 4,
#        units = "in")

# Julian date
c2 <- model$BUGSoutput$sims.list$c2

c2s <- data.frame(mean = apply(c2, 2, mean), 
           lo = apply(c2, 2, quantile, 0.125), 
           hi = apply(c2, 2, quantile, 0.875)) %>%
  mutate(sig = case_when(lo > 0 | hi < 0 ~ "Yes",
                         TRUE ~ ""))
rownames(c2s) <- colnames(site.occ)[-c(1)]

ggplot(data = c2s, aes(x = rownames(c2s), y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = rownames(c2s), y = 1, color = sig), shape = 8)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Ectoparasite Species", y = "Julian Date Coefficient")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1.1),
        panel.grid = element_blank(), legend.position = "None")

# ggsave(filename = "datecov.jpeg", width = 9, height = 4,
#        units = "in")

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
ecto.specs <- names(ecto.list2)#[-c(1,4)]

# Site-level r-squared
mean.site.r2 <- apply(site.r2, 1, mean)
median.site.r2 <- apply(site.r2, 1, median)
site.df <- data.frame(Ecto = ecto.specs, site.r2 = mean.site.r2,
                      median.r2 = median.site.r2)

# Get family and order names for each species
site.df <- site.df %>%
  left_join(n.host, by = c("Ecto" = "ecto"))

# If counting host by subfamily:
# site.df <- site.df %>%
#   left_join(n.host.sub, by = c("Ecto" = "ecto"))

# If counting host by family:
# site.df <- site.df %>%
#   left_join(n.host.fam, by = c("Ecto" = "ecto"))

# Flea classifications
fleaclass <- read.csv("FleaClassifications.csv") %>%
  unite(col = "Ecto", Genus, Species, sep = "_") %>%
  select(Ecto, Classification)
# fleaclass <- fleaclass[-2,]

# Merge with site df
site.df <- site.df %>%
  left_join(fleaclass, by = "Ecto") %>%
  mutate(Order = case_when(startsWith(Ecto, "Ixodes") ~ "Ixodida",
                           TRUE ~ Order)) %>%
  mutate(Classification = case_when(Order == "Ixodida"  ~ "Ephemeral",
                                    Order == "Mesostigmata" 
                                    ~ "Nest",
                                    Order == "Diptera" ~ "Diptera",
                                    Order == "Psocodea" ~ "Fur",
                                    TRUE ~ Classification)) %>%
  mutate(hosts = as.numeric(hosts)) %>%
  mutate(hosts = case_when(Ecto == "Ixodes_scapularis_larva" ~ 5,
                           Ecto == "Ixodes_scapularis_nymph" ~ 5,
                           TRUE ~ hosts))

# Stat test
kruskal.test(formula = site.r2~Classification, 
             data = site.df[site.df$Classification != "Diptera",])
# Kruskal because variances and sample sizes are different

# Eta-sq (effect size)
kruskal_effsize(data = site.df[site.df$Classification != "Diptera",],
                formula = site.r2~Classification)

# dunn test
dunn_test(data = site.df[site.df$Classification != "Diptera",],
          formula = site.r2~Classification)

# full fig
siteR2fig <- ggplot(data = site.df, aes(x = Classification, 
                                        y = site.r2))+
  geom_boxplot(fill = "lightgray", outlier.shape = NA)+
  geom_jitter(aes(color = Order))+
  scale_color_viridis_d(end = 0.95)+
  labs(y = bquote("Site"~R^2))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), legend.position = "none")

# ggsave(filename = "siter2.jpeg", height = 4, width = 6,
#        units = "in", dpi = 600)

# Host-level r-squared
mean.host.r2 <- apply(host.r2, 1, mean)
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
  left_join(fleaclass, by = "Ecto") %>%
  mutate(Order = case_when(startsWith(Ecto, "Ixodes") ~ "Ixodida",
                           TRUE ~ Order)) %>%
  mutate(Classification = case_when(Order == "Ixodida"  ~ "Ephemeral",
                                    Order == "Mesostigmata" ~ 
                                      "Nest",
                                    Order == "Diptera" ~ "Diptera",
                                    Order == "Psocodea" ~ "Fur",
                                    TRUE ~ Classification)) %>%
  mutate(hosts = as.numeric(hosts)) %>%
  mutate(hosts = case_when(Ecto == "Ixodes_scapularis_larva" ~ 5,
                           Ecto == "Ixodes_scapularis_nymph" ~ 5,
                           TRUE ~ hosts))

# Stat test
kruskal.test(formula = host.r2~Classification, 
             data = host.df[host.df$Classification != "Diptera",])
# Kruskal because variances and sample sizes are different
# No Diptera because there is 1 species, 2 total

# Eta-sq (effect size)
kruskal_effsize(data = host.df[host.df$Classification != "Diptera",],
                formula = host.r2~Classification)
dunn_test(data = host.df[host.df$Classification != "Diptera",],
          formula = host.r2~Classification)

# Final fig
hostR2fig <- ggplot(data = host.df, aes(x = Classification, 
                                        y = host.r2))+
  geom_boxplot(fill = 'lightgray', outlier.shape = NA)+
  geom_jitter(aes(color = Order))+
  scale_color_viridis_d(end = 0.95)+
  labs(y = bquote("Host"~R^2))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "hostr2.jpeg", height = 4, width = 6,
#        units = 'in', dpi = 600)

r2fig <- siteR2fig + hostR2fig +
  plot_annotation(tag_levels = "a")

# ggsave(filename = "r2fig.jpeg", height = 3, width = 10, 
#        units = 'in',dpi = 600)

# R-squared vs. host specificity -------------------
summary(lm(data = site.df, formula = hosts~site.r2))
summary(lm(data = host.df, formula = hosts~host.r2))
# This is wild

r2site <- ggplot(data = site.df, aes(x = site.r2, y = hosts))+
  geom_point(aes(color = Classification), size = 2)+
  # geom_smooth(se = F, method = 'lm', color = "black")+
  scale_color_viridis_d(end = 0.95)+
  labs(x = bquote("Site"~R^2), y = "Number of Hosts")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("siter2hosts.jpeg", height =3, width = 5,
#        units = "in")
# ggsave("siter2hostsub.jpeg", height =3, width = 5, units = "in")
# ggsave("siter2hostfam.jpeg", height =3, width = 5, units = "in")

r2host <- ggplot(data = host.df, aes(x = host.r2, y = hosts))+
  geom_point(aes(color = Classification), size = 2)+
  geom_smooth(se = F, method = 'lm', color = "black")+
  scale_color_viridis_d(end = 0.95)+
  labs(x = bquote("Host"~R^2), y = "Number of Hosts")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("hostr2hosts.jpeg", height = 3, width = 5,
#        units = "in")
# ggsave("hostr2hostsub.jpeg", height = 3, width = 5, units = "in")
# ggsave("hostr2hostfam.jpeg", height = 3, width = 5, units = "in")

r2site + r2host +
  plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "a")

# ggsave(filename = "r2nhost.jpeg", width = 10, height = 3, 
#        units = "in", dpi = 600)

# Hosts per classification ------------------
# Stat test
kruskal.test(formula = hosts~Classification, 
             data = site.df[site.df$Classification != "Diptera",])
# Kruskal because variances and sample sizes are different

# Eta-sq (effect size)
kruskal_effsize(data = site.df[site.df$Classification != "Diptera",],
                formula = hosts~Classification)

# Dunn test
dunn_test(data = site.df[site.df$Classification != "Diptera",],
          formula = hosts~Classification)

ggplot(data = site.df, aes(x = Classification, y = hosts))+
  geom_boxplot(fill = 'lightgray')+
  geom_jitter(aes(color = Order))+
  geom_text(data = data.frame(labs = unique(site.df$Classification),
                              groups = c("b", " ", "a", "b")),
            aes(x = labs, y = 6.75, label = groups))+
  scale_color_viridis_d(end = 0.95)+
  labs(y = "Number of Host Species")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("numhosts.jpeg", height = 3, width = 5,
#        units = "in")

# Get ecto abundances on hosts for prefs -----------------
# Get parasite abund for each host
ecto.host.abund <- mamm.ecto %>%
  select(Abbrev, Order, Family, Genus, Species, Other) %>%
  filter(is.na(Order) == F) %>%
  mutate(Species = case_when(Genus == "Ixodes" & 
                               Species == "scapularis" ~
                               paste(Species, Other, sep = "_"),
                             TRUE ~ Species)) %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
  select(Abbrev, Ecto) %>%
  group_by(Abbrev, Ecto) %>%
  summarise(spec_abund = n())

# Make figs for raw abundances
for(i in 1:length(unique(ecto.host.abund$Ecto))){
  sub <- ecto.host.abund[ecto.host.abund$Ecto == unique(ecto.host.abund$Ecto)[i],]
  
  plt <- ggplot(data = sub, 
                aes(x = fct_reorder(Abbrev, desc(spec_abund)), 
                         y = spec_abund))+
    geom_col()+
    labs(x = "Hosts", y = "Abundance", 
         title = unique(ecto.host.abund$Ecto)[i])+
    theme_bw()+
    theme(panel.grid = element_blank())
  
  print(plt)
}  

# Get host relative abundances
ecto.host.rel <- mamm.clean %>%
  mutate(Abund = rowSums(select(., `1`:`9`), na.rm = T)) %>%
  group_by(Abbrev) %>%
  summarise(Abund = sum(Abund)) %>%
  ungroup() %>%
  mutate(RelAbund = Abund/sum(Abund)) %>%
  select(-Abund) %>%
  right_join(y = ecto.host.abund, by = "Abbrev") %>%
  mutate(rel_spec_abund = (1-RelAbund)*spec_abund)

for(i in 1:length(unique(ecto.host.rel$Ecto))){
  sub <- ecto.host.rel[ecto.host.rel$Ecto == unique(ecto.host.rel$Ecto)[i],]
  
  plt <- ggplot(data = sub, 
                aes(x = fct_reorder(Abbrev, desc(rel_spec_abund)), 
                    y = rel_spec_abund))+
    geom_col()+
    labs(x = "Hosts", y = "Relative Abundance", 
         title = unique(ecto.host.rel$Ecto)[i])+
    theme_bw()+
    theme(panel.grid = element_blank())
  
  print(plt)
}

# Host occupancy model -------------------------
# Coerce mammal data to wide format
mamm.clean <- mamm.ecto %>%
  select(Site, Day, Abbrev) %>%
  mutate(Occ = 1) %>%
  complete(Site, Day, Abbrev) %>%
  pivot_wider(names_from = "Day", values_from = "Occ", 
              values_fn = list(Occ = "sum"), values_fill = 0) 

# Coerce to array
mamm.array <- abind(split(mamm.clean, mamm.clean$Abbrev), 
      along = 3)
mamm.array <- mamm.array[,-c(1:2),]
rownames(mamm.array) <- unique(mamm.clean$Site)

# Clean it up
mamm.array[is.na(mamm.array)] <- 0
mamm.array <- array(as.numeric(mamm.array), dim = dim(mamm.array))
mamm.array[mamm.array > 1] <- 1

# New date covariate
missing.dates <- read.csv("Nocap_sites.csv") %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y"))

mamm.date <- mamm.ecto %>%
  full_join(y = missing.dates, by = c("Site", "Date", "Day")) %>%
  mutate(julian = format(Date, "%j")) %>%
  mutate(julian = as.vector(scale(as.numeric(julian)))) %>%
  dplyr::select(Site, Day, julian) %>%
  distinct(Site, Day, .keep_all = T) %>%
  pivot_wider(names_from = Day, values_from = julian) %>%
  arrange(Site) %>%
  relocate(Site, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`) 

date.cov <- as.matrix(mamm.date[,-1], nrow = nrow(mamm.date),
                      ncol = ncol(mamm.date)-1)

# Write MSOM
cat("
    model{
    # Hyperpriors
    a0.mean ~ dunif(0,1)
    mean.a0 <- log(a0.mean)-log(1-a0.mean)
    tau.a0 ~ dgamma(0.1, 0.1)
    
    a1.mean ~ dunif(0,1)
    mean.a1 <- log(a1.mean)-log(1-a1.mean)
    tau.a1 ~ dgamma(0.1, 0.1)
    
    a2.mean ~ dunif(0,1)
    mean.a2 <- log(a2.mean)-log(1-a2.mean)
    tau.a2 ~ dgamma(0.1, 0.1)
    
    b0.mean ~ dunif(0,1)
    mean.b0 <- log(b0.mean)-log(1-b0.mean)
    tau.b0 ~ dgamma(0.1, 0.1)
    
    b1.mean ~ dunif(0,1)
    mean.b1 <- log(b1.mean)-log(1-b1.mean)
    tau.b1 ~ dgamma(0.1, 0.1)
    
    # Species-level priors
    for(i in 1:nspec){
      
      a0[i] ~ dnorm(mean.a0, tau.a0)
      a1[i] ~ dnorm(mean.a1, tau.a1)
      a2[i] ~ dnorm(mean.a2, tau.a2)
      
      b0[i] ~ dnorm(mean.b0, tau.b0)
      b1[i] ~ dnorm(mean.b1, tau.b1)
      
      # Occupancy
      for(j in 1:nsite){
        logit(psi[j,i]) <- a0[i] + a1[i]*deadveg[j] + 
                            a2[i]*grassforb[j]
        Z[j,i] ~ dbern(psi[j,i])
        
       # Detection
       for(k in 1:nday){
         logit(p[j,k,i]) <- b0[i] + b1[i]*date[j,k]
         obs[j,k,i] ~ dbern(p[j,k,i]*Z[j,i])
        }
      }
    }
    }", file = "HostMSOM.txt")

# Set up and run
nspec <- dim(mamm.array)[3]
nsite <- dim(mamm.array)[1]
nday <- dim(mamm.array)[2]

data <- list(nspec = nspec, nsite = nsite, nday = nday, 
             deadveg = comp.cov1, grassforb = comp.cov2,
             date = date.cov, obs = mamm.array)

# Get observed presence-absence matrix
maxobs <- apply(mamm.array, c(1,3), max)

params <- c('a1', 'a2', 'b1', 'psi', 'Z')

# Init values
inits <- function(){
  inits <- list(Z = maxobs)
}

# smamm.model <- jags(model.file = 'HostMSOM.txt', data = data,
#               n.chains = 3, parameters.to.save = params,
#               inits = inits, n.burnin = 1000, n.iter = 5000,
#               n.thin = 5)
# 
# saveRDS(smamm.model, file = "HostMod.rds")

smamms <- readRDS(file = "HostMod.rds")

# Get dead veg coefs
mamma1 <- as.data.frame(smamms$BUGSoutput$sims.list$a1)
colnames(mamma1) <- unique(mamm.clean$Abbrev)

mamm.plt1 <- mamma1 %>%
  pivot_longer(cols = everything(), names_to = "MammSpec", 
               values_to = "Output") %>%
  group_by(MammSpec) %>%
  summarise(meana1 = mean(Output), lo = quantile(Output, 0.125),
            hi = quantile(Output, 0.875))

ggplot(data = mamm.plt1, aes(x = MammSpec, y = meana1))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(x = "Host Species", y = "Dead Vegetation Coefficient")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# get grass/forb coef
mamma2 <- as.data.frame(smamms$BUGSoutput$sims.list$a2)
colnames(mamma2) <- unique(mamm.clean$Abbrev)

mamm.plt2 <- mamma2 %>%
  pivot_longer(cols = everything(), names_to = "MammSpec", 
               values_to = "Output") %>%
  group_by(MammSpec) %>%
  summarise(meana2 = mean(Output), lo = quantile(Output, 0.125),
            hi = quantile(Output, 0.875))

ggplot(data = mamm.plt2, aes(x = MammSpec, y = meana2))+
  geom_point()+
  geom_errorbar(aes(ymin = lo, ymax = hi))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(x = "Host Species", y = "Grass/Forb Coefficient")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())
# mamms mostly differ in forest vs not preferences

# merge ectos and smamms
a1s$Parasite <- rownames(a1s)
ecto.prefs <- read.csv(file = "HostPrefs_raw.csv")
# ecto.prefs <- read.csv(file = "HostPrefs_Adj.csv")

cov.resps <- full_join(mamm.plt1, ecto.prefs, 
          by = c("MammSpec" = "PrefMammal")) %>%
  filter(is.na(Parasite) == F) %>%
  left_join(y = a1s, by = "Parasite") %>%
  select(-c(lo.x, hi.x, lo.y, hi.y, sig)) %>%
  left_join(y = site.df, by = c("Parasite" = "Ecto")) %>%
  select(MammSpec:mean, Classification) %>%
  mutate(Habitat = case_when(MammSpec %in% c("ZAHU", "MIPE") ~ "Open",
                             MammSpec %in% c("PEMA", "MYGA")~"Forest",
                             TRUE ~ "Generalist"))

# Test host/parasite coef associations
summary(lm(data = cov.resps, mean~meana1))
kruskal.test(data = cov.resps, mean~Habitat)
dunn_test(data = cov.resps, mean~Habitat)
ggplot(data = cov.resps, aes(x = Habitat, y = mean))+
  geom_boxplot(fill = "lightgray")+
  # geom_point()+
  # geom_smooth(method = 'lm', se = F, color = 'black')+
  labs(x = "Host Habitat", y = "Parasite Coefficient")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "hostcoefbox_adj.jpeg", width = 5, height = 3,
#        units = "in")

# Do it with primary hosts only
summary(lm(data = cov.resps[cov.resps$Secondary != "Secondary_host",], mean~meana1))
t.test(data = cov.resps[cov.resps$Secondary != "Secondary_host",], mean~Habitat)
cohens_d(data = cov.resps[cov.resps$Secondary != "Secondary_host",],
          mean~Habitat)

ggplot(data = cov.resps[cov.resps$Secondary != "Secondary_host",], 
       aes(x = meana1, y = mean))+
  # geom_boxplot(fill = "lightgray")+
  geom_point()+
  # geom_smooth(method = 'lm', se = F, color = 'black')+
  labs(x = "Host Coefficient", y = "Parasite Coefficient")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "hostcoefregress_prim_raw.jpeg", width = 5,
#        height = 3, units = "in")

# Getting P. leucopus ectos ---------------------
pele.ectos <- mamm.ecto %>%
  select(Abbrev, Site, Order, Family, Genus, Species, Other) %>%
  filter(is.na(Order) == F) %>%
  mutate(Species = case_when(Genus == "Ixodes" & 
                               Species == "scapularis" ~
                               paste(Species, Other, sep = "_"),
                             TRUE ~ Species)) %>%
  unite(col = "Ecto", Genus, Species, sep = "_", na.rm = T) %>%
  select(Abbrev, Site, Ecto) %>%
  group_by(Abbrev, Site, Ecto) %>%
  summarise(spec_abund = n()) %>%
  filter(Abbrev == "PELE")
