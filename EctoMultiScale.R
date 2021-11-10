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
  filter(Tag != "RECAP" & Tag != "DESC")

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

# Vector of hosts per site
hostvec <- apply(tagmat, 2, function(x) length(na.omit(x)))

ecto.list2 <- lapply(ecto.list, 
                    function(x) x <- select(x, -c(Site, Ecto, Tag)))

# Coerce to array
ecto.ar <- array(as.numeric(unlist(ecto.list2)), 
                 dim=c(nrow(ecto.list2[[1]]), ncol(ecto.list2[[1]]), 
                       length(ecto.list2)))

# Site covariate: Vertical complexity -----------------------

# Site covariate: Cover types -------------------------

# Host covariate: Species -----------------------
# Get species abbrev of each host
host.spec <- mecto.clean %>%
  select(Tag, Abbrev) %>%
  distinct() 

# Convert to factor
hostspec <- factor(host.spec$Abbrev, levels = unique(host.spec$Abbrev))

# Host covariate: Adj. body mass -----------------------

# Host covariate: Sex -------------------------

# Detection covariate: Julian date? -----------------------

# Write model -------------------------
cat("
    model{
    # Define hyperpriors
    
    mean.a0 ~ dunif(0,1)
    a0.mean <- log(mean.a0)-log(1-mean.a0)
    tau.a0 ~ dgamma(0.1, 0.1)
    
    mean.b0 ~ dunif(0,1)
    b0.mean <- log(mean.b0)-log(1-mean.b0)
    tau.b0 ~ dgamma(0.1, 0.1)
    
    mean.c0 ~ dunif(0,1)
    c0.mean <- log(mean.c0)-log(1-mean.c0)
    tau.c0 ~ dgamma(0.1, 0.1)
    
    mean.c1 ~ dunif(0,1)
    c1.mean <- log(mean.c1)-log(1-mean.c1)
    tau.c1 ~ dgamma(0.1, 0.1)
    
    # Priors
    for(i in 1:necto){
    
      a0[i] ~ dnorm(a0.mean, tau.a0)
      
      b0[i] ~ dnorm(b0.mean, tau.b0)
      b1[i,1] <- 0
      for(s in 2:K){
        b1[i,s] ~ dnorm(0, 1)
      }
      
      c0[i] ~ dnorm(c0.mean, tau.c0)
      c1[i] ~ dnorm(c1.mean, tau.c1)
    
      # Occupancy model: Site
      for(j in 1:nsite){
      
        logit(psi[j,i]) <- a0[i]
        Z[j,i] ~ dbern(psi[j,i])

        # Occupancy model: Host
        for(k in tagmat[1:hostvec[j],j]){
        
          logit(theta[k,i]) <- b0[i] + b1[i,host[k]]
          Y[k,i] ~ dbern(theta[k,i]*Z[j,i])
    
          # Detection model
          for(l in 1:ncap[k]){
          
          logit(p[k,l,i]) <- c0[i] + c1[i]*l
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

# Define data to send and params to keep
datalist <- list(necto = necto, nsite = nsite, ncap = ncap,
                 obs = ecto.ar, tagmat = tagmat, hostvec = hostvec, 
                 host = hostspec, K = length(unique(hostspec)))

params <- c('a0', 'b0', 'b1', 'c0', 'c1')

# Init values
site.occ <- mecto.smol %>%
  dplyr::select(Site, Ecto, Occ) %>%
  filter(is.na(Ecto) == F) %>%
  distinct() %>%
  arrange(Site) %>%
  pivot_wider(names_from = Ecto, values_from = Occ)

site.occ[is.na(site.occ)] <- 0
site.occ <- site.occ[,order(colnames(site.occ))]
sitemax <- as.matrix(site.occ[,-16])

# Known presences per host
hostmax <- apply(ecto.ar, c(1,3), max, na.rm = T)

inits <- function(){
  inits <- list(
    Z = sitemax,
    Y = hostmax
  )
}

# Send model to JAGS
model <- jags(model.file = 'ectomod.txt', data = datalist, n.chains = 3,
              parameters.to.save = params, inits = inits, n.burnin = 4000,
              n.iter = 10000, n.thin = 3)

# Save model
saveRDS(model, file = "ectomod.rds")

# Prelim results ----------------------
# Effect of host spec
b1 <- model$BUGSoutput$sims.list$b1

b1lo <- as.data.frame(apply(b1, c(2,3), quantile, 0.025))
colnames(b1lo) <- levels(hostspec)
b1lo$Ecto <- sort(unique(mecto.caps$Ecto))

b1lo <- pivot_longer(b1lo, cols = -Ecto, 
                     names_to = "host", values_to = "lo")

b1hi <- as.data.frame(apply(b1, c(2,3), quantile, 0.975))
colnames(b1hi) <- levels(hostspec)
b1hi$Ecto <- sort(unique(mecto.caps$Ecto))

b1hi <- pivot_longer(b1hi, cols = -Ecto, 
                     names_to = "host", values_to = "hi")

hilo <- inner_join(b1hi, b1lo, by = c("host", "Ecto"))
# will need to re-run model for each host spec but at least
# prelim stuff makes sense

# Effect of capture no.
c1 <- model$BUGSoutput$sims.list$c1
quantile(c1, c(0.025, 0.975))
capno <- data.frame(mean = apply(c1, 2, mean), 
                    lo = apply(c1, 2, quantile, 0.025), 
                    hi = apply(c1, 2, quantile, 0.975))
# doesn't appear to affect detectability
