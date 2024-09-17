Experimental study (Tadpoles)
================
Rodolfo Pelinson
2024-09-17

# Mass, stage and individual performance analysis of generalist amphibians for different detritus type.

``` r
library(glmmTMB)
library(bbmle)
library(DHARMa)
```

Loading data

``` r
path <- c("C:/Users/rodol/OneDrive/repos/Detritus_generalist_amphibians/")

data_table <- read.csv(paste(path, "data/tadpoles_data.csv", sep = ""))

data_table$Treatment <- factor(data_table$Leaf.Litter.Treatment, levels = c("Cerrado","Pasture","Sugarcane"))
data_table$Aquarium <- factor(data_table$Aquarium)


Treatment_Pasture_Sugarcane <- rep(NA, nrow(data_table))
for(i in 1:nrow(data_table)){
  if(data_table$Treatment[i] == "Pasture" | data_table$Treatment[i] == "Sugarcane"){
    Treatment_Pasture_Sugarcane[i] <- "Pasture-Sugarcane"
  }else{Treatment_Pasture_Sugarcane[i] <- as.character(data_table$Treatment[i])}
}

Treatment_Pasture_Cerrado <- rep(NA, nrow(data_table))
for(i in 1:nrow(data_table)){
  if(data_table$Treatment[i] == "Pasture" | data_table$Treatment[i] == "Cerrado"){
    Treatment_Pasture_Cerrado[i] <- "Pasture-Cerrado"
  }else{Treatment_Pasture_Cerrado[i] <- as.character(data_table$Treatment[i])}
}

Treatment_Sugarcane_Cerrado <- rep(NA, nrow(data_table))
for(i in 1:nrow(data_table)){
  if(data_table$Treatment[i] == "Sugarcane" | data_table$Treatment[i] == "Cerrado"){
    Treatment_Sugarcane_Cerrado[i] <- "Sugarcane-Cerrado"
  }else{Treatment_Sugarcane_Cerrado[i] <- as.character(data_table$Treatment[i])}
}

data_table <- data.frame(data_table, Treatment_Pasture_Sugarcane, Treatment_Pasture_Cerrado,Treatment_Sugarcane_Cerrado)


data_table_nattereri <- data_table[which(data_table$Species == "P. nattereri"),]
data_table_scinax <- data_table[which(data_table$Species == "S. fuscovarius"),]
data_table_cuvieri <- data_table[which(data_table$Species == "P. cuvieri"),]
data_table_centralis <- data_table[which(data_table$Species == "P. centralis"),]

data_table <- rbind(data_table_centralis ,data_table_cuvieri, data_table_nattereri, data_table_scinax)
```

# Mass

## S. fuscovarius

``` r
mass_scinax_no_effect <- glmmTMB(Tadpole.Mass~ 1 + (1|Aquarium), data=data_table_scinax, family=Gamma(link = "log"))
mass_scinax_land <- glmmTMB(Tadpole.Mass~ Treatment + (1|Aquarium), data=data_table_scinax, family=Gamma(link = "log"))
mass_scinax_Pasture_Sugarcane_model <- glmmTMB(Tadpole.Mass ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=data_table_scinax, family=Gamma(link = "log"))
mass_scinax_Pasture_Cerrado_model <- glmmTMB(Tadpole.Mass~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table_scinax, family=Gamma(link = "log"))
mass_scinax_Sugarcane_Cerrado_model <- glmmTMB(Tadpole.Mass~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table_scinax, family=Gamma(link = "log"))


AICTAB_mass_scinax <- AICctab(mass_scinax_no_effect, mass_scinax_land, mass_scinax_Pasture_Sugarcane_model,mass_scinax_Pasture_Cerrado_model,mass_scinax_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                 mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_mass_scinax
```

    ##                               logLik  AICc    dLogLik dAICc   df weight
    ## No Effect                     -1479.7  2965.5     0.0     0.7 3  0.251 
    ## Savanna # Pasture # Sugarcane -1478.3  2966.8     1.4     1.9 5  0.133 
    ## Savanna # Pasture = Sugarcane -1478.4  2964.9     1.4     0.0 4  0.348 
    ## Savanna = Pasture # Sugarcane -1479.6  2967.4     0.1     2.5 4  0.099 
    ## Savanna = Sugarcane # Pasture -1479.1  2966.3     0.6     1.5 4  0.168

## P. cuvieri

``` r
mass_cuvieri_no_effect <- glmmTMB(Tadpole.Mass~ 1 + (1|Aquarium), data=data_table_cuvieri, family=Gamma(link = "log"))
mass_cuvieri_land <- glmmTMB(Tadpole.Mass~ Treatment + (1|Aquarium), data=data_table_cuvieri, family=Gamma(link = "log"))
mass_cuvieri_Pasture_Sugarcane_model <- glmmTMB(Tadpole.Mass ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=data_table_cuvieri, family=Gamma(link = "log"))
mass_cuvieri_Pasture_Cerrado_model <- glmmTMB(Tadpole.Mass~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table_cuvieri, family=Gamma(link = "log"))
mass_cuvieri_Sugarcane_Cerrado_model <- glmmTMB(Tadpole.Mass~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table_cuvieri, family=Gamma(link = "log"))


AICTAB_mass_cuvieri <- AICctab(mass_cuvieri_no_effect, mass_cuvieri_land, mass_cuvieri_Pasture_Sugarcane_model,mass_cuvieri_Pasture_Cerrado_model,mass_cuvieri_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                  mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_mass_cuvieri
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -902.2 1810.6    0.0    15.0 3  <0.001
    ## Savanna # Pasture # Sugarcane -892.7 1795.6    9.6     0.0 5  0.904 
    ## Savanna # Pasture = Sugarcane -896.3 1800.8    6.0     5.1 4  0.070 
    ## Savanna = Pasture # Sugarcane -902.2 1812.7    0.0    17.0 4  <0.001
    ## Savanna = Sugarcane # Pasture -897.3 1802.8    4.9     7.2 4  0.025

## P. nattereri

``` r
mass_nattereri_no_effect <- glmmTMB(Tadpole.Mass~ 1 + (1|Aquarium), data=data_table_nattereri, family=Gamma(link = "log"))
mass_nattereri_land <- glmmTMB(Tadpole.Mass~ Treatment + (1|Aquarium), data=data_table_nattereri, family=Gamma(link = "log"))
mass_nattereri_Pasture_Sugarcane_model <- glmmTMB(Tadpole.Mass ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=data_table_nattereri, family=Gamma(link = "log"))
mass_nattereri_Pasture_Cerrado_model <- glmmTMB(Tadpole.Mass~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table_nattereri, family=Gamma(link = "log"))
mass_nattereri_Sugarcane_Cerrado_model <- glmmTMB(Tadpole.Mass~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table_nattereri, family=Gamma(link = "log"))


AICTAB_mass_nattereri <- AICctab(mass_nattereri_no_effect, mass_nattereri_land, mass_nattereri_Pasture_Sugarcane_model,mass_nattereri_Pasture_Cerrado_model,mass_nattereri_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                    mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_mass_nattereri
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -957.7 1921.6    0.0    11.1 3  0.002 
    ## Savanna # Pasture # Sugarcane -950.1 1910.5    7.7     0.0 5  0.513 
    ## Savanna # Pasture = Sugarcane -951.2 1910.7    6.5     0.2 4  0.472 
    ## Savanna = Pasture # Sugarcane -957.5 1923.2    0.3    12.6 4  <0.001
    ## Savanna = Sugarcane # Pasture -954.9 1918.1    2.8     7.6 4  0.012

## P. centralis

``` r
mass_centralis_no_effect <- glmmTMB(Tadpole.Mass~ 1 + (1|Aquarium), data=data_table_centralis, family=Gamma(link = "log"))
mass_centralis_land <- glmmTMB(Tadpole.Mass~ Treatment + (1|Aquarium), data=data_table_centralis, family=Gamma(link = "log"))
mass_centralis_Pasture_Sugarcane_model <- glmmTMB(Tadpole.Mass ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=data_table_centralis, family=Gamma(link = "log"))
mass_centralis_Pasture_Cerrado_model <- glmmTMB(Tadpole.Mass~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table_centralis, family=Gamma(link = "log"))
mass_centralis_Sugarcane_Cerrado_model <- glmmTMB(Tadpole.Mass~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table_centralis, family=Gamma(link = "log"))


AICTAB_mass_centralis <- AICctab(mass_centralis_no_effect, mass_centralis_land, mass_centralis_Pasture_Sugarcane_model,mass_centralis_Pasture_Cerrado_model,mass_centralis_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                    mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_mass_centralis 
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -287.2  581.0    0.0     2.3 3  0.114 
    ## Savanna # Pasture # Sugarcane -284.1  579.8    3.1     1.1 5  0.210 
    ## Savanna # Pasture = Sugarcane -285.1  579.2    2.1     0.5 4  0.277 
    ## Savanna = Pasture # Sugarcane -284.8  578.7    2.4     0.0 4  0.364 
    ## Savanna = Sugarcane # Pasture -287.2  583.4    0.0     4.7 4  0.034

# Stage

## S. fuscovarius

``` r
stage_scinax_no_effect <- glmmTMB(Gosner.Stage~ 1 + (1|Aquarium), data=data_table_scinax, family=gaussian(link = "identity"))
stage_scinax_land <- glmmTMB(Gosner.Stage~ Treatment + (1|Aquarium), data=data_table_scinax, family=gaussian(link = "identity"))
stage_scinax_Pasture_Sugarcane_model <- glmmTMB(Gosner.Stage ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=data_table_scinax, family=gaussian(link = "identity"))
stage_scinax_Pasture_Cerrado_model <- glmmTMB(Gosner.Stage~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table_scinax, family=gaussian(link = "identity"))
stage_scinax_Sugarcane_Cerrado_model <- glmmTMB(Gosner.Stage~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table_scinax, family=gaussian(link = "identity"))


AICTAB_stage_scinax <- AICctab(stage_scinax_no_effect, stage_scinax_land, stage_scinax_Pasture_Sugarcane_model,stage_scinax_Pasture_Cerrado_model,stage_scinax_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                    mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_stage_scinax
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -460.9  928.0    0.0     3.0 3  0.094 
    ## Savanna # Pasture # Sugarcane -457.8  926.0    3.1     1.0 5  0.260 
    ## Savanna # Pasture = Sugarcane -458.4  925.0    2.5     0.0 4  0.418 
    ## Savanna = Pasture # Sugarcane -460.9  930.0    0.0     5.0 4  0.035 
    ## Savanna = Sugarcane # Pasture -459.2  926.5    1.8     1.5 4  0.193

## P. cuvieri

``` r
stage_cuvieri_no_effect <- glmmTMB(Gosner.Stage~ 1 + (1|Aquarium), data=data_table_cuvieri, family=gaussian(link = "identity"))
stage_cuvieri_land <- glmmTMB(Gosner.Stage~ Treatment + (1|Aquarium), data=data_table_cuvieri, family=gaussian(link = "identity"))
stage_cuvieri_Pasture_Sugarcane_model <- glmmTMB(Gosner.Stage ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=data_table_cuvieri, family=gaussian(link = "identity"))
stage_cuvieri_Pasture_Cerrado_model <- glmmTMB(Gosner.Stage~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table_cuvieri, family=gaussian(link = "identity"))
stage_cuvieri_Sugarcane_Cerrado_model <- glmmTMB(Gosner.Stage~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table_cuvieri, family=gaussian(link = "identity"))


AICTAB_stage_cuvieri <- AICctab(stage_cuvieri_no_effect, stage_cuvieri_land, stage_cuvieri_Pasture_Sugarcane_model,stage_cuvieri_Pasture_Cerrado_model,stage_cuvieri_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                     mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_stage_cuvieri
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -339.4  684.9    0.0    20.1 3  <0.001
    ## Savanna # Pasture # Sugarcane -327.3  664.8   12.1     0.0 5  0.9442
    ## Savanna # Pasture = Sugarcane -334.9  678.0    4.5    13.2 4  0.0013
    ## Savanna = Pasture # Sugarcane -339.3  686.8    0.1    22.0 4  <0.001
    ## Savanna = Sugarcane # Pasture -331.2  670.5    8.2     5.7 4  0.0545

## P. nattereri

``` r
stage_nattereri_no_effect <- glmmTMB(Gosner.Stage~ 1 + (1|Aquarium), data=data_table_nattereri, family=gaussian(link = "identity"))
stage_nattereri_land <- glmmTMB(Gosner.Stage~ Treatment + (1|Aquarium), data=data_table_nattereri, family=gaussian(link = "identity"))
stage_nattereri_Pasture_Sugarcane_model <- glmmTMB(Gosner.Stage ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=data_table_nattereri, family=gaussian(link = "identity"))
stage_nattereri_Pasture_Cerrado_model <- glmmTMB(Gosner.Stage~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table_nattereri, family=gaussian(link = "identity"))
stage_nattereri_Sugarcane_Cerrado_model <- glmmTMB(Gosner.Stage~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table_nattereri, family=gaussian(link = "identity"))


AICTAB_stage_nattereri <- AICctab(stage_nattereri_no_effect, stage_nattereri_land, stage_nattereri_Pasture_Sugarcane_model,stage_nattereri_Pasture_Cerrado_model,stage_nattereri_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                       mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_stage_nattereri
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -379.7  765.5    0.0     7.4 3  0.0148
    ## Savanna # Pasture # Sugarcane -374.4  759.2    5.3     1.1 5  0.3463
    ## Savanna # Pasture = Sugarcane -374.9  758.1    4.7     0.0 4  0.5936
    ## Savanna = Pasture # Sugarcane -379.4  767.1    0.3     8.9 4  0.0068
    ## Savanna = Sugarcane # Pasture -377.7  763.6    2.0     5.5 4  0.0385

## P. centralis

``` r
stage_centralis_no_effect <- glmmTMB(Gosner.Stage~ 1 + (1|Aquarium), data=data_table_centralis, family=gaussian(link = "identity"))
stage_centralis_land <- glmmTMB(Gosner.Stage~ Treatment + (1|Aquarium), data=data_table_centralis, family=gaussian(link = "identity"))
stage_centralis_Pasture_Sugarcane_model <- glmmTMB(Gosner.Stage ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=data_table_centralis, family=gaussian(link = "identity"))
stage_centralis_Pasture_Cerrado_model <- glmmTMB(Gosner.Stage~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table_centralis, family=gaussian(link = "identity"))
stage_centralis_Sugarcane_Cerrado_model <- glmmTMB(Gosner.Stage~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table_centralis, family=gaussian(link = "identity"))


AICTAB_stage_centralis <- AICctab(stage_centralis_no_effect, stage_centralis_land, stage_centralis_Pasture_Sugarcane_model,stage_centralis_Pasture_Cerrado_model,stage_centralis_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                       mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_stage_centralis 
```

    ##                               logLik AICc  dLogLik dAICc df weight
    ## No Effect                     -93.8  194.3   0.0     3.5 3  0.103 
    ## Savanna # Pasture # Sugarcane -90.9  193.3   3.0     2.5 5  0.170 
    ## Savanna # Pasture = Sugarcane -93.0  195.1   0.8     4.3 4  0.069 
    ## Savanna = Pasture # Sugarcane -90.9  190.7   3.0     0.0 4  0.605 
    ## Savanna = Sugarcane # Pasture -93.3  195.6   0.5     4.9 4  0.052

# BCI

``` r
bci_table <- data_table[-which(is.na(data_table$Gosner.Stage)), 
                        c("Tadpole.Mass","Gosner.Stage","Leaf.Litter.Treatment","Species", "Aquarium")]
bci_table$Treatment <- factor(bci_table$Leaf.Litter.Treatment, levels = c("Cerrado","Pasture","Sugarcane"))
bci_table$Aquarium <- factor(bci_table$Aquarium)


quadratic <- lm(log(bci_table$Tadpole.Mass) ~ bci_table$Gosner.Stage + 
                  I(bci_table$Gosner.Stage^2))
bci_table$BCI <- residuals(quadratic)

Treatment_Pasture_Sugarcane <- rep(NA, nrow(bci_table))
for(i in 1:nrow(bci_table)){
  if(bci_table$Treatment[i] == "Pasture" | bci_table$Treatment[i] == "Sugarcane"){
    Treatment_Pasture_Sugarcane[i] <- "Pasture-Sugarcane"
  }else{Treatment_Pasture_Sugarcane[i] <- as.character(bci_table$Treatment[i])}
}

Treatment_Pasture_Cerrado <- rep(NA, nrow(bci_table))
for(i in 1:nrow(bci_table)){
  if(bci_table$Treatment[i] == "Pasture" | bci_table$Treatment[i] == "Cerrado"){
    Treatment_Pasture_Cerrado[i] <- "Pasture-Cerrado"
  }else{Treatment_Pasture_Cerrado[i] <- as.character(bci_table$Treatment[i])}
}

Treatment_Sugarcane_Cerrado <- rep(NA, nrow(bci_table))
for(i in 1:nrow(bci_table)){
  if(bci_table$Treatment[i] == "Sugarcane" | bci_table$Treatment[i] == "Cerrado"){
    Treatment_Sugarcane_Cerrado[i] <- "Sugarcane-Cerrado"
  }else{Treatment_Sugarcane_Cerrado[i] <- as.character(bci_table$Treatment[i])}
}


bci_table <- data.frame(bci_table, Treatment_Pasture_Sugarcane, Treatment_Pasture_Cerrado,Treatment_Sugarcane_Cerrado)


bci_table_nattereri <- bci_table[which(bci_table$Species == "P. nattereri"),]
bci_table_scinax <- bci_table[which(bci_table$Species == "S. fuscovarius"),]
bci_table_cuvieri <- bci_table[which(bci_table$Species == "P. cuvieri"),]
bci_table_centralis <- bci_table[which(bci_table$Species == "P. centralis"),]

bci_table <- rbind(bci_table_centralis ,bci_table_cuvieri, bci_table_nattereri, bci_table_scinax)
```

## S. fuscovarius

``` r
bci_scinax_no_effect <- glmmTMB(BCI~ 1 + (1|Aquarium), data=bci_table_scinax, family=gaussian(link = "identity"))
bci_scinax_land <- glmmTMB(BCI~ Treatment + (1|Aquarium), data=bci_table_scinax, family=gaussian(link = "identity"))
bci_scinax_Pasture_Sugarcane_model <- glmmTMB(BCI ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=bci_table_scinax, family=gaussian(link = "identity"))
bci_scinax_Pasture_Cerrado_model <- glmmTMB(BCI~ Treatment_Pasture_Cerrado + (1|Aquarium), data=bci_table_scinax, family=gaussian(link = "identity"))
bci_scinax_Sugarcane_Cerrado_model <- glmmTMB(BCI~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=bci_table_scinax, family=gaussian(link = "identity"))


AICTAB_bci_scinax <- AICctab(bci_scinax_no_effect, bci_scinax_land, bci_scinax_Pasture_Sugarcane_model,bci_scinax_Pasture_Cerrado_model,bci_scinax_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                               mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_bci_scinax
```

    ##                               logLik AICc dLogLik dAICc df weight
    ## No Effect                     0.4    5.4  0.0     3.6   3  0.088 
    ## Savanna # Pasture # Sugarcane 3.4    3.5  3.0     1.8   5  0.226 
    ## Savanna # Pasture = Sugarcane 1.5    5.2  1.2     3.4   4  0.099 
    ## Savanna = Pasture # Sugarcane 0.6    7.0  0.2     5.2   4  0.040 
    ## Savanna = Sugarcane # Pasture 3.2    1.7  2.9     0.0   4  0.547

## P. cuvieri

``` r
bci_cuvieri_no_effect <- glmmTMB(BCI~ 1 + (1|Aquarium), data=bci_table_cuvieri, family=gaussian(link = "identity"))
bci_cuvieri_land <- glmmTMB(BCI~ Treatment + (1|Aquarium), data=bci_table_cuvieri, family=gaussian(link = "identity"))
bci_cuvieri_Pasture_Sugarcane_model <- glmmTMB(BCI ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=bci_table_cuvieri, family=gaussian(link = "identity"))
bci_cuvieri_Pasture_Cerrado_model <- glmmTMB(BCI~ Treatment_Pasture_Cerrado + (1|Aquarium), data=bci_table_cuvieri, family=gaussian(link = "identity"))
bci_cuvieri_Sugarcane_Cerrado_model <- glmmTMB(BCI~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=bci_table_cuvieri, family=gaussian(link = "identity"))


AICTAB_bci_cuvieri <- AICctab(bci_cuvieri_no_effect, bci_cuvieri_land, bci_cuvieri_Pasture_Sugarcane_model,bci_cuvieri_Pasture_Cerrado_model,bci_cuvieri_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_bci_cuvieri
```

    ##                               logLik AICc  dLogLik dAICc df weight
    ## No Effect                      17.5  -28.9   0.0     0.0 3  0.30  
    ## Savanna # Pasture # Sugarcane  18.6  -27.0   1.1     1.9 5  0.11  
    ## Savanna # Pasture = Sugarcane  18.3  -28.3   0.7     0.6 4  0.22  
    ## Savanna = Pasture # Sugarcane  17.5  -26.8   0.0     2.1 4  0.11  
    ## Savanna = Sugarcane # Pasture  18.4  -28.6   0.9     0.3 4  0.26

## P. nattereri

``` r
bci_nattereri_no_effect <- glmmTMB(BCI~ 1 + (1|Aquarium), data=bci_table_nattereri, family=gaussian(link = "identity"))
bci_nattereri_land <- glmmTMB(BCI~ Treatment + (1|Aquarium), data=bci_table_nattereri, family=gaussian(link = "identity"))
bci_nattereri_Pasture_Sugarcane_model <- glmmTMB(BCI ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=bci_table_nattereri, family=gaussian(link = "identity"))
bci_nattereri_Pasture_Cerrado_model <- glmmTMB(BCI~ Treatment_Pasture_Cerrado + (1|Aquarium), data=bci_table_nattereri, family=gaussian(link = "identity"))
bci_nattereri_Sugarcane_Cerrado_model <- glmmTMB(BCI~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=bci_table_nattereri, family=gaussian(link = "identity"))


AICTAB_bci_nattereri <- AICctab(bci_nattereri_no_effect, bci_nattereri_land, bci_nattereri_Pasture_Sugarcane_model,bci_nattereri_Pasture_Cerrado_model,bci_nattereri_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                  mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_bci_nattereri
```

    ##                               logLik AICc  dLogLik dAICc df weight
    ## No Effect                     -22.1   50.3   0.0     2.1 3  0.159 
    ## Savanna # Pasture # Sugarcane -19.8   50.0   2.3     1.8 5  0.182 
    ## Savanna # Pasture = Sugarcane -21.9   52.1   0.1     3.9 4  0.064 
    ## Savanna = Pasture # Sugarcane -21.1   50.5   1.0     2.3 4  0.145 
    ## Savanna = Sugarcane # Pasture -20.0   48.2   2.1     0.0 4  0.451

## P. centralis

``` r
bci_centralis_no_effect <- glmmTMB(BCI~ 1 + (1|Aquarium), data=bci_table_centralis, family=gaussian(link = "identity"))
bci_centralis_land <- glmmTMB(BCI~ Treatment + (1|Aquarium), data=bci_table_centralis, family=gaussian(link = "identity"))
bci_centralis_Pasture_Sugarcane_model <- glmmTMB(BCI ~ Treatment_Pasture_Sugarcane + (1|Aquarium), data=bci_table_centralis, family=gaussian(link = "identity"))
bci_centralis_Pasture_Cerrado_model <- glmmTMB(BCI~ Treatment_Pasture_Cerrado + (1|Aquarium), data=bci_table_centralis, family=gaussian(link = "identity"))
bci_centralis_Sugarcane_Cerrado_model <- glmmTMB(BCI~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=bci_table_centralis, family=gaussian(link = "identity"))


AICTAB_bci_centralis <- AICctab(bci_centralis_no_effect, bci_centralis_land, bci_centralis_Pasture_Sugarcane_model,bci_centralis_Pasture_Cerrado_model,bci_centralis_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                  mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_bci_centralis 
```

    ##                               logLik AICc  dLogLik dAICc df weight
    ## No Effect                     -11.9   30.4   0.0     5.8 3  0.0291
    ## Savanna # Pasture # Sugarcane  -7.0   25.5   4.9     0.9 5  0.3367
    ## Savanna # Pasture = Sugarcane  -7.8   24.7   4.1     0.0 4  0.5179
    ## Savanna = Pasture # Sugarcane  -9.4   27.8   2.5     3.2 4  0.1069
    ## Savanna = Sugarcane # Pasture -11.8   32.7   0.1     8.0 4  0.0094
