Experimental study (Water)
================
Rodolfo Pelinson
2024-09-17

# Water and nutrients as a function of different detritus type.

``` r
library(glmmTMB)
library(bbmle)
library(DHARMa)
```

Loading data

## Water Parameters

``` r
path <- c("C:/Users/rodol/OneDrive/repos/Detritus_generalist_amphibians/")

data_table <- read.csv(paste(path,"data/Water_Quality_data.csv", sep = "") )

data_table$Treatment <- factor(data_table$Leaf.Litter.Treatment, levels = c("Cerrado", "Pasture", "Sugarcane"))


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
```

### Conductivity

``` r
Conductivity_no_effect <- glmmTMB(Conductivity~ 1 + (1|Aquarium), family=gaussian(link = "identity"),  data=data_table)
Conductivity_land <- glmmTMB(Conductivity~ Treatment + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
Conductivity_Date <-glmmTMB(Conductivity~ Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
Conductivity_land_Date <-glmmTMB(Conductivity~ Treatment + Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
Conductivity_land_Date_interaction <-glmmTMB(Conductivity~ Treatment * Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))


AICTAB_Conductivity_int <- AICctab(Conductivity_no_effect, Conductivity_land, Conductivity_Date, Conductivity_land_Date,Conductivity_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                                    mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_Conductivity_int
```

    ##                                      logLik AICc   dLogLik dAICc  df weight
    ## No Effect                            -741.2 1488.6    0.0   175.1 3  <0.001
    ## Land Use                             -668.6 1347.7   72.6    34.2 5  <0.001
    ## Time Effect                          -731.2 1472.9   10.0   159.4 5  <0.001
    ## Land Use and Time Effect             -658.6 1332.1   82.6    18.6 7  <0.001
    ## Land Use and Time Effect Interaction -644.7 1313.5   96.5     0.0 11 1

Post-hoc comparisons

``` r
Conductivity_Pasture_Sugarcane_model <- glmmTMB(Conductivity~ Treatment_Pasture_Sugarcane *Date + (1|Aquarium), family=gaussian(link = "identity"), data=data_table)
Conductivity_Pasture_Cerrado_model <- glmmTMB(Conductivity~ Treatment_Pasture_Cerrado* Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
Conductivity_Sugarcane_Cerrado_model <- glmmTMB(Conductivity~ Treatment_Sugarcane_Cerrado*Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))


AICTAB_Conductivity <- AICctab(Conductivity_no_effect, Conductivity_land_Date_interaction, Conductivity_Pasture_Sugarcane_model,Conductivity_Pasture_Cerrado_model,Conductivity_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_Conductivity
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -741.2 1488.6    0.0   175.1 3  <0.001
    ## Savanna # Pasture # Sugarcane -644.7 1313.5   96.5     0.0 11 1     
    ## Savanna # Pasture = Sugarcane -702.7 1422.5   38.5   109.1 8  <0.001
    ## Savanna = Pasture # Sugarcane -728.1 1473.3   13.1   159.8 8  <0.001
    ## Savanna = Sugarcane # Pasture -674.0 1365.1   67.2    51.6 8  <0.001

### Temperature

``` r
Temperature_no_effect <- glmmTMB(Temperature~ 1 + (1|Aquarium), family=Gamma(link = "log"),  data=data_table)
Temperature_land <- glmmTMB(Temperature~ Treatment + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Temperature_Date <-glmmTMB(Temperature~ Date + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Temperature_land_Date <-glmmTMB(Temperature~ Treatment + Date+ (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Temperature_land_Date_interaction <-glmmTMB(Temperature~ Treatment * Date+ (1|Aquarium), data=data_table, family=Gamma(link = "log"))

AICTAB_Temperature_int <- AICctab(Temperature_no_effect, Temperature_land, Temperature_Date, Temperature_land_Date,Temperature_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                                  mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_Temperature_int
```

    ##                                      logLik AICc   dLogLik dAICc  df weight
    ## No Effect                            -223.5  453.2    0.0   113.6 3  <0.001
    ## Land Use                             -222.2  454.8    1.4   115.2 5  <0.001
    ## Time Effect                          -164.6  339.6   59.0     0.0 5  0.60  
    ## Land Use and Time Effect             -163.2  341.2   60.3     1.6 7  0.26  
    ## Land Use and Time Effect Interaction -159.3  342.5   64.3     2.9 11 0.14

Post-hoc comparisons

``` r
Temperature_Pasture_Sugarcane_model <- glmmTMB(Temperature~ Treatment_Pasture_Sugarcane + (1|Aquarium), family=Gamma(link = "log"), data=data_table)
Temperature_Pasture_Cerrado_model <- glmmTMB(Temperature~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Temperature_Sugarcane_Cerrado_model <- glmmTMB(Temperature~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table, family=Gamma(link = "log"))


AICTAB_Temperature <- AICctab(Temperature_no_effect, Temperature_land, Temperature_Pasture_Sugarcane_model,Temperature_Pasture_Cerrado_model,Temperature_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                              mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_Temperature
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -223.5  453.2    0.0     0.4 3  0.266 
    ## Savanna # Pasture # Sugarcane -222.2  454.8    1.4     1.9 5  0.124 
    ## Savanna # Pasture = Sugarcane -223.5  455.2    0.1     2.4 4  0.098 
    ## Savanna = Pasture # Sugarcane -222.3  452.9    1.2     0.0 4  0.320 
    ## Savanna = Sugarcane # Pasture -222.8  453.9    0.7     1.0 4  0.192

### pH

``` r
pH_no_effect <- glmmTMB(pH~ 1 + (1|Aquarium), family=gaussian(link = "identity"),  data=data_table)
pH_land <- glmmTMB(pH~ Treatment + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
pH_Date <-glmmTMB(pH~ Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
pH_land_Date <-glmmTMB(pH~ Treatment + Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
pH_land_Date_interaction <-glmmTMB(pH~ Treatment * Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))


AICTAB_pH_int <- AICctab(pH_no_effect, pH_land, pH_Date, pH_land_Date,pH_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                         mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_pH_int
```

    ##                                      logLik AICc  dLogLik dAICc df weight
    ## No Effect                            -43.5   93.2   0.0    13.8 3  <0.001
    ## Land Use                             -36.3   83.0   7.3     3.5 5  0.1438
    ## Time Effect                          -39.6   89.6   3.9    10.1 5  0.0052
    ## Land Use and Time Effect             -32.3   79.5  11.2     0.0 7  0.8254
    ## Land Use and Time Effect Interaction -31.2   86.5  12.3     7.0 11 0.0247

Post-hoc comparisons

``` r
pH_Pasture_Sugarcane_model <- glmmTMB(pH~ Treatment_Pasture_Sugarcane +Date+ (1|Aquarium), family=gaussian(link = "identity"), dispformula = ~Treatment, data=data_table)
pH_Pasture_Cerrado_model <- glmmTMB(pH~ Treatment_Pasture_Cerrado+Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"), dispformula = ~Treatment)
pH_Sugarcane_Cerrado_model <- glmmTMB(pH~ Treatment_Sugarcane_Cerrado+Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"), dispformula = ~Treatment)


AICTAB_pH <- AICctab(pH_no_effect, pH_land_Date, pH_Pasture_Sugarcane_model,pH_Pasture_Cerrado_model,pH_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                     mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_pH
```

    ##                               logLik AICc  dLogLik dAICc df weight
    ## No Effect                     -43.5   93.2   0.0    73.5 3  <0.001
    ## Savanna # Pasture # Sugarcane -32.3   79.5  11.2    59.8 7  <0.001
    ## Savanna # Pasture = Sugarcane  -1.3   19.7  42.2     0.0 8  0.9871
    ## Savanna = Pasture # Sugarcane  -6.1   29.2  37.5     9.5 8  0.0085
    ## Savanna = Sugarcane # Pasture  -6.7   30.5  36.8    10.8 8  0.0044

### Turbidity

``` r
Turbidity_no_effect <- glmmTMB(Turbidity~ 1 + (1|Aquarium), family=Gamma(link = "log"),  data=data_table)
Turbidity_land <- glmmTMB(Turbidity~ Treatment + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Turbidity_Date <-glmmTMB(Turbidity~ Date + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Turbidity_land_Date <-glmmTMB(Turbidity~ Treatment + Date+ (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Turbidity_land_Date_interaction <-glmmTMB(Turbidity~ Treatment * Date+ (1|Aquarium), data=data_table, family=Gamma(link = "log"))


AICTAB_Turbidity_int <- AICctab(Turbidity_no_effect, Turbidity_land, Turbidity_Date, Turbidity_land_Date,Turbidity_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                               mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_Turbidity_int
```

    ##                                      logLik AICc   dLogLik dAICc  df weight
    ## No Effect                            -213.9  434.0    0.0    42.2 3  <0.001
    ## Land Use                             -198.1  406.6   15.8    14.9 5  <0.001
    ## Time Effect                          -204.1  418.7    9.8    27.0 5  <0.001
    ## Land Use and Time Effect             -188.5  391.8   25.4     0.0 7  0.66  
    ## Land Use and Time Effect Interaction -184.5  393.1   29.4     1.3 11 0.34

Post-hoc comparisons

``` r
Turbidity_Pasture_Sugarcane_model <- glmmTMB(Turbidity~ Treatment_Pasture_Sugarcane+Date + (1|Aquarium), family=Gamma(link = "log"), data=data_table)
Turbidity_Pasture_Cerrado_model <- glmmTMB(Turbidity~ Treatment_Pasture_Cerrado+Date + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Turbidity_Sugarcane_Cerrado_model <- glmmTMB(Turbidity~ Treatment_Sugarcane_Cerrado+Date + (1|Aquarium), data=data_table, family=Gamma(link = "log"))


AICTAB_Turbidity <- AICctab(Turbidity_no_effect, Turbidity_land_Date, Turbidity_Pasture_Sugarcane_model,Turbidity_Pasture_Cerrado_model,Turbidity_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                           mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_Turbidity
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -213.9  434.0    0.0    43.8 3  <0.001
    ## Savanna # Pasture # Sugarcane -188.5  391.8   25.4     1.6 7  0.31  
    ## Savanna # Pasture = Sugarcane -202.3  417.2   11.6    27.0 6  <0.001
    ## Savanna = Pasture # Sugarcane -199.5  411.5   14.4    21.4 6  <0.001
    ## Savanna = Sugarcane # Pasture -188.8  390.2   25.1     0.0 6  0.69

### DOp

``` r
DOp_no_effect <- glmmTMB(DOp~ 1 + (1|Aquarium), family=gaussian(link = "identity"),  data=data_table)
DOp_land <- glmmTMB(DOp~ Treatment + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOp_Date <-glmmTMB(DOp~ Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOp_land_Date <-glmmTMB(DOp~ Treatment + Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOp_land_Date_interaction <-glmmTMB(DOp~ Treatment * Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))


AICTAB_DOp_int <- AICctab(DOp_no_effect, DOp_land, DOp_Date, DOp_land_Date,DOp_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                          mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_DOp_int
```

    ##                                      logLik AICc   dLogLik dAICc  df weight
    ## No Effect                            -587.5 1181.2    0.0     0.0 3  0.428 
    ## Land Use                             -586.5 1183.5    1.0     2.3 5  0.137 
    ## Time Effect                          -585.7 1181.8    1.8     0.6 5  0.314 
    ## Land Use and Time Effect             -584.7 1184.2    2.8     3.0 7  0.095 
    ## Land Use and Time Effect Interaction -581.4 1186.8    6.1     5.6 11 0.026

Post-hoc comparisons

``` r
DOp_Pasture_Sugarcane_model <- glmmTMB(DOp~ Treatment_Pasture_Sugarcane + (1|Aquarium), family=gaussian (link = "identity"), data=data_table)
DOp_Pasture_Cerrado_model <- glmmTMB(DOp~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table, family=gaussian (link = "identity"))
DOp_Sugarcane_Cerrado_model <- glmmTMB(DOp~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table, family=gaussian (link = "identity"))


AICTAB_DOp <- AICctab(DOp_no_effect, DOp_land, DOp_Pasture_Sugarcane_model,DOp_Pasture_Cerrado_model,DOp_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                      mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_DOp
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -587.5 1181.2    0.0     0.0 3  0.32  
    ## Savanna # Pasture # Sugarcane -586.5 1183.5    1.0     2.3 5  0.10  
    ## Savanna # Pasture = Sugarcane -586.6 1181.4    0.9     0.2 4  0.29  
    ## Savanna = Pasture # Sugarcane -587.4 1183.1    0.1     1.9 4  0.12  
    ## Savanna = Sugarcane # Pasture -587.1 1182.4    0.4     1.2 4  0.17

### DOmg

``` r
DOmg_no_effect <- glmmTMB(DOmg~ 1 + (1|Aquarium), family=gaussian(link = "identity"),  data=data_table)
DOmg_land <- glmmTMB(DOmg~ Treatment + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOmg_Date <-glmmTMB(DOmg~ Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOmg_land_Date <-glmmTMB(DOmg~ Treatment + Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOmg_land_Date_interaction <-glmmTMB(DOmg~ Treatment * Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))


AICTAB_DOmg_int <- AICctab(DOmg_no_effect, DOmg_land, DOmg_Date, DOmg_land_Date,DOmg_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                           mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_DOmg_int
```

    ##                                      logLik AICc   dLogLik dAICc  df weight
    ## No Effect                            -231.3  468.7    0.0     0.0 3  0.6909
    ## Land Use                             -230.5  471.5    0.8     2.7 5  0.1763
    ## Time Effect                          -231.0  472.5    0.2     3.8 5  0.1040
    ## Land Use and Time Effect             -230.3  475.4    1.0     6.6 7  0.0249
    ## Land Use and Time Effect Interaction -227.6  479.1    3.7    10.4 11 0.0039

Post-hoc comparisons

``` r
DOmg_Pasture_Sugarcane_model <- glmmTMB(DOmg~ Treatment_Pasture_Sugarcane + (1|Aquarium), family=gaussian(link = "identity"), data=data_table)
DOmg_Pasture_Cerrado_model <- glmmTMB(DOmg~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOmg_Sugarcane_Cerrado_model <- glmmTMB(DOmg~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))

AICTAB_DOmg <- AICctab(DOmg_no_effect, DOmg_land, DOmg_Pasture_Sugarcane_model,DOmg_Pasture_Cerrado_model,DOmg_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                       mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))

AICTAB_DOmg
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -231.3  468.7    0.0     0.0 3  0.35  
    ## Savanna # Pasture # Sugarcane -230.5  471.5    0.8     2.7 5  0.09  
    ## Savanna # Pasture = Sugarcane -230.6  469.4    0.7     0.7 4  0.25  
    ## Savanna = Pasture # Sugarcane -231.2  470.7    0.1     2.0 4  0.13  
    ## Savanna = Sugarcane # Pasture -230.9  470.2    0.4     1.4 4  0.17

## Nutrients

Loading data

``` r
path <- c("C:/Users/rodol/OneDrive/repos/Detritus_generalist_amphibians/")

#data_table <- read_delim(paste(path,"data/Water_Nutrients_data.csv", sep = ""), 
#                                   delim = "\t", escape_double = FALSE, 
#                                   col_types = cols(TOC = col_number(), 
#                                                    TC = col_number(), IC = col_number(), 
#                                                    TN = col_number(), TDP = col_number()), 
#                                   locale = locale(decimal_mark = ","), 
#                                   trim_ws = TRUE)

data_table <- read.csv(paste(path,"data/Water_Nutrients_data.csv", sep = ""))



data_table$Treatment <- factor(data_table$Leaf.Litter.Treatment, levels = c("Cerrado", "Pasture", "Sugarcane"))


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
```

### TOC

``` r
TOC_no_effect <- glmmTMB(TOC~ 1 , family=Gamma(link="log"),  data=data_table)
TOC_land <- glmmTMB(TOC~ Treatment  ,data=data_table, family=Gamma(link="log"))
TOC_Pasture_Sugarcane_model <- glmmTMB(TOC~ Treatment_Pasture_Sugarcane,  family=Gamma(link="log"), data=data_table)
TOC_Pasture_Cerrado_model <- glmmTMB(TOC~ Treatment_Pasture_Cerrado , data=data_table, family=Gamma(link="log"))
TOC_Sugarcane_Cerrado_model <- glmmTMB(TOC~ Treatment_Sugarcane_Cerrado , data=data_table, family=Gamma(link="log"))


AICTAB_TOC <- AICctab(TOC_no_effect, TOC_land, TOC_Pasture_Sugarcane_model,TOC_Pasture_Cerrado_model,TOC_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                      mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))

AICTAB_TOC
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -171.2  346.7    0.0    58.4 2  <0.001
    ## Savanna # Pasture # Sugarcane -140.4  289.8   30.8     1.4 4  0.33  
    ## Savanna # Pasture = Sugarcane -164.1  334.7    7.1    46.3 3  <0.001
    ## Savanna = Pasture # Sugarcane -167.3  341.2    3.9    52.8 3  <0.001
    ## Savanna = Sugarcane # Pasture -140.9  288.4   30.3     0.0 3  0.67

### TC

``` r
TC_no_effect <- glmmTMB(TC~ 1 , family=Gamma(link="log"),  data=data_table)
TC_land <- glmmTMB(TC~ Treatment  ,data=data_table, family=Gamma(link="log"))
TC_Pasture_Sugarcane_model <- glmmTMB(TC~ Treatment_Pasture_Sugarcane,  family=Gamma(link="log"), data=data_table)
TC_Pasture_Cerrado_model <- glmmTMB(TC~ Treatment_Pasture_Cerrado , data=data_table, family=Gamma(link="log"))
TC_Sugarcane_Cerrado_model <- glmmTMB(TC~ Treatment_Sugarcane_Cerrado , data=data_table, family=Gamma(link="log"))


AICTAB_TC <- AICctab(TC_no_effect, TC_land, TC_Pasture_Sugarcane_model,TC_Pasture_Cerrado_model,TC_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                     mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))

AICTAB_TC
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -189.6  383.4    0.0    83.5 2  <0.001
    ## Savanna # Pasture # Sugarcane -145.5  299.9   44.1     0.0 4  0.9978
    ## Savanna # Pasture = Sugarcane -175.7  357.9   13.9    58.0 3  <0.001
    ## Savanna = Pasture # Sugarcane -187.5  381.6    2.1    81.6 3  <0.001
    ## Savanna = Sugarcane # Pasture -152.8  312.2   36.8    12.3 3  0.0022

### IC

``` r
IC_no_effect <- glmmTMB(IC~ 1 , family=Gamma(link="log"),  data=data_table)
IC_land <- glmmTMB(IC~ Treatment  ,data=data_table, family=Gamma(link="log"))
IC_Pasture_Sugarcane_model <- glmmTMB(IC~ Treatment_Pasture_Sugarcane,  family=Gamma(link="log"), data=data_table)
IC_Pasture_Cerrado_model <- glmmTMB(IC~ Treatment_Pasture_Cerrado , data=data_table, family=Gamma(link="log"))
IC_Sugarcane_Cerrado_model <- glmmTMB(IC~ Treatment_Sugarcane_Cerrado , data=data_table, family=Gamma(link="log"))


AICTAB_IC <- AICctab(IC_no_effect, IC_land, IC_Pasture_Sugarcane_model,IC_Pasture_Cerrado_model,IC_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                     mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))

AICTAB_IC
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -143.7  291.6    0.0    87.8 2  <0.001
    ## Savanna # Pasture # Sugarcane  -97.4  203.8   46.2     0.0 4  1     
    ## Savanna # Pasture = Sugarcane -115.0  236.6   28.7    32.8 3  <0.001
    ## Savanna = Pasture # Sugarcane -143.6  293.7    0.1    89.9 3  <0.001
    ## Savanna = Sugarcane # Pasture -123.3  253.1   20.4    49.3 3  <0.001

### TN

``` r
TN_no_effect <- glmmTMB(TN~ 1 , family=Gamma(link="log"),  data=data_table)
TN_land <- glmmTMB(TN~ Treatment  ,data=data_table, family=Gamma(link="log"))
TN_Pasture_Sugarcane_model <- glmmTMB(TN~ Treatment_Pasture_Sugarcane,  family=Gamma(link="log"), data=data_table)
TN_Pasture_Cerrado_model <- glmmTMB(TN~ Treatment_Pasture_Cerrado , data=data_table, family=Gamma(link="log"))
TN_Sugarcane_Cerrado_model <- glmmTMB(TN~ Treatment_Sugarcane_Cerrado , data=data_table, family=Gamma(link="log"))


AICTAB_TN <- AICctab(TN_no_effect, TN_land, TN_Pasture_Sugarcane_model,TN_Pasture_Cerrado_model,TN_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                     mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))

AICTAB_TN
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -154.9  314.1    0.0    35.0 2  <0.001
    ## Savanna # Pasture # Sugarcane -135.1  279.1   19.8     0.0 4  0.77  
    ## Savanna # Pasture = Sugarcane -153.3  313.2    1.6    34.1 3  <0.001
    ## Savanna = Pasture # Sugarcane -145.7  298.0    9.2    18.9 3  <0.001
    ## Savanna = Sugarcane # Pasture -137.5  281.5   17.4     2.4 3  0.23

### TDP

``` r
TDP_no_effect <- glmmTMB(TDP~ 1 , family=gaussian(link = "identity"),  data=data_table)
TDP_land <- glmmTMB(TDP~ Treatment  ,data=data_table, family=gaussian(link = "identity"))
TDP_Pasture_Sugarcane_model <- glmmTMB(TDP~ Treatment_Pasture_Sugarcane,  family=gaussian(link = "identity"), data=data_table)
TDP_Pasture_Cerrado_model <- glmmTMB(TDP~ Treatment_Pasture_Cerrado , data=data_table, family=gaussian(link = "identity"))
TDP_Sugarcane_Cerrado_model <- glmmTMB(TDP~ Treatment_Sugarcane_Cerrado , data=data_table, family=gaussian(link = "identity"))


AICTAB_TDP <- AICctab(TDP_no_effect, TDP_land, TDP_Pasture_Sugarcane_model,TDP_Pasture_Cerrado_model,TDP_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                      mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))

AICTAB_TDP
```

    ##                               logLik AICc   dLogLik dAICc  df weight
    ## No Effect                     -381.7  767.6    0.0    51.6 2  <0.001
    ## Savanna # Pasture # Sugarcane -353.5  716.0   28.1     0.0 4  0.89  
    ## Savanna # Pasture = Sugarcane -356.8  720.1   24.9     4.1 3  0.11  
    ## Savanna = Pasture # Sugarcane -371.6  749.8   10.1    33.8 3  <0.001
    ## Savanna = Sugarcane # Pasture -380.5  767.5    1.2    51.6 3  <0.001
