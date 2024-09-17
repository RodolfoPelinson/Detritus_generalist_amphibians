library(glmmTMB)
library(bbmle)
library(DHARMa)


data_table <- read.csv("Water_Quality_data.csv" )

dir()
data_table$Treatment <- factor(data_table$Leaf.Litter.Treatment, levels = c("Cerrado", "Pasture", "Sugarcane"))


Treatment_Pasture_Sugarcane <- rep(NA, nrow(data_table))
for(i in 1:nrow(data_table)){
  if(data_table$Treatment[i] == "Pasture" | data_table$Treatment[i] == "Sugarcane"){
    Treatment_Pasture_Sugarcane[i] <- "Pasture-Sugarcane"
  }else{Treatment_Pasture_Sugarcane[i] <- as.character(data_table$Treatment[i])}
}
Treatment_Pasture_Sugarcane

Treatment_Pasture_Cerrado <- rep(NA, nrow(data_table))
for(i in 1:nrow(data_table)){
  if(data_table$Treatment[i] == "Pasture" | data_table$Treatment[i] == "Cerrado"){
    Treatment_Pasture_Cerrado[i] <- "Pasture-Cerrado"
  }else{Treatment_Pasture_Cerrado[i] <- as.character(data_table$Treatment[i])}
}
Treatment_Pasture_Cerrado

Treatment_Sugarcane_Cerrado <- rep(NA, nrow(data_table))
for(i in 1:nrow(data_table)){
  if(data_table$Treatment[i] == "Sugarcane" | data_table$Treatment[i] == "Cerrado"){
    Treatment_Sugarcane_Cerrado[i] <- "Sugarcane-Cerrado"
  }else{Treatment_Sugarcane_Cerrado[i] <- as.character(data_table$Treatment[i])}
}
Treatment_Sugarcane_Cerrado

data_table <- data.frame(data_table, Treatment_Pasture_Sugarcane, Treatment_Pasture_Cerrado,Treatment_Sugarcane_Cerrado)


##


Conductivity_no_effect <- glmmTMB(Conductivity~ 1 + (1|Aquarium), family=gaussian(link = "identity"),  data=data_table)
Conductivity_land <- glmmTMB(Conductivity~ Treatment + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
Conductivity_Date <-glmmTMB(Conductivity~ Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
Conductivity_land_Date <-glmmTMB(Conductivity~ Treatment + Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
Conductivity_land_Date_interaction <-glmmTMB(Conductivity~ Treatment * Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))


AICTAB_Conductivity_int <- AICctab(Conductivity_no_effect, Conductivity_land, Conductivity_Date, Conductivity_land_Date,Conductivity_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                                    mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_Conductivity_int

Conductivity_Pasture_Sugarcane_model <- glmmTMB(Conductivity~ Treatment_Pasture_Sugarcane *Date + (1|Aquarium), family=gaussian(link = "identity"), data=data_table)
Conductivity_Pasture_Cerrado_model <- glmmTMB(Conductivity~ Treatment_Pasture_Cerrado* Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
Conductivity_Sugarcane_Cerrado_model <- glmmTMB(Conductivity~ Treatment_Sugarcane_Cerrado*Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))


AICTAB_Conductivity <- AICctab(Conductivity_no_effect, Conductivity_land_Date_interaction, Conductivity_Pasture_Sugarcane_model,Conductivity_Pasture_Cerrado_model,Conductivity_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                                mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_Conductivity
###

Temperature_no_effect <- glmmTMB(Temperature~ 1 + (1|Aquarium), family=Gamma(link = "log"),  data=data_table)
Temperature_land <- glmmTMB(Temperature~ Treatment + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Temperature_Date <-glmmTMB(Temperature~ Date + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Temperature_land_Date <-glmmTMB(Temperature~ Treatment + Date+ (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Temperature_land_Date_interaction <-glmmTMB(Temperature~ Treatment * Date+ (1|Aquarium), data=data_table, family=Gamma(link = "log"))

AICTAB_Temperature_int <- AICctab(Temperature_no_effect, Temperature_land, Temperature_Date, Temperature_land_Date,Temperature_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                                  mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_Temperature_int

Temperature_Pasture_Sugarcane_model <- glmmTMB(Temperature~ Treatment_Pasture_Sugarcane + (1|Aquarium), family=Gamma(link = "log"), data=data_table)
Temperature_Pasture_Cerrado_model <- glmmTMB(Temperature~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Temperature_Sugarcane_Cerrado_model <- glmmTMB(Temperature~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table, family=Gamma(link = "log"))


AICTAB_Temperature <- AICctab(Temperature_no_effect, Temperature_land, Temperature_Pasture_Sugarcane_model,Temperature_Pasture_Cerrado_model,Temperature_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                              mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_Temperature

####


pH_no_effect <- glmmTMB(pH~ 1 + (1|Aquarium), family=gaussian(link = "identity"),  data=data_table)
pH_land <- glmmTMB(pH~ Treatment + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
pH_Date <-glmmTMB(pH~ Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
pH_land_Date <-glmmTMB(pH~ Treatment + Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
pH_land_Date_interaction <-glmmTMB(pH~ Treatment * Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))


AICTAB_pH_int <- AICctab(pH_no_effect, pH_land, pH_Date, pH_land_Date,pH_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                         mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_pH_int


pH_Pasture_Sugarcane_model <- glmmTMB(pH~ Treatment_Pasture_Sugarcane +Date+ (1|Aquarium), family=gaussian(link = "identity"), dispformula = ~Treatment, data=data_table)
pH_Pasture_Cerrado_model <- glmmTMB(pH~ Treatment_Pasture_Cerrado+Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"), dispformula = ~Treatment)
pH_Sugarcane_Cerrado_model <- glmmTMB(pH~ Treatment_Sugarcane_Cerrado+Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"), dispformula = ~Treatment)


AICTAB_pH <- AICctab(pH_no_effect, pH_land_Date, pH_Pasture_Sugarcane_model,pH_Pasture_Cerrado_model,pH_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                     mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_pH
####

Turbidity_no_effect <- glmmTMB(Turbidity~ 1 + (1|Aquarium), family=Gamma(link = "log"),  data=data_table)
Turbidity_land <- glmmTMB(Turbidity~ Treatment + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Turbidity_Date <-glmmTMB(Turbidity~ Date + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Turbidity_land_Date <-glmmTMB(Turbidity~ Treatment + Date+ (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Turbidity_land_Date_interaction <-glmmTMB(Turbidity~ Treatment * Date+ (1|Aquarium), data=data_table, family=Gamma(link = "log"))


AICTAB_Turbidity_int <- AICctab(Turbidity_no_effect, Turbidity_land, Turbidity_Date, Turbidity_land_Date,Turbidity_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                               mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_Turbidity_int


Turbidity_Pasture_Sugarcane_model <- glmmTMB(Turbidity~ Treatment_Pasture_Sugarcane+Date + (1|Aquarium), family=Gamma(link = "log"), data=data_table)
Turbidity_Pasture_Cerrado_model <- glmmTMB(Turbidity~ Treatment_Pasture_Cerrado+Date + (1|Aquarium), data=data_table, family=Gamma(link = "log"))
Turbidity_Sugarcane_Cerrado_model <- glmmTMB(Turbidity~ Treatment_Sugarcane_Cerrado+Date + (1|Aquarium), data=data_table, family=Gamma(link = "log"))


AICTAB_Turbidity <- AICctab(Turbidity_no_effect, Turbidity_land_Date, Turbidity_Pasture_Sugarcane_model,Turbidity_Pasture_Cerrado_model,Turbidity_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                           mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_Turbidity
###

DOp_no_effect <- glmmTMB(DOp~ 1 + (1|Aquarium), family=gaussian(link = "identity"),  data=data_table)
DOp_land <- glmmTMB(DOp~ Treatment + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOp_Date <-glmmTMB(DOp~ Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOp_land_Date <-glmmTMB(DOp~ Treatment + Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOp_land_Date_interaction <-glmmTMB(DOp~ Treatment * Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))


AICTAB_DOp_int <- AICctab(DOp_no_effect, DOp_land, DOp_Date, DOp_land_Date,DOp_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                          mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_DOp_int



DOp_Pasture_Sugarcane_model <- glmmTMB(DOp~ Treatment_Pasture_Sugarcane + (1|Aquarium), family=gaussian (link = "identity"), data=data_table)
DOp_Pasture_Cerrado_model <- glmmTMB(DOp~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table, family=gaussian (link = "identity"))
DOp_Sugarcane_Cerrado_model <- glmmTMB(DOp~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table, family=gaussian (link = "identity"))


AICTAB_DOp <- AICctab(DOp_no_effect, DOp_land, DOp_Pasture_Sugarcane_model,DOp_Pasture_Cerrado_model,DOp_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                      mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))
AICTAB_DOp
####

DOmg_no_effect <- glmmTMB(DOmg~ 1 + (1|Aquarium), family=gaussian(link = "identity"),  data=data_table)
DOmg_land <- glmmTMB(DOmg~ Treatment + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOmg_Date <-glmmTMB(DOmg~ Date + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOmg_land_Date <-glmmTMB(DOmg~ Treatment + Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOmg_land_Date_interaction <-glmmTMB(DOmg~ Treatment * Date+ (1|Aquarium), data=data_table, family=gaussian(link = "identity"))


AICTAB_DOmg_int <- AICctab(DOmg_no_effect, DOmg_land, DOmg_Date, DOmg_land_Date,DOmg_land_Date_interaction ,weights = T, base = T, logLik = T, sort = F,
                           mnames = c("No Effect", "Land Use", "Time Effect", "Land Use and Time Effect", "Land Use and Time Effect Interaction" ))
AICTAB_DOmg_int


DOmg_Pasture_Sugarcane_model <- glmmTMB(DOmg~ Treatment_Pasture_Sugarcane + (1|Aquarium), family=gaussian(link = "identity"), data=data_table)
DOmg_Pasture_Cerrado_model <- glmmTMB(DOmg~ Treatment_Pasture_Cerrado + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))
DOmg_Sugarcane_Cerrado_model <- glmmTMB(DOmg~ Treatment_Sugarcane_Cerrado + (1|Aquarium), data=data_table, family=gaussian(link = "identity"))

AICTAB_DOmg <- AICctab(DOmg_no_effect, DOmg_land, DOmg_Pasture_Sugarcane_model,DOmg_Pasture_Cerrado_model,DOmg_Sugarcane_Cerrado_model,weights = T, base = T, logLik = T, sort = F,
                       mnames = c("No Effect", "Savanna # Pasture # Sugarcane", "Savanna # Pasture = Sugarcane", "Savanna = Pasture # Sugarcane", "Savanna = Sugarcane # Pasture" ))

AICTAB_DOmg


