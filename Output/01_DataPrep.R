rm(list = ls())

library(tidyverse) #always tidyverse
library(lubridate) #convert dates to readable form

dat <- read.csv("Data/Chapter2_Data_Clean.csv") #clean final data
geochem <- read.csv("Data/soil_pca.csv") #best fit
dat.ts <- read.csv("Data/Ch2.TimeSeries.Data .csv") #census series

#make the proper variables numeric

mod_data <- transform(dat, Init_Height = as.numeric(as.character(Init_Height)), 
                      Init_LeafNo = as.numeric(as.character(Init_LeafNo)), 
                      Init_Dia_1 = as.numeric(as.character(Init_Dia_1)), 
                      Init_Dia_2 = as.numeric(as.character(Init_Dia_2)), 
                      Final_Height = as.numeric(as.character(Final_Height)),
                      Root_LG = as.numeric(as.character(Root_LG)),
                      Final_LeafNo = as.numeric(as.character(Final_LeafNo)),
                      Final_Dia_1 = as.numeric(as.character(Final_Dia_1)),
                      Final_Dia_2 = as.numeric(as.character(Final_Dia_2)),
                      Stem_Mass = as.numeric(as.character(Stem_Mass)),
                      Root_Mass = as.numeric(as.character(Root_Mass)),
                      Tot_Biomass = as.numeric(as.character(Tot_Biomass))
)

#add harvest date column
mod_data$begin <- dmy(mod_data$begin)
mod_data$end<- dmy(mod_data$end)

exp.interval<- interval(mod_data$begin, mod_data$end) #calculate interval between dates

mod_data$duration <- time_length(exp.interval, unit = "day") #days within interval

#Combine Table and Position for random effect
mod_data <-  mutate( mod_data, Location = str_c( Table, Position, sep = "_" ))

#Fix misnamed Inoc_Sp

mod_data$Inoc_Sp <- ifelse(mod_data$Inoc_Tag == "19844",
                           "LACIAG", mod_data$Inoc_Sp)

#Create a column that identifies the control pots
mod_data <- mod_data %>% mutate(soil = case_when(Inoc_Sp_No == "C" ~ "control",
                                                 (Seedling == Inoc_Sp) & (Inoc_Sp_No != "C") ~ "conspecific",
                                                 Seedling != Inoc_Sp ~ "heterospecific"))

#Create a column that identifies a pot as CW, CD, HW, HD
mod_data<- mutate(mod_data, trt = str_c(soil, Treatment, sep = "_"))

#Create a "Status" column that denotes if a seedling is alive or dead at the end of the experiment
mod_data$status <- ifelse(!is.na(mod_data$Final_Height), "1", "0")

#Match Geochem and mod_data
#first clean up to make sure everything matches
geochem$Inoc_Sp<-gsub(" $","",geochem$Inoc_Sp,perl=T)

## Match tag names between geochem and perfdat where necessary
mod_data$Inoc_Tag <- ifelse(mod_data$Inoc_Tag == "", "22401", mod_data$Inoc_Tag)

mod_data$Inoc_Tag <- ifelse(mod_data$Inoc_Tag == "Psidfr", "xxx", mod_data$Inoc_Tag)

geochem <- geochem %>% mutate(Inoc_Tag = case_when(Inoc_Tag == "46965" ~ "59,46;1,1",
                                                   Inoc_Tag == "48187" ~ "57,43;4,3",
                                                   Inoc_Tag == "47905" ~ "57,55;3,3",
                                                   Inoc_Tag == "control" ~ "Control",
                                                   TRUE ~ Inoc_Tag))
# Group pots that shared the same soil inoculum
EUGENE1 <- c("C_15", "G_21", "K_31", "H_23", "C_24", "L_18", "E_29", "K_6", "L_24", "B_6")
EUGENE2 <- c("J_26", "I_24", "F_28", "J_31", "H_10", "J_32", "H_8", "J_3", "I_32", "H_25")
EUGENE3 <- c("L_19", "E_4", "F_16", "H_29", "A_18", "G_22", "I_31", "G_9", "G_25", "E_6")

GUAPST1 <- c("L_20", "F_18", "J_2", "E_10", "E_5", "I_8", "J_20", "L_12", "K_8", 
             "F_24", "F_2", "G_29", "G_12", "I_2")
GUAPST2 <- c("G_32", "H_27", "J_19", "F_8", "F_23", "G_16", "J_11", "K_26", "H_33", 
             "L_29", "F_26", "E_20", "H_22", "F_7")
GUAPST3 <- c("G_24", "L_2", "A_9", "K_24", "F_5", "B_26", "K_2", "F_11", "H_7",
             "I_17", "H_16", "H_15", "L_23", "G_8")

HEISCO1 <- c("E_25", "H_3", "G_33", "F_25", "J_14", "H_2", "G_13", "G_10")
HEISCO2 <- c("F_33", "I_25", "L_14", "J_6", "K_32", "K_3", "J_8", "I_19")
HEISCO3 <- c("J_4", "K_18", "H_26", "F_32", "L_5", "J_30", "H_17", "L_11")

# Assign groups to specific inoculum (pass 2)

mod_data$Inoc_Tag <- ifelse(mod_data$Location %in% EUGENE3, "20846",
                            mod_data$Inoc_Tag)
mod_data$Inoc_Tag <- ifelse(mod_data$Location %in% EUGENE1, "21420",
                            mod_data$Inoc_Tag)
mod_data$Inoc_Tag <- ifelse(mod_data$Location %in% EUGENE2, "24466",
                            mod_data$Inoc_Tag)

mod_data$Inoc_Tag <- ifelse(mod_data$Location %in% GUAPST2, "22575",
                            mod_data$Inoc_Tag)
mod_data$Inoc_Tag <- ifelse(mod_data$Location %in% GUAPST3, "24460",
                            mod_data$Inoc_Tag)
mod_data$Inoc_Tag <- ifelse(mod_data$Location %in% GUAPST1, "20507_2",
                            mod_data$Inoc_Tag)

mod_data$Inoc_Tag <- ifelse(mod_data$Location %in% HEISCO1, "20732",
                            mod_data$Inoc_Tag)
mod_data$Inoc_Tag <- ifelse(mod_data$Location %in% HEISCO2, "46666",
                            mod_data$Inoc_Tag)
mod_data$Inoc_Tag <- ifelse(mod_data$Location %in% HEISCO3, "47183_2",
                            mod_data$Inoc_Tag)

#Now it's time to bring them together
mod_data <- left_join(mod_data, geochem, by = c("Inoc_Sp", "Inoc_Tag"))

#Remove the Inoc_Tag with misidentified moisture treatments
mod_data <- mod_data %>% filter(Inoc_Tag != "22401")

#save the file
write.csv(mod_data, file = "Data/mod_data_allpots.csv")

#### Time Series Data Prep ------------------------------------------------------
#Fix a misnamed value
dat.ts[1696,1] <- "A"

dat.ts <- transform(dat.ts, Height = as.numeric(as.character(Height)), 
                    LeafNo = as.numeric(as.character(LeafNo)) 
                    
)
#Add Duration Column
dat.ts$Begin <- dmy(dat.ts$Begin)
dat.ts$CensusDate<- dmy(dat.ts$CensusDate)

exp.interval<- interval(dat.ts$Begin, dat.ts$CensusDate) #calculate interval from start of experiment

dat.ts$Duration <- time_length(exp.interval, unit = "day") #days within interval

#Add Treatment

# Create a location column with both table and position accounted for (random effect)
dat.ts <- mutate(dat.ts, Location = str_c(Table, Position, sep = "_"))

Inoc_Sp <- data.frame(mod_data$Location, mod_data$Inoc_Sp)#extract Inoc_Sp data from other dataset
colnames(Inoc_Sp) <- c("Location", "Inoc_Sp") #make sure column names are the same
dat.ts <- left_join(dat.ts, Inoc_Sp, by = "Location")
dat.ts <- dat.ts %>% mutate(soil = case_when(Inoc_Tag == "Control" ~ "control",
                                             (Seedling == Inoc_Sp) & (Inoc_Tag != "Control") ~ "conspecific",
                                             Seedling != Inoc_Sp ~ "heterospecific"))


# Census Interval
ts.ints <- dat.ts %>% group_by(Location) %>%
  arrange(CensusDate, .by_group = TRUE) %>%
  mutate(censusint = CensusDate - lag(CensusDate, default = first(CensusDate)))

ts.ints$censusint <- as.numeric(ts.ints$censusint)
#Create a column of height and LeafNo from the previous census:
prevht <- ts.ints %>% group_by(Location) %>%
  arrange(CensusDate, .by_group = TRUE) %>%
  mutate(htprevcensus = lag(Height))

prevlf <- prevht %>% group_by(Location) %>%
  arrange(CensusDate, .by_group = TRUE) %>%
  mutate(lfprevcensus = lag(LeafNo))

#Survival between census period
prevlf <- mutate(prevlf, status = ifelse(Height != 999, 0, 1))

prevsurv <- prevlf %>% group_by(Location) %>%
  arrange(CensusDate, .by_group = TRUE) %>%
  mutate(statusprevcensus = lag(status))

# remove mis identified soil source

prevsurv <- prevsurv %>% filter(Inoc_Tag != "22401")

write.csv(prevsurv, file = "Data/TimeSeriesData.csv")
