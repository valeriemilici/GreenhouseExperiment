rm(list = ls())

library(tidyverse) #always tidyverse
library(lubridate) #convert dates to readable form

dat <- read.csv("Data/Chapter2_Data_Clean.csv")
dat.ts <- read.csv("Data/Ch2.TimeSeries.Data .csv")

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

#Create a column that identifies the control pots
mod_data <- mod_data %>% mutate(soil = case_when(Inoc_Sp_No == "C" ~ "control",
                                                 (Seedling == Inoc_Sp) & (Inoc_Sp_No != "C") ~ "conspecific",
                                                 Seedling != Inoc_Sp ~ "heterospecific"))

#relevel so that control is reference
#mod_data$soil <- relevel(as.factor(mod_data$soil), ref = "control")
## VM notes: perhaps better to retain conspecific as the control to see how control and het differ from 
## cons soil treatment. 

#Create a column that identifies a pot as CW, CD, HW, HD
mod_data<- mutate(mod_data, trt = str_c(soil, Treatment, sep = "_"))

#Create a "Status" column that denotes if a seedling is alive or dead at the end of the experiment
mod_data$status <- ifelse(!is.na(mod_data$Final_Height), "1", "0")

#save the file
write.csv(mod_data, file = "Data/mod_data_allpots.csv")

###############################
#### Time Series Data Prep ####
###############################
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


write.csv(prevsurv, file = "Data/TimeSeriesData.csv")
