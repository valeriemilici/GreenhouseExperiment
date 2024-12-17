# This files uses the geochemistry data from the soil samples to create PCAs
# of the variation in soil geochemistry. The PCA that maximizes the variation
# explained in Dim1 and Dim2 will be used for downstream analyses

# Initialize Workspace
rm(list =ls())

library(tidyverse) #clean data manipulation
library(lme4) #mixed effects models
library(lmerTest) #p-values from linear mixed models
library(performance) #model diagnostics
library(corrr) #for PCAs
library(ggcorrplot) #PCA
library(FactoMineR) #PCA
library(factoextra) #PCA
library(ggfortify) #PCA 

# load data
geodat <- read.csv("Data/Soil_Geochem.csv")

# fix one column (for aesthetics)
geodat$Sample.ID <- toupper(geodat$Sample.ID)
geodat$Sample.ID <- ifelse(geodat$Sample.ID == "PSIDIUM", "PSIDFR", geodat$Sample.ID)

# convert random character columns to numeric

geodat$m_B <- as.numeric(geodat$m_B)
geodat$m_P <- as.numeric(geodat$m_P)
geodat$m_B[is.na(geodat$m_B)] <- 0
geodat$m_P[is.na(geodat$m_P)] <- 0
# Give duplicate tags unique identifiers

geodat[3,2] <- "20507_2"
geodat[12,2] <- "47183_2"

# Fix missing variables

geodat[1,8] <- 2.08
geodat[1,9] <- 98.8

### PCA #1 - likely most important variables -------
# pH, N, P, K, C

geodat3 <- geodat %>% dplyr::select(1,2, LOI., pHCaCl2, m_P, m_K, NO2.NO3, NH4) 
colSums(is.na(geodat3))

geodat4 <- geodat3 %>% dplyr::select(!c(1:2))

# Make the PCA

geochem_pca2 <- PCA(geodat4)

# Plot the PCA
geo2 <- fviz_pca_biplot(geochem_pca2, 
                        col.ind = geodat3$Sample.ID, palette = "jco", 
                        addEllipses = F, label = "all",  repel = TRUE,
                        legend.title = "Species") +
  theme_classic(12)

geo2

### PCA w/o Nitrogen (true measured values for all samples) ---------------

geodat5 <- geodat %>% dplyr::select(LOI., pHCaCl2, m_P, m_K) 
colSums(is.na(geodat5))

# Make the PCA

geochem_pca3 <- PCA(geodat5)

# Plot the PCA
geo3 <- fviz_pca_biplot(geochem_pca3, 
                        col.ind = geodat3$Sample.ID, palette = "jco", 
                        addEllipses = F, label = "all",  repel = TRUE,
                        legend.title = "Species") +
  theme_classic(12)

geo3

#geo3 explains more variation than geo2

ggsave(filename = "Figures/soil_pca2.png", plot = geo3)
# extract the value of dim1 for the 20 individuals in the PCA
soil2 <- data.frame(geochem_pca3$ind$coord[,1])

info <- geodat3 %>% dplyr::select(1:2)

soil2 <- cbind(soil2, info)
colnames(soil2) <- c("soil.dim1", "Inoc_Sp", "Inoc_Tag")

centroids2 <- soil2 %>% group_by(Inoc_Sp) %>%
  summarize(soil.dim1 = mean(soil.dim1)) %>%
  mutate(Inoc_Tag = "control")

soil3 <- rbind(soil2, centroids2)

write.csv(soil3, "Data/soil_pca.csv")
