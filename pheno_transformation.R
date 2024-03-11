
# Read in file and set up work directory
setwd()

library(dplyr)
library(ggplot2)

getwd()

# Read in phenotypic file
pheno = read.table("./PAHmerge_phenotype.txt", header=TRUE)

# Data exploration-------------------------------------------------------------- 
str(pheno)
summary(pheno)


# Histogram of MPAPdiff
png('hist_density_MPAPdiff.png', width=6, height=6, units ='in', res=1000)
hist(pheno$MPAPdiff, 
     breaks=50,
     col="peachpuff",
     probability=TRUE,
     ylim=c(0, 0.1),
     xlab="MPAPdiff",
     main="Histogram and Density plot") # Shows density instead of frequency
lines(density(pheno$MPAPdiff),
      lwd=2,
      col = "black") # col = "chocolate3"
dev.off()

# QQplot
png('QQ plot of MPAPdiff.png', width=6, height=6, units ='in', res=1000)
qqnorm(pheno$MPAPdiff, cex=0.7, main="QQ plot of MPAPdiff") # Checking QQplot
qqline(pheno$MPAPdiff, col="red")
dev.off()




# Data transformation-----------------------------------------------------------

# Remove outliers 
cleaned_pheno = pheno[pheno$MPAPdiff<60 & pheno$MPAPdiff> -40,]

# Standardize data
cleaned_standard_pheno <- cleaned_pheno %>% 
  mutate(centered_MPAPdiff = MPAPdiff - mean(MPAPdiff)) %>% 
  mutate(standardized_MPAPdiff = centered_MPAPdiff / sd(MPAPdiff)) %>%
  mutate(normalized_MPAPdiff = (MPAPdiff - min(MPAPdiff)) / (max(MPAPdiff) - min(MPAPdiff)))

# Histogram
png('hist_density_standard_outliersRemoved_MPAPdiff.png', width=6, height=6, units ='in', res=1000)
hist(cleaned_standard_pheno$standardized_MPAPdiff, 
     breaks=50,
     col="peachpuff",
     probability=TRUE,
     ylim=c(0, 1.0),
     xlab="Standardized MPAPdiff after removing outliers",
     main="Histogram and Density plot of standardized data")
# Shows density instead of frequency
lines(density(cleaned_standard_pheno$standardized_MPAPdiff),
      lwd=2,
      col = "black")
dev.off()

# QQ plot
png('QQ plot of standardized MPAPdiff.png', width=6, height=6, units ='in', res=1000)
qqnorm(cleaned_standard_pheno$standardized_MPAPdiff, cex=1) # QQplot of the cleaned and standardized data
qqline(cleaned_standard_pheno$standardized_MPAPdiff, col="red")
dev.off()