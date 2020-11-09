# -------------------------------------------------------------------------------------------------
# Dataset to reproduce the results on paper:
# Testing the equality of normal distributed groups' means under unequal variances by doex package

# Created by: Mustafa Cavus in 09/11/2020
# Contact: mustafacavus@eskisehir.edu.tr
# -------------------------------------------------------------------------------------------------

# Dataset
# Example: Litter Weight Dose Response
# A data set analyzed by Westfall and Rom (1990) involces
# litter weights of mice born from mothers assigned to three 
# different dosage groups and a control. For the low dose
# group the dose metameter is 5, for the medium dose 
# group it is 50, and for the high dose group it is 500.

# Verify the packages
library(doex)
library(car)

# average litter weight data
weight_data <- c(22.69, 26.59, 28.85, 28.03, 29.05,
                 23.61, 22.21, 26.81, 26.01, 25.98,
                 24.75, 26.60, 24.10, 23.09, 26.56,
                 26.58, 23.65, 26.19, 25.11, 28.18,
                 27.84, 21.45, 19.85, 30.95, 22.40,
                 26.95, 20.23, 26.46, 28.64, 21.48,
                 25.04, 24.18, 21.74, 25.64, 26.86,
                 17.39, 20.73, 16.34, 22.75, 24.80,
                 28.25, 22.33, 26.43, 24.50, 24.04,
                 21.71, 25.43, 29.21, 22.84, 17.54,
                 24.69, 24.44, 22.18, 18.79, 23.58,
                 24.18, 23.30, 19.55, 26.90, 26.38,
                 20.53, 24.10, 16.13, 21.11, 23.03,
                 16.26, 26.19, 20.99, 26.33, 26.31, 
                 30.61, 26.48, 24.31, 27.98)

# dose 
dose <- as.factor(c(rep("0", 20), rep("5", 19),
                    rep("50", 18), rep("500", 17)))

born_weight_data <- data.frame(weight_data, dose) 

# testing the variance homogeneity assumption
car::leveneTest(weight_data ~ dose)

# F-value = 3.3819, p-value: 0.0229*

# testing the equality of means of the weight data for dose treatments
doex::AF(weight_data, dose)




