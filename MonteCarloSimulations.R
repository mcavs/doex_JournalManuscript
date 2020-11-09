# -------------------------------------------------------------------------------------------------
# Script to reproduce the results on paper:
# Testing the equality of normal distributed groups' means under unequal variances by doex package

# Created by: Mustafa Cavus in 03/11/2020
# Contact: mustafacavus@eskisehir.edu.tr
# -------------------------------------------------------------------------------------------------
# Verify the packages
library(doex)
library(ggplot2)
#################################################################
# Main simulation functions for k = 3, 5, 7 groups
#################################################################
comp_anova_k3 = function(n1 = 10,
                         n2 = 10,
                         n3 = 10,
                         m1 = 0,
                         m2 = 0,
                         m3 = 0,
                         v1 = 0.1,
                         v2 = 0.2,
                         v3 = 0.3,
                         rept = 10){
  
  # Vectors for storing the p-values in for loop
  p.ag  <- numeric(rept)
  p.agf <- numeric(rept)
  p.af  <- numeric(rept)
  p.bx  <- numeric(rept)
  p.bf  <- numeric(rept)
  p.b2  <- numeric(rept)
  p.cf  <- numeric(rept)
  p.fa  <- numeric(rept)
  p.gf  <- numeric(rept)
  p.jf  <- numeric(rept)
  p.mbf <- numeric(rept)
  p.mw  <- numeric(rept)
  p.pb  <- numeric(rept)
  p.pf  <- numeric(rept)
  p.ss  <- numeric(rept)
  p.we  <- numeric(rept)
  p.wa  <- numeric(rept)
  
  for(i in 1:rept){
    
    # Setting a seed for producible results
    set.seed(i) 
    
    # Generating samples for k groups
    sample <- c(rnorm(n1, m1, (v1) ^ 2),
                rnorm(n2, m2, (v2) ^ 2),
                rnorm(n3, m3, (v3) ^ 2))
    
    # Creating the group labels for k groups
    group  <- c(rep(1, n1),
                rep(2, n2),
                rep(3, n3))
    
    # p-values of the tests
    p.ag[i]  <- AG(sample, group)[3]        # Alexander-Govern
    p.agf[i] <- AGF(sample, group,rept)[1]  # Alvandi et al. GF
    p.af[i]  <- AF(sample, group)[4]        # Approximate F
    p.bx[i]  <- BX(sample, group)[4]        # Box F
    p.bf[i]  <- BF(sample, group)[3]        # Brown-Forsythe 
    p.b2[i]  <- B2(0.05, sample, group)[3]  # B2
    p.cf[i]  <- CF(sample, group)[3]        # Cochran F
    p.fa[i]  <- FA(sample, group, rept)[1]  # Fiducial Approach 
    p.gf[i]  <- GF(sample, group, rept)[1]  # Generalized F
    p.jf[i]  <- JF(sample, group)[4]        # Johansen F
    p.mbf[i] <- MBF(sample, group)[4]       # Modified Brown-Forsythe
    p.mw[i]  <- MW(sample, group)[3]        # Modified Welch
    p.pb[i]  <- PB(sample, group, rept)[1]  # Parametric bootstrap
    p.pf[i]  <- PF(sample, group, rept)[1]  # Permutation F
    p.ss[i]  <- SS(sample, group)[3]        # Scott-Smith
    p.we[i]  <- WE(sample, group)[3]        # Welch
    p.wa[i]  <- WA(sample, group)[3]        # Welch-Aspin
  }
  
  # Emprical Type I error probability / power of the tests 
  
  return(c("ag" = mean(p.ag  < 0.05),
           "agf"= mean(p.agf < 0.05),
           "af" = mean(p.af  < 0.05),
           "bx" = mean(p.bx  < 0.05),
           "bf" = mean(p.bf  < 0.05),
           "b2" = mean(p.b2  < 0.05),
           "cf" = mean(p.cf  < 0.05),
           "fa" = mean(p.fa  < 0.05),
           "gf" = mean(p.gf  < 0.05),
           "jf" = mean(p.jf  < 0.05),
           "mbf"= mean(p.mbf < 0.05),
           "mw" = mean(p.mw  < 0.05),
           "pb" = mean(p.pb  < 0.05),
           "pf" = mean(p.pf  < 0.05),
           "ss" = mean(p.ss  < 0.05),
           "we" = mean(p.we  < 0.05),
           "wa" = mean(p.wa  < 0.05)))
}

comp_anova_k5 = function(n1 = 10,
                         n2 = 10,
                         n3 = 10,
                         n4 = 10,
                         n5 = 10,
                         m1 = 0,
                         m2 = 0,
                         m3 = 0,
                         m4 = 0,
                         m5 = 0,
                         v1 = 0.1,
                         v2 = 0.2,
                         v3 = 0.3,
                         v4 = 0.4,
                         v5 = 0.5,
                         rept = 10){
  
  # Vectors for storing the p-values in for loop
  p.ag  <- numeric(rept)
  p.agf <- numeric(rept)
  p.af  <- numeric(rept)
  p.bx  <- numeric(rept)
  p.bf  <- numeric(rept)
  p.b2  <- numeric(rept)
  p.cf  <- numeric(rept)
  p.fa  <- numeric(rept)
  p.gf  <- numeric(rept)
  p.jf  <- numeric(rept)
  p.mbf <- numeric(rept)
  p.mw  <- numeric(rept)
  p.pb  <- numeric(rept)
  p.pf  <- numeric(rept)
  p.ss  <- numeric(rept)
  p.we  <- numeric(rept)
  p.wa  <- numeric(rept)
  
  for(i in 1:rept){
    
    # Setting a seed for producible results
    set.seed(i)
    
    # Generating samples for k groups
    sample <- c(rnorm(n1, m1, (v1) ^ 2),
                rnorm(n2, m2, (v2) ^ 2),
                rnorm(n3, m3, (v3) ^ 2),
                rnorm(n4, m4, (v4) ^ 2),
                rnorm(n5, m5, (v5) ^ 2))
    
    # Creating the group labels for k groups
    group  <- c(rep(1, n1),
                rep(2, n2),
                rep(3, n3),
                rep(4, n4),
                rep(5, n5))
    
    # p-values of the tests
    p.ag[i]  <- AG(sample, group)[3]        # Alexander-Govern
    p.agf[i] <- AGF(sample, group,rept)[1]  # Alvandi et al. GF
    p.af[i]  <- AF(sample, group)[4]        # Approximate F
    p.bx[i]  <- BX(sample, group)[4]        # Box F
    p.bf[i]  <- BF(sample, group)[3]        # Brown-Forsythe 
    p.b2[i]  <- B2(0.05, sample, group)[3]  # B2
    p.cf[i]  <- CF(sample, group)[3]        # Cochran F
    p.fa[i]  <- FA(sample, group, rept)[1]  # Fiducial Approach 
    p.gf[i]  <- GF(sample, group, rept)[1]  # Generalized F
    p.jf[i]  <- JF(sample, group)[4]        # Johansen F
    p.mbf[i] <- MBF(sample, group)[4]       # Modified Brown-Forsythe
    p.mw[i]  <- MW(sample, group)[3]        # Modified Welch
    p.pb[i]  <- PB(sample, group, rept)[1]  # Parametric bootstrap
    p.pf[i]  <- PF(sample, group, rept)[1]  # Permutation F
    p.ss[i]  <- SS(sample, group)[3]        # Scott-Smith
    p.we[i]  <- WE(sample, group)[3]        # Welch
    p.wa[i]  <- WA(sample, group)[3]        # Welch-Aspin
  }
  
  # Emprical Type I error probability / power of the tests 
  return(c("ag" = mean(p.ag  < 0.05),
           "agf"= mean(p.agf < 0.05),
           "af" = mean(p.af  < 0.05),
           "bx" = mean(p.bx  < 0.05),
           "bf" = mean(p.bf  < 0.05),
           "b2" = mean(p.b2  < 0.05),
           "cf" = mean(p.cf  < 0.05),
           "fa" = mean(p.fa  < 0.05),
           "gf" = mean(p.gf  < 0.05),
           "jf" = mean(p.jf  < 0.05),
           "mbf"= mean(p.mbf < 0.05),
           "mw" = mean(p.mw  < 0.05),
           "pb" = mean(p.pb  < 0.05),
           "pf" = mean(p.pf  < 0.05),
           "ss" = mean(p.ss  < 0.05),
           "we" = mean(p.we  < 0.05),
           "wa" = mean(p.wa  < 0.05)))
}

comp_anova_k7 = function(n1 = 10,
                         n2 = 10,
                         n3 = 10,
                         n4 = 10,
                         n5 = 10,
                         n6 = 10,
                         n7 = 10,
                         m1 = 0,
                         m2 = 0,
                         m3 = 0,
                         m4 = 0,
                         m5 = 0,
                         m6 = 0,
                         m7 = 0,
                         v1 = 0.1,
                         v2 = 0.2,
                         v3 = 0.3,
                         v4 = 0.4,
                         v5 = 0.5,
                         v6 = 0.6,
                         v7 = 0.7,
                         rept = 10){
  
  # Vectors for storing the p-values in for loop
  p.ag  <- numeric(rept)
  p.agf <- numeric(rept)
  p.af  <- numeric(rept)
  p.bx  <- numeric(rept)
  p.bf  <- numeric(rept)
  p.b2  <- numeric(rept)
  p.cf  <- numeric(rept)
  p.fa  <- numeric(rept)
  p.gf  <- numeric(rept)
  p.jf  <- numeric(rept)
  p.mbf <- numeric(rept)
  p.mw  <- numeric(rept)
  p.pb  <- numeric(rept)
  p.pf  <- numeric(rept)
  p.ss  <- numeric(rept)
  p.we  <- numeric(rept)
  p.wa  <- numeric(rept)
  
  for(i in 1:rept){
    
    # Setting a seed for producible results
    set.seed(i)
    
    # Generating samples for k groups
    sample <- c(rnorm(n1, m1, (v1) ^ 2),
                rnorm(n2, m2, (v2) ^ 2),
                rnorm(n3, m3, (v3) ^ 2),
                rnorm(n4, m4, (v4) ^ 2),
                rnorm(n5, m5, (v5) ^ 2),
                rnorm(n6, m6, (v6) ^ 2),
                rnorm(n7, m7, (v7) ^ 2))
    
    # Creating the group labels for k groups
    group  <- c(rep(1, n1),
                rep(2, n2),
                rep(3, n3),
                rep(4, n4),
                rep(5, n5),
                rep(6, n6),
                rep(7, n7))
    
    # p-values of the tests
    p.ag[i]  <- AG(sample, group)[3]        # Alexander-Govern
    p.agf[i] <- AGF(sample, group,rept)[1]  # Alvandi et al. GF
    p.af[i]  <- AF(sample, group)[4]        # Approximate F
    p.bx[i]  <- BX(sample, group)[4]        # Box F
    p.bf[i]  <- BF(sample, group)[3]        # Brown-Forsythe 
    p.b2[i]  <- B2(0.05, sample, group)[3]  # B2
    p.cf[i]  <- CF(sample, group)[3]        # Cochran F
    p.fa[i]  <- FA(sample, group, rept)[1]  # Fiducial Approach 
    p.gf[i]  <- GF(sample, group, rept)[1]  # Generalized F
    p.jf[i]  <- JF(sample, group)[4]        # Johansen F
    p.mbf[i] <- MBF(sample, group)[4]       # Modified Brown-Forsythe
    p.mw[i]  <- MW(sample, group)[3]        # Modified Welch
    p.pb[i]  <- PB(sample, group, rept)[1]  # Parametric bootstrap
    p.pf[i]  <- PF(sample, group, rept)[1]  # Permutation F
    p.ss[i]  <- SS(sample, group)[3]        # Scott-Smith
    p.we[i]  <- WE(sample, group)[3]        # Welch
    p.wa[i]  <- WA(sample, group)[3]        # Welch-Aspin
  }
  
  # Emprical Type I error probability / power of the tests 
  return(c("ag" = mean(p.ag  < 0.05),
           "agf"= mean(p.agf < 0.05),
           "af" = mean(p.af  < 0.05),
           "bx" = mean(p.bx  < 0.05),
           "bf" = mean(p.bf  < 0.05),
           "b2" = mean(p.b2  < 0.05),
           "cf" = mean(p.cf  < 0.05),
           "fa" = mean(p.fa  < 0.05),
           "gf" = mean(p.gf  < 0.05),
           "jf" = mean(p.jf  < 0.05),
           "mbf"= mean(p.mbf < 0.05),
           "mw" = mean(p.mw  < 0.05),
           "pb" = mean(p.pb  < 0.05),
           "pf" = mean(p.pf  < 0.05),
           "ss" = mean(p.ss  < 0.05),
           "we" = mean(p.we  < 0.05),
           "wa" = mean(p.wa  < 0.05)))
}
#################################################################
#################################################################
# Figure 1: Type I error probability of the tests for k = 3
#################################################################
error_k3 <- matrix(0, 24, 17)
error_k3[1,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3)
error_k3[2,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7)
error_k3[3,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 2,   v3 = 3)
error_k3[4,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 4,   v3 = 7)

error_k3[5,]  <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3)
error_k3[6,]  <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7)
error_k3[7,]  <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 2,   v3 = 3)
error_k3[8,]  <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 4,   v3 = 7)

error_k3[9,]  <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3)
error_k3[10,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7)
error_k3[11,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 2,   v3 = 3)
error_k3[12,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 4,   v3 = 7)

error_k3[13,] <- comp_anova_k3(n1 = 5 , n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3)
error_k3[14,] <- comp_anova_k3(n1 = 5 , n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7)
error_k3[15,] <- comp_anova_k3(n1 = 5 , n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 2,   v3 = 3)
error_k3[16,] <- comp_anova_k3(n1 = 5 , n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 4,   v3 = 7)

error_k3[17,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3)
error_k3[18,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7)
error_k3[19,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 2,   v3 = 3)
error_k3[20,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 4,   v3 = 7)

error_k3[21,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3)
error_k3[22,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7)
error_k3[23,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 2,   v3 = 3)
error_k3[24,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0, v1 = 1,   v2 = 4,   v3 = 7)
# Creating Figure 1
colnames(error_k3) <- c("AF", "AG", "AGF", "B2", "BF", "BX", "CF", "FA", "GF", "JF", "MBF", "MW", "PB", "PF", "SS", "WA", "WE")
boxplot(error_k3, ylab = "Type I error probability", xlab = "Tests", main = "Type I error probability of the tests for k = 3")
abline(h = 0.045, col = "Red", lty = 2) # Bradley (1970)'s robustness lower limit
abline(h = 0.055, col = "Red", lty = 2) # Bradley (1970)'s robustness upper limit
#################################################################
#################################################################
# Figure 2: Type I error probability of the tests for k = 5
#################################################################
error_k5 <- matrix(0, 24, 17)
error_k5[1,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
error_k5[2,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
error_k5[3,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
error_k5[4,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)

error_k5[5,]  <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
error_k5[6,]  <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
error_k5[7,]  <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
error_k5[8,]  <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)

error_k5[9,]  <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
error_k5[10,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
error_k5[11,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
error_k5[12,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)

error_k5[13,] <- comp_anova_k5(n1 = 4 , n2 = 6 , n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
error_k5[14,] <- comp_anova_k5(n1 = 4 , n2 = 6 , n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
error_k5[15,] <- comp_anova_k5(n1 = 4 , n2 = 6 , n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
error_k5[16,] <- comp_anova_k5(n1 = 4 , n2 = 6 , n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)

error_k5[17,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
error_k5[18,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
error_k5[19,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
error_k5[20,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)

error_k5[21,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
error_k5[22,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
error_k5[23,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
error_k5[24,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
# Creating Figure 2
colnames(error_k5) <- c("AF", "AG", "AGF", "B2", "BF", "BX", "CF", "FA", "GF", "JF", "MBF", "MW", "PB", "PF", "SS", "WA", "WE")
boxplot(error_k5, ylab = "Type I error probability", xlab = "Tests", main = "Type I error probability of the tests for k = 5")
abline(h = 0.045, col = "Red", lty = 2) # Bradley (1970)'s robustness lower limit
abline(h = 0.055, col = "Red", lty = 2) # Bradley (1970)'s robustness upper limit
#################################################################
#################################################################
# Figure 3: Type I error probability of the tests for k = 7
#################################################################
error_k7 <- matrix(0, 24, 17)
error_k7[1,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
error_k7[2,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
error_k7[3,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
error_k7[4,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)

error_k7[5,]  <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
error_k7[6,]  <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
error_k7[7,]  <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
error_k7[8,]  <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)

error_k7[9,]  <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
error_k7[10,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
error_k7[11,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
error_k7[12,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)

error_k7[13,] <- comp_anova_k7(n1 = 4,  n2 = 6,  n3 = 8,  n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
error_k7[14,] <- comp_anova_k7(n1 = 4,  n2 = 6,  n3 = 8,  n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
error_k7[15,] <- comp_anova_k7(n1 = 4,  n2 = 6,  n3 = 8,  n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
error_k7[16,] <- comp_anova_k7(n1 = 4,  n2 = 6,  n3 = 8,  n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)

error_k7[17,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
error_k7[18,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
error_k7[19,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
error_k7[20,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)

error_k7[21,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
error_k7[22,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
error_k7[23,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
error_k7[24,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
# Creating Figure 3
colnames(error_k7) <- c("AF", "AG", "AGF", "B2", "BF", "BX", "CF", "FA", "GF", "JF", "MBF", "MW", "PB", "PF", "SS", "WA", "WE")
boxplot(error_k7, ylab = "Type I error probability", xlab = "Tests", main = "Type I error probability of the tests for k = 7")
abline(h = 0.045, col = "Red", lty = 2) # Bradley (1970)'s robustness lower limit
abline(h = 0.055, col = "Red", lty = 2) # Bradley (1970)'s robustness upper limit
#################################################################
#################################################################
# Table 1: Penalized powers for k = 3
#################################################################
## n_i = (10, 10, 10)
#################################################################
power_k3 <- matrix(0, 72, 17)
power_k3[1,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[2,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[3,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3)

power_k3[4,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[5,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[6,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7)

power_k3[7,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[8,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[9,]  <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 2,   v3 = 3)

power_k3[10,] <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[11,] <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[12,] <- comp_anova_k3(n1 = 10, n2 = 10, n3 = 10, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 4,   v3 = 7)
#################################################################
## n_i = (30, 30, 30)
#################################################################
power_k3[13,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[14,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[15,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3)

power_k3[16,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[17,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[18,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7)

power_k3[19,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[20,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[21,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 2,   v3 = 3)

power_k3[22,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[23,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[24,] <- comp_anova_k3(n1 = 30, n2 = 30, n3 = 30, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 4,   v3 = 7)
#################################################################
## n_i = (50, 50, 50)
#################################################################
power_k3[25,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[26,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[27,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3)

power_k3[28,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[29,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[30,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7)

power_k3[31,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[32,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[33,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 2,   v3 = 3)

power_k3[34,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[35,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[36,] <- comp_anova_k3(n1 = 50, n2 = 50, n3 = 50, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 4,   v3 = 7)
#################################################################
## n_i = (5, 10, 15)
#################################################################
power_k3[37,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[38,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[39,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3)

power_k3[40,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[41,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[42,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7)

power_k3[43,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[44,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[45,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 2,   v3 = 3)

power_k3[46,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[47,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[48,] <- comp_anova_k3(n1 = 5, n2 = 10, n3 = 15, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 4,   v3 = 7)
#################################################################
## n_i = (20, 30, 40)
#################################################################
power_k3[49,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[50,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[51,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3)

power_k3[52,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[53,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[54,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7)

power_k3[55,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[56,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[57,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 2,   v3 = 3)

power_k3[58,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[59,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[60,] <- comp_anova_k3(n1 = 20, n2 = 30, n3 = 40, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 4,   v3 = 7)
#################################################################
## n_i = (25, 50, 75)
#################################################################
power_k3[61,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[62,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3)
power_k3[63,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3)

power_k3[64,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[65,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7)
power_k3[66,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7)

power_k3[67,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[68,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 2,   v3 = 3)
power_k3[69,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 2,   v3 = 3)

power_k3[70,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0.3, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[71,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 0.8, v1 = 1,   v2 = 4,   v3 = 7)
power_k3[72,] <- comp_anova_k3(n1 = 25, n2 = 50, n3 = 75, m1 = 0, m2 = 0, m3 = 1.5, v1 = 1,   v2 = 4,   v3 = 7)
# Creating Table 1
n_col_k3 <- c(rep(c("(10, 10, 10)", "(30, 30, 30)", "(50, 50, 50)",
                    "(5, 10, 15)", "(20, 30, 40)", "(25, 50, 75)"), each = 12))

var_col_k3 <- c(rep(rep(c("(0.1, 0.2, 0.3)",
                      "(0.1, 0.4, 0.7)",
                      "(1, 2, 3)",
                      "(1, 4, 7)"), each = 3), 6))

es_col_k3 <- c(rep(c("0.3", "0.8", "1.5"), 24))

power_k3 <- as.data.frame(power_k3)
power_k3_table <- data.frame(n_col_k3, var_col_k3, es_col_k3, power_k3)
colnames(power_k3_table) <- c("n_i", "sigma_i^2", "/Delta_i", "AF", "AG", "AGF", "B2", "BF", "BX", "CF", "FA", "GF", "JF", "MBF", "MW", "PB", "PF", "SS", "WA", "WE")
power_k3_table
#################################################################
#################################################################
# Table 2: Penalized powers for k = 5
#################################################################
## n_i = (10, 10, 10, 10, 10)
#################################################################
power_k5 <- matrix(0, 72, 17)
power_k5[1,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[2,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[3,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)

power_k5[4,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[5,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[6,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)

power_k5[7,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[8,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[9,]  <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)

power_k5[10,] <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[11,] <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[12,] <- comp_anova_k5(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
#################################################################
## n_i = (30, 30, 30, 30, 30)
#################################################################
power_k5[13,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[14,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[15,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)

power_k5[16,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[17,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[18,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)

power_k5[19,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[20,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[21,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)

power_k5[22,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[23,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[24,] <- comp_anova_k5(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
#################################################################
## n_i = (50, 50, 50, 50, 50)
#################################################################
power_k5[25,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[26,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[27,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)

power_k5[28,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[29,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[30,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)

power_k5[31,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[32,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[33,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)

power_k5[34,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[35,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[36,] <- comp_anova_k5(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
#################################################################
## n_i = (4, 6, 10, 14, 16)
#################################################################
power_k5[37,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[38,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[39,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)

power_k5[40,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[41,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[42,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)

power_k5[43,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[44,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[45,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)

power_k5[46,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[47,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[48,] <- comp_anova_k5(n1 = 4, n2 = 6, n3 = 10, n4 = 14, n5 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
#################################################################
## n_i = (12, 18, 30, 42, 48)
#################################################################
power_k5[49,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[50,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[51,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)

power_k5[52,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[53,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[54,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)

power_k5[55,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[56,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[57,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)

power_k5[58,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[59,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[60,] <- comp_anova_k5(n1 = 12, n2 = 18, n3 = 30, n4 = 42, n5 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
#################################################################
## n_i = (20, 30, 50, 70, 80)
#################################################################
power_k5[61,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[62,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)
power_k5[63,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5)

power_k5[64,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[65,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)
power_k5[66,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5)

power_k5[67,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[68,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)
power_k5[69,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5)

power_k5[70,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[71,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
power_k5[72,] <- comp_anova_k5(n1 = 20, n2 = 30, n3 = 50, n4 = 70, n5 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15)
# Creating Table 2
n_col_k5 <- c(rep(c("(10, 10, 10, 10, 10)", 
                    "(30, 30, 30, 30, 30)", 
                    "(50, 50, 50, 50, 50)",
                    "(4, 6, 10, 14, 16)", 
                    "(12, 18, 30, 42, 48)", 
                    "(20, 30, 50, 70, 80)"), each = 12))

var_col_k5 <- c(rep(rep(c("(0.1, 0.2, 0.3, 0.4, 0.5)",
                          "(0.1, 0.4, 0.7, 1.1, 1.5)",
                          "(1, 2, 3, 4, 5)",
                          "(1, 4, 7, 11, 15)"), each = 3), 6))

es_col_k5 <- c(rep("0.3", "0.8", "1.5"), 24)

power_k5 <- as.data.frame(power_k5)
power_k5_table <- data.frame(n_col_k5, var_col_k5, es_col_k5, power_k5)
colnames(power_k5_table) <- c("n_i", "sigma_i^2", "/Delta_i", "AF", "AG", "AGF", "B2", "BF", "BX", "CF", "FA", "GF", "JF", "MBF", "MW", "PB", "PF", "SS", "WA", "WE")
power_k5_table
#################################################################
#################################################################
# Table 3: Penalized powers for k = 7
#################################################################
## n_i = (10, 10, 10, 10, 10, 10, 10)
#################################################################
power_k7 <- matrix(0, 72, 17)
power_k7[1,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[2,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[3,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)

power_k7[4,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[5,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[6,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)

power_k7[7,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[8,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[9,]  <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)

power_k7[10,] <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[11,] <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[12,] <- comp_anova_k7(n1 = 10, n2 = 10, n3 = 10, n4 = 10, n5 = 10, n6 = 10, n7 = 10, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
#################################################################
## n_i = (30, 30, 30, 30, 30, 30, 30)
#################################################################
power_k7[13,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[14,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[15,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)

power_k7[16,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[17,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[18,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)

power_k7[19,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[20,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[21,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)

power_k7[22,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[23,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[24,] <- comp_anova_k7(n1 = 30, n2 = 30, n3 = 30, n4 = 30, n5 = 30, n6 = 30, n7 = 30, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
#################################################################
## n_i = (50, 50, 50, 50, 50, 50, 50)
#################################################################
power_k7[25,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[26,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[27,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)

power_k7[28,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[29,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[30,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)

power_k7[31,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[32,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[33,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)

power_k7[34,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[35,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[36,] <- comp_anova_k7(n1 = 50, n2 = 50, n3 = 50, n4 = 50, n5 = 50, n6 = 50, n7 = 50, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
#################################################################
## n_i = (4, 6, 8, 10, 12, 14, 16)
#################################################################
power_k7[37,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[38,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[39,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)

power_k7[40,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[41,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[42,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)

power_k7[43,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[44,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[45,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)

power_k7[46,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[47,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[48,] <- comp_anova_k7(n1 = 4, n2 = 6, n3 = 8, n4 = 10, n5 = 12, n6 = 14, n7 = 16, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
#################################################################
## n_i = (12, 18, 24, 30, 36, 42, 48)
#################################################################
power_k7[49,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[50,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[51,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)

power_k7[52,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[53,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[54,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)

power_k7[55,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[56,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[57,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)

power_k7[58,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[59,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[60,] <- comp_anova_k7(n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36, n6 = 42, n7 = 48, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
#################################################################
## n_i = (20, 30, 40, 50, 60, 70, 80)
#################################################################
power_k7[61,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[62,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)
power_k7[63,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.2, v3 = 0.3, v4 = 0.4, v5 = 0.5, v6 = 0.6, v7 = 0.7)

power_k7[64,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[65,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)
power_k7[66,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 0.1, v2 = 0.4, v3 = 0.7, v4 = 1.1, v5 = 1.5, v6 = 1.9, v7 = 2.3)

power_k7[67,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[68,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)
power_k7[69,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 2,   v3 = 3,   v4 = 4,   v5 = 5,   v6 = 6,   v7 = 7)

power_k7[70,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.3, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[71,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 0.8, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
power_k7[72,] <- comp_anova_k7(n1 = 20, n2 = 30, n3 = 40, n4 = 50, n5 = 60, n6 = 70, n7 = 80, m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, m7 = 1.5, v1 = 1,   v2 = 4,   v3 = 7,   v4 = 11,  v5 = 15,  v6 = 19,  v7 = 23)
# Creating Table 3
n_col_k7 <- c(rep(c("(10, 10, 10, 10, 10, 10, 10)", 
                    "(30, 30, 30, 30, 30, 30, 30)", 
                    "(50, 50, 50, 50, 50, 50, 50)",
                    "(4, 6, 8, 10, 12, 14, 16)", 
                    "(12, 18, 24, 30, 36, 42, 48)", 
                    "(20, 30, 40, 50, 60, 70, 80)"), each = 12))

var_col_k7 <- c(rep(rep(c("(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.6)",
                          "(0.1, 0.4, 0.7, 1.1, 1.5, 1.9, 2.3)",
                          "(1, 2, 3, 4, 5, 6, 7)",
                          "(1, 4, 7, 11, 15, 19, 23)"), each = 3), 6))

es_col_k7 <- c(rep("0.3", "0.8", "1.5"), 24)

power_k7 <- as.data.frame(power_k7)

power_k7 <- as.data.frame(power_k7)
power_k7_table <- data.frame(n_col_k7, var_col_k7, es_col_k7, power_k7)
colnames(power_k7_table) <- c("n_i", "sigma_i^2", "/Delta_i", "AF", "AG", "AGF", "B2", "BF", "BX", "CF", "FA", "GF", "JF", "MBF", "MW", "PB", "PF", "SS", "WA", "WE")
power_k7_table
#################################################################
#################################################################