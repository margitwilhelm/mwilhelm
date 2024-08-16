### Analysis of west coast steenbras (Lithugnathus aureti) body growth
### using bootstrapped length-at-age (otolith) analyses
#
# This script contains analyses conducted for the manuscript:
#
# "Extremely slow somatic growth and intermittent recruitment of west coast 
# steenbras, an over-exploited, longevous Sparid, analysed with novel bootstrapped 
# methods"
#
# Authors of the manuscript:
#
# Margit Wilhelm 
# Arariky Shikongo 
# Angelika Veii
# Ralf Schwamborn
#
#
## Fishmethods GrowthTraject with Land crab ("Guaiamum") data ### -> Changed for steenbras data
# Version 10 (September 2023)
# New resampling function
# By Margit Wilhelm: mwilhelm@unam.na

install.packages('TMB', type = 'source') 
library(TMB) #Template Model Builder: A General Random Effect Tool Inspired by 'ADMB'
library(fishmethods)
library(TropFishR)
library(tidyverse) #For ggsave

#To clean memory
gc()  
gc(reset=T)
rm(list = ls())

#Data : West coast steenbras - southern population - TR growth increments 2004 - 2009
path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# import data# 4 columns for fishmethods GrowthTraject function
guaiamum1 <- read.csv("steenTR1b.csv", header=T) 

View(guaiamum1) 
# TR1: 75 growth increments 2004-2009 (at least 28 days, and leave all positives and negatives in); 
# Assume lengths were measured in total length USED for LFD also 
# Later -> convert to fork length with equation from current data
# TR1b = added 2020-2022 data and removed -13 and - 7.5 growth = 80 growth increments -> USE THIS

#growth_tagging function in TropFishR
guaiamum1$delta_t = guaiamum1$T2 - guaiamum1$T1
guaiamum1$delta_L = guaiamum1$L2 - guaiamum1$L1

#gumi2 =subset(guaiamum1, guaiamum1$delta_L>-1)

#View(gumi2)
#View(guaiamum1)

Munro2=growth_tagging(param = guaiamum1,"Munro", time_unit = "year", Linf_range=c(70,140))
Munro2
Munro2$reg_coeffs
Munro2$r2
Munro2$Linf
Munro2$K
Munro2$conf_int_K

#Function growth_tagging:
#growth_tagging(param, method, Linf_range = c(5, 600), time_unit = "year")
#param	
#a list consisting of following parameters:
#L1: length at tagging [cm],
#L2: length at recapture [cm],
#delta_t: time interval between tagging and recapture 
#(instead two vectors with t1 (age at tagging) and t2 (age at recapture) can be provided.
                                                      
#method	
#indicating which of following methods should be applied: "GullandHolt" or "Munro".
                                                      
#Linf_range	
#two values indicating the lower and upper limits of the range, in which the optimise searches for the Linf value with the best fit (lowest CV value ),
                                                      
#time_unit	
#indicating the unit of the time interval, either "year", "month", "week", or "day

#-------------------------------------------------------
# estimates the parameters of the VBGF and model AIC
# 
#Model 4 of Francis (1988)
#grotag
#List specifying the design of the model to estimate. 
#Use 1 to designate whether a parameter(s) should be estimated. 
#Type of parameters are: 
#nu=growth variability (1 parameter), 
#m=bias parameter of measurement error (1 parameter), 
#p=outlier probability (1 parameter), and
#sea=seasonal variation (2 parameters: u and w). 
#Model 1 of Francis is the default settings of 0 for nu, m, p and sea.
#sigma = v*mu, v = growth variability, 
#Mean growth is made up of galpha and gbeta

alp <- min(guaiamum1$L1)
bet <- max(guaiamum1$L2)
alp #24 + 3
bet #61.5 => fit failed -> Use 40, 50


vbfg.mod.1 <- with (guaiamum1,
        grotag(L1=L1, L2=L2, T1=T1, T2=T2,
               alpha= 40,
               beta = 50,
             design=list(nu=1,m=1,p=1,sea=0), #seasonal variation = 0
             stvalue=list(sigma=0.9,nu=0.4,m=-1,p=0.1,u=0,w=0),
             upper=list(sigma=5,nu=1,m=2,p=0.5,u=0,w=0),#keep seasonal, u, and w at 0
             lower=list(sigma=0,nu=0,m=-2,p=0,u=0,w=0),control=list(maxit=2e9)))

vbfg.mod.1

vbfg.mod.1$VBparms

## plot histogram of residuals
vbfg.mod.1$residuals

plot(vbfg.mod.1$residuals)
length(vbfg.mod.1$residuals) # 75 residuals #80
hist(vbfg.mod.1$residuals) #most residuals are negative
summary(vbfg.mod.1$residuals)
sd(vbfg.mod.1$residuals)

shapiro.test(vbfg.mod.1$residuals) #NOT normally distributed (alpha 40, beta 50) #residuals are not normal BUT normal with second fit #mostly negative residuals

hist(vbfg.mod.1$residuals, 
     xlab = "Residuals of the VBGF (cm/year)",
     ylab = "Frequency (N)", main ="",
     col= "grey")

hist(vbfg.mod.1$residuals, 
     xlab = "Residuals of the VBGF (cm/year)",
     ylab = "Frequency (N)", main ="",
     col= "grey",
     ylim= c(0,15),
     breaks =c(seq (-11, 5, by = 1))
     )
 text(4,27,"N= 80")

#-------------------------------------------------------
# vbgf curve - Steenbras data

guaiamum1m <- guaiamum1

temp<-guaiamum1 #[c(guaiamum1m$T2-guaiamum1m$T1)>0,]

#-------PLOT
#mar = c(3,3,1,1)#bottom, left, top, right Default = c(5.1, 4.1, 4.1, 2.1).
#oma = c(0,0,0,0) #outer margins
#Parameters         NA
#Linf       51.7422574
#K           0.1232942

opar <- par(mfrow = c(2,1), mar = c(1,1,1,1), oma = c(0,0,0,0), mgp = c(2,0.5,0), 
            tcl = -0.25,
            cex = 1)
growthTraject(0.1233, 62,
              lentag=temp$L1, lenrec=temp$L2,timelib=c(temp$T2-temp$T1), ylim = c(20,65), 
              xlim = c(2,45), 
              main = "Growth increments & one fitted VBGF curve", 
              ylab = "Fork length (cm)", xlab = "Relative age (years)" ) #K, Linf
growthTraject(0.069802155, 71.9587101,
              lentag=temp$L1, lenrec=temp$L2,timelib=c(temp$T2-temp$T1), ylim = c(20,65), 
              xlim = c(2,45), 
              main = "Growth increments & Holtzhausen and Kirchner, 2001a (S) fitted VBGF curve", ylab = "Fork length (cm)", xlab = "Relative age (years)" ) #K, Linf
dev.off()

#Holtzhausen & Kirchner South
#LS <- VBGF(list(Linf=71.9587101, K=0.069802155,  t0=-3.757304163), t=t) 
#Holtzhausen & Kirchner N
#LN <- VBGF(list(Linf=84.30894408, K=0.092343421,  t0=-2.44717737), t=t) 

ggsave(filename = 'growthIncrements_Curve.png') 

plot = growthTraject(0.15, 62,
lentag=temp$L1, lenrec=temp$L2,timelib=c(temp$T2-temp$T1), ylim = c(20,65), xlim = c(2,25), 
main = "Growth increments & one fitted VBGF curve", ylab = "Fork length (cm)", xlab = "Relative age (years)" ) #K, Linf
#device = 'png', width = 32, height = 18 ,units = 'cm', dpi = 300)

#75 growth increments (at least 28 days, and leave all positives and negatives in), otherwise don't fit
#ALPHA = 48, BETA = 53
#Linf       44.7103751
#K           0.4431283
#Seasonal = 0 makes no difference

#Linf       48.612876
#K           1.128268 -> including recent

#ALPHA = 40, BETA = 50 -> Redo bootstrap with this alpha and beta (but doubtful result)
#Linf       51.7422574
#K           0.1232942

####################
#### BOOTSTRAP #### 
# now bootstrap

# ?grotag
gmodelm <- with (guaiamum1m,
              grotag(L1=L1, L2=L2, T1=T1, T2=T2,
              alpha= 40,
              beta = 50,
              design=list(nu=1,m=1,p=1,sea=0), #seasonal variation = 0
              stvalue=list(sigma=0.9,nu=0.4,m=-1,p=0.1,u=0,w=0),
              upper=list(sigma=5,nu=1,m=2,p=0.5,u=0,w=0),#keep seasonal, u, and w at 0
              lower=list(sigma=0,nu=0,m=-2,p=0,u=0,w=0),control=list(maxit=2e9)))
gmodelm$VBparms

Linf <- gmodelm$VBparms$Estimate[2]  # Linf
Linf
K <- gmodelm$VBparms$Estimate[3]  # K
K

# now sample 80/2 = 40 increments 
ssize = nrow(guaiamum1m)
samplerows <- sample(ssize , 40)
newsubset <- guaiamum1m[samplerows,]
summary(newsubset)
nrow(newsubset)

# now apply grotag on the subset
alp <- 40 #min(newsubset$L1)+3
bet <- 50 #max(newsubset$L2)-3

gmodel2m <- with (newsubset,
                 grotag(L1=L1, L2=L2, T1=T1, T2=T2,alpha=alp,beta=bet,
                        design=list(nu=0,m=1,p=0,sea=0),
                        stvalue=list(sigma=0.9,nu=0.4,m=-1,p=0,u=0,w=0),
                        upper=list(sigma=5,nu=1,m=2,p=0.5,u=0,w=0),
                        lower=list(sigma=0,nu=0,m=-2,p=0,u=0,w=0),control=list(maxit=1e4)))

gmodel2m


Linf <- gmodel2m$VBparms$Estimate[2]  # Linf -> Doesn't work if fit failed
# Linf <- (- 1) * Linf
Linf

K <- gmodel2m$VBparms$Estimate[3]  # K
K

#vbgf curve - 
temp<-newsubset[c(newsubset$T2-newsubset$T1)>0,]
growthTraject(K,Linf,lentag=temp$L1, lenrec=temp$L2,timelib=c(temp$T2-temp$T1), ylim = c(0,80), 
              main = "subset, n = 40") #Doesn't work if fit failed. Which is a lot

## 3.a Define the full sampling function
resample.fun.full <- function(input.data, sssize) { 
  
  ssize = nrow(input.data)
  samplerows <- sample(ssize , sssize, replace = TRUE)
  newsubset <- input.data[samplerows,]
  summary(newsubset)
  nrow(newsubset)
  
  # now apply grotag on the subset
  alp <- min(newsubset$L1)+3
  bet <- max(newsubset$L2)-3
  
  
  gmodel2 <- with (newsubset,
                   grotag(L1=L1, L2=L2, T1=T1, T2=T2,alpha=alp,beta= bet,
                          design=list(nu=0,m=1,p=0,sea=0),
                          stvalue=list(sigma=0.9,nu=0.4,m=-1,p=0,u=0,w=0.1),
                          upper=list(sigma=5,nu=1,m=2,p=0.5,u=0,w=0.9),
                          lower=list(sigma=0,nu=0,m=-2,p=0,u=0,w=0),control=list(maxit=1e4)))
  
  Linf <- gmodel2$VBparms$Estimate[2]  # Linf
   # Linf <- (- 1) * Linf
  Linf <- round(Linf, digits = 1)
  Linf
  
  K <- gmodel2$VBparms$Estimate[3]  # K
  K
  K <- round(K, digits = 3)
  
  
  # vbgf curve - Guaiamum data
  temp<-newsubset[c(newsubset$T2-newsubset$T1)>0,]
  growthTraject(K,Linf,lentag=temp$L1, lenrec=temp$L2,timelib=c(temp$T2-temp$T1), ylim = c(0,80), 
                main = "subset, n = 40")

  outp <- c(K,Linf)
    
  ; outp # K , Linf
  
}

# resample.fun.full(guaiamum1m, 130)
# for (i in 1:3)  {  results <- resample.fun.full(guaiamum1m, 42) }
#  with error suppression
# for (i in 1:10) {results <- tryCatch( resample.fun.full(guaiamum1m, 42) ) }
# for (i in 1:10) { tryCatch( resample.fun.full(guaiamum1m, 42) ) }

### BOOTSTRAP OK ### 
## generate output matrix

N.rep <- 1000 # N of repeated resamplings  40 = sample size (1/2 of 80)

outp.mat1 <- data.frame(c(1:(N.rep)),rep.int(0,(N.rep)))
names(outp.mat1) <- c("K", "Linf")
summary(outp.mat1)

for(line.no in 1:(N.rep)){
  
  tryCatch({ results <-  resample.fun.full(guaiamum1m, 40)

    outp.mat1$K[line.no] <- results[1]
    outp.mat1$Linf[line.no] <- results[2]
    
  }, error=function(e){})
    
}

View(outp.mat1)

# write results file ----------------
# first remove 0 values for failed fit
outp.mat.clean <- subset(outp.mat1, Linf != 0)

View(outp.mat.clean)

outp <- data.frame(outp.mat.clean)

write.csv(outp, file = "1000Stoutput_run2_C1.csv")
#Do this for all runs (9 in total), and paste together
#Plot at end after putting all together. 

##
#To clean memory
gc()  
gc(reset=T)
rm(list = ls())

outp.mat.clean = read.csv("1000StoutputB1_9_from_Margit.csv", header =TRUE)

# mean and 95% confidence interval for K
mean(outp.mat.clean$K) # mean
median(outp.mat.clean$K) # median
quantile(outp.mat.clean$K, c(0.025, 0.975))  # 95% confidence interval
nrow(outp.mat.clean) # N successful runs

# mean and 95% confidence interval for Linf
mean(outp.mat.clean$Linf)    # mean
median(outp.mat.clean$Linf) # median
quantile(outp.mat.clean$Linf, c(0.025, 0.975)) # 95% confidence interval

hist(outp.mat.clean$Linf)
hist(outp.mat.clean$Linf, xlim = c(30, 200), ylim = c(0, 50), breaks = seq(0, 800000, by = 0.5) )
hist(outp.mat.clean$K)
hist(outp.mat.clean$K, xlim = c(0, 0.8), breaks = seq(0,1.2, by = 0.005) )

boxplot(outp.mat.clean$Linf)
boxplot(outp.mat.clean$Linf, ylim = c(40, 200), ylab= "Linf (cm)")

library(beanplot)
par <- par(mfrow = c(1,1), mar = c(3,3,2,1), oma = c(2,0,0,0))
beanplot(outp.mat.clean$Linf, beanlines = "quantile", ylab = "log Linf (cm)")
beanplot(outp.mat.clean$K, beanlines = "quantile", ylab = "K (/year)")

# xy plots
plot(outp.mat.clean$K , outp.mat.clean$Linf)
plot(outp.mat.clean$K , outp.mat.clean$Linf, ylim = c(40,200), xlab = "K (/year)", ylab="Linf (cm)")

