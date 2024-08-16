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
# Steenbras (Lithognatus aureti) growth ----------
#  -----------------
# v31 --------------
# July 19, 2024 -------------
# reading in and analyzing posteriors 
# Ralf Schwamborn: ralf.schwamborn@ufpe.br
# seed = Sys.Time!
# MA 3 vs MA 5
# Linfmin50 vs Linfmin30
# Tag-Recapture
# Otoliths: Northern vs Southern Stock

# Data: 
# Fish lengths (fork length) from 2004-08-25 to 2009-12-30 
# (approx. five and a half years) 
# Length range: 12 cm to 79 cm


# Contents:     ---------------------------------------------------------

# I. Chapter I: read in the data and make fishboot plots ---------------------
# I. A: MA3 (Linf 50 to 200cm) 
# I. B: MA5 (Linf 50 to 200cm)
# I. C. MA5 (Linf 30 to 200cm) 
# I. D. MA3 (Linf 30 to 200cm) 
# I. E. Tag-recapture
# I. F. Otoliths


# II. Chapter II: Comparisons and tests ------------------------------
# II.1 comparing and testing differences in medians
# II.1 comparing and testing differences in confidence interval width 


# III. Chapter III: Plots, customized boxplots------------------------
# Plots of 95 % confidence intervals (customized boxplots)

# IV. Chapter IV: S vs N stock, males vs females vs herm. (otoliths only) ------------------------
# S vs N stock (otoliths only)
# males vs females vs hermaphrodites (otoliths only)


# V. Chapter V: simple LFD plots for paper, with example "optimum curve" ----------------------

# VI. Chapter VI: simple T&R plots for paper,  --------------------------------------
#      with "curve swarms", increments  and "optimum curve" --------------------


# VII. Chapter VII:  LAA (otolith reading data) plots for paper,  --------------------------------------
#      with "curve swarms",  and "optimum curve" --------------------
#      "optimum curve": Median and Max. density  
#       Northern  vs Southern stock 
#       colors by sex (J., M., H., F.) 

# VIII. Chapter VIII: simple LAA (otolith reading data) plots for paper,  --------------------------------------
#      without "curve swarms", ages by sex and "optimum curve" --------------------
#       Northern  vs Southern stock 
#       Direct testing of medion length,  Northern  vs Southern stock
#       Permutation test Northern  vs Southern stock, median length at age
#       Simple boxplots
#       Simple points plot (symbols and cols by region)

# IX. Recruitment analysis (numbers of small indiv.) -------------
#          Numbers of small fish, and description of "first peak"


###########
##########
#########

# 0. Initialization ----------------------------------------------------------

# clean memory
gc()  
gc(reset=T)
rm(list = ls())

opar <- par()  


# required packages =====
library(parallel)
library(TropFishR)
library(fishboot)
library(ks)


# I. Chapter I: read in the data and make fishboot plots ---------------------
#  Reading in and analyzing posteriors (MA 3 vs MA 5, Linfmin30 vs Linfmin50) ---------------


# I.A.  MA 3 (Linf 50 to 200 cm) ---------------------------

# A1 set wd -------------------

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
getwd()

#path2 = path +'/' + "LFA_S_stock_New_analyses_Sys_time_results"


# A2 read data, MA3, Linfmin = 50 cm -------------------

# 3 new analyses (seed = Sys.time
res_seed_Sys_time_exp8_MA3_n200_Linfmin50_8a <- read.csv("res_seed_Sys_time_exp8_MA3_n200_Linfmin50_8a.csv", sep=";")
res_seed_Sys_time_exp8_MA3_n300_Linfmin50_8b <- read.csv("res_seed_Sys_time_exp8_MA3_n300_Linfmin50_8b.csv", sep=";")
res_seed_Sys_time_exp8A_MA3_n300_Linfmin50_8c <- read.csv("res_seed_Sys_time_exp8A_MA3_n300_Linfmin50_8c.csv", sep=";")

# 1 old analysis (OK)
res_EXP_S_slow_n200_1c <- read.csv("res_EXP_S_slow_n200_1c.csv", sep="")


# check for identical results ------
head(res_seed_Sys_time_exp8_MA3_n200_Linfmin50_8a) 
head(res_seed_Sys_time_exp8_MA3_n300_Linfmin50_8b )
head(res_seed_Sys_time_exp8A_MA3_n300_Linfmin50_8c )
head(res_EXP_S_slow_n200_1c) 

# OK< not identical
 
# check for median Link ad K (should be similar)
 summary(res_seed_Sys_time_exp8_MA3_n200_Linfmin50_8a) # median : Linf 91.31   Median K:0.13131
 summary(res_seed_Sys_time_exp8_MA3_n300_Linfmin50_8b )#Median : 95.28   Median K:0.14143
 summary(res_seed_Sys_time_exp8A_MA3_n300_Linfmin50_8c )# Median : 86.50   Median K:0.09941
 summary(res_EXP_S_slow_n200_1c) # Median : 86.66   Median :0.10162

               
# code to check for identical results --------
table(res_seed_Sys_time_exp8_MA3_n300_Linfmin50_8b == res_seed_Sys_time_exp8A_MA3_n300_Linfmin50_8c, useNA = 'ifany')
# not identical!


# A 3 merge the data, MA3, Linfimin 50  ------------------------

 expA_MA3_Linfmin50_Sys.Time.res_1000 <- rbind(res_EXP_S_slow_n200_1c,
                res_seed_Sys_time_exp8_MA3_n200_Linfmin50_8a,
                res_seed_Sys_time_exp8_MA3_n300_Linfmin50_8b,
                res_seed_Sys_time_exp8A_MA3_n300_Linfmin50_8c)
 summary(expA_MA3_Linfmin50_Sys.Time.res_1000)
 length(expA_MA3_Linfmin50_Sys.Time.res_1000$Linf) # 1000 runs, OK
 
                                          
 
View(MA3_res_EXP_S_slow_exp1_ALLSizes) # 1000 rows, OK

# A 4 Save and Analyse the results, MA3 ----------------------------------

write.table(expA_MA3_Linfmin50_Sys.Time.res_1000, 
            "expA_MA3_Linfmin50_Sys.Time.res_1000.csv", sep = ";") 

res_MA3_Linf50to200 <- expA_MA3_Linfmin50_Sys.Time.res_1000

summary(res_MA3_Linf50to200)

# > summary(res_MA3_Linf50to200)
# Linf              K              t_anchor            C                  ts          
# Min.   : 51.66   Min.   :0.02409   Min.   :0.0601   Min.   :0.001018   Min.   :0.008406  
# 1st Qu.: 78.33   1st Qu.:0.07447   1st Qu.:0.3816   1st Qu.:0.405840   1st Qu.:0.348751  
# Median : 89.57   Median :0.11628   Median :0.4974   Median :0.528622   Median :0.466486  
# Mean   : 94.48   Mean   :0.14573   Mean   :0.5074   Mean   :0.534049   Mean   :0.469838  
# 3rd Qu.:105.68   3rd Qu.:0.18815   3rd Qu.:0.6362   3rd Qu.:0.663533   3rd Qu.:0.584086  
# Max.   :181.27   Max.   :1.02473   Max.   :0.9911   Max.   :0.997138   Max.   :0.999698  
# phiL      
# Min.   :2.672  
# 1st Qu.:2.882  
# Median :2.978  
# Mean   :3.002  
# 3rd Qu.:3.110  
# Max.   :3.598 

#  comment: Ctrl+Shift+C

# 95% CI and CIW, CIW% ---------------------------------

#define functions for lower and upper 95% CI limits ----------
quantile_lowlim <- function(x) {quantile(x,  0.025)}
quantile_upperlim <- function(x) {quantile(x,  0.975)}

#lower 95%CI limit
quantile(res_MA3_Linf50to200$Linf, 0.025)# OK, 63.31596 cm = NEW,  68.6 cm, Old
quantile_lowlim(res_MA3_Linf50to200$Linf)# # OK, 63.31596 cm = NEW , tested, OK, 68.6 cm OK

# apply functions by columns, for for lower and upper 95% CI limits -------------

lapply(res_MA3_Linf50to200, quantile_lowlim)
lapply(res_MA3_Linf50to200, quantile_upperlim)

# 
#
# $Linf
# 97.5% 
# 151.9674 
# 
# $K
# 97.5% 
# 0.4438945 
# 
# $t_anchor
# 97.5% 
# 0.8659509 
# 
# $C
# 97.5% 
# 0.8736716 
# 
# $ts
# 97.5% 
# 0.8790932 
# 
# $phiL
# 97.5% 
# 3.31553 


# make a table with median and 95%CI --------------

vec1 <- unlist(lapply(res_MA3_Linf50to200,median))
vec2 <- unlist(lapply(res_MA3_Linf50to200, quantile_lowlim))
vec3 <- unlist(lapply(res_MA3_Linf50to200, quantile_upperlim))

class(vec2)

table.1a.MA3_Linfmin50.res_median_CI <- data.frame( median = vec1,
                                                   lowCI = vec2,
                                                   upperCI = vec3)

table.1a.MA3_Linfmin50.res_median_CI
View(table1a.MA3_Linfmin50.res_median_CI)

# additional columns (text95%CI, CIW, CIW%,text95%CI)

attach(table.1a.MA3_Linfmin50.res_median_CI)

table.1a.MA3_Linfmin50.res_median_CI$CIW <-  upperCI -lowCI      

attach(table.1a.MA3_Linfmin50.res_median_CI)

table.1a.MA3_Linfmin50.res_median_CI$CIW_perc <-   ( CIW /  median  ) * 100        

attach(table.1a.MA3_Linfmin50.res_median_CI)


# round values (3 digits)                                                       
table.1a.MA3_Linfmin50.res_median_CI <- round (table.1a.MA3_Linfmin50.res_median_CI[,1:5], 3)                                                    
                                                       
#  create a new column `text95CI` ----------------                                                       
                                                    
# columns to paste together for text
    cols2 <- c("lowCI", "upperCI")
                                                       
# create a new column `text95CI` with the three columns collapsed together (character string)
table.1a.MA3_Linfmin50.res_median_CI$text95CI <- apply(table.1a.MA3_Linfmin50.res_median_CI[ , cols2 ] , 1 , paste , collapse = "_to_" )
                                                       
  
  table.1a.MA3_Linfmin50.res_median_CI

# export (write to hard drive)   ---------------------
  
  path = dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(path)

  write.table(table.1a.MA3_Linfmin50.res_median_CI, file = "table.1a.MA3_Linfmin50.res_median_CI.csv")
  
  write.table(res_MA3_Linf50to200, file = "res_MA3_Linfmin50_Sys.TimeOK.csv")

  
# create  an empty fishboot "res" object (for plotting) ----------------------------
  library(TropFishR)
  
  library(fishboot)
  
  data(alba)

# ?ELEFAN_SA_boot

# settings (these settings may not be optimal - for demo only)
MA <- 7
init_par <- NULL
low_par <- list(Linf = 9, K = 0.3, t_anchor = 0)
up_par <- list(Linf = 11, K = 1, t_anchor = 1)
SA_time <- 1
SA_temp <- 1e5
nresamp <- 2


## parallel version
library(parallel)
t1 <- Sys.time()
res <- ELEFAN_SA_boot(lfq=alba, MA = MA, seasonalised = FALSE,
                      init_par = init_par, up_par = up_par, low_par = low_par,
                      SA_time = SA_time, SA_temp = SA_temp,
                      nresamp = nresamp, parallel = TRUE, no_cores = detectCores()-1,
                      seed = 1
)
t2 <- Sys.time()
t2 - t1
res

# insert  MA 3 Steenbras data into the res object ------------

res_ok <- res

res_ok$bootRaw <- res_MA3_Linf50to200


# plotting with fishboot (MA3)  -----------------------


# par(opar)


# univariate density plot of bootstrapped# pars

res <- res_ok
univariate_density(res)

# par(opar)


# Linf / K scatterhist GA
LinfK_scatterhist(res, Linf.breaks=50, K.breaks=50, 
                  phi.contour = TRUE, phi.contour.lty = 1, phi.contour.col = 8,
                  pt.col = rgb(0,0,0,0.1), pt.pch = "."
)


# par(opar)


# VBGF by time growth curve plot
CIinfo <- vbgfCI_time(
  res = res, 
  agemax = 12, CI = 95,
  perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
  ci.col = 1, ci.lty = 2, ci.lwd = 1, 
  maxd.col = 1, maxd.lty = 1, maxd.lwd = 2
)


# par(opar)


###############
##########
#######
##
# I.B. MA 5 (Linf 50 to 200 cm) -----------------------------------

# I. B 1:  set wd -------------------

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# B2 read data, MA5 -------------------
# NEW
# MA 5, Linfmin 50
# 3 new analyses (seed = Sys.time)
res_seed_Sys_time_exp8_MA5_n200_Linfmin50_8a <- read.csv("res_seed_Sys_time_exp8_MA5_n200_Linfmin50_8a.csv", sep=";")
res_seed_Sys_time_exp8_MA5_n300_Linfmin50_8b <- read.csv("res_seed_Sys_time_exp8_MA5_n300_Linfmin50_8b.csv", sep=";")
res_seed_Sys_time_expB9_MA5_n200_Linfmin50_9a <- read.csv("res_seed_Sys_time_expB9_MA5_n200_Linfmin50_9a.csv", sep=";")

# 1 old analysis (seed=1)
res_EXP_S_MA5_slow_n300_2c <- read.csv("res_EXP_S_MA5_slow_n300_2c.csv", sep="")


# check for identical results ------
head(res_seed_Sys_time_exp8_MA5_n200_Linfmin50_8a)
head(res_seed_Sys_time_exp8_MA5_n300_Linfmin50_8b)
head(res_seed_Sys_time_expB9_MA5_n200_Linfmin50_9a)
head(res_EXP_S_MA5_slow_n300_2c)


# B 3 merge the data, MA5, Linfmin 50 ------------------------

expB_MA5_Linfmin50_Sys.Time.res_1000 <- rbind(res_EXP_S_MA5_slow_n300_2c,
                                        res_seed_Sys_time_exp8_MA5_n200_Linfmin50_8a,
                                        res_seed_Sys_time_exp8_MA5_n300_Linfmin50_8b,
                                        res_seed_Sys_time_expB9_MA5_n200_Linfmin50_9a)
                                                   

View(MA5_res_EXP_S_slow_exp1_ALLSizes) # 1000 rows, OK

# B4 Save and Analyse the results, MA5 ----------------------------------

write.table(expB_MA5_Linfmin50_Sys.Time.res_1000, 
            "expB_MA5_Linfmin50_Sys.Time.res_1000.csv", sep = ";") 


summary(expB_MA5_Linfmin50_Sys.Time.res_1000)


# MA5, Linfmin 50 , NEW (seed = Sys.time, OK)

# > summary(expB_MA5_Linfmin50_Sys.Time.res_1000)
# Linf              K              t_anchor             C                  ts         
# Min.   : 51.66   Min.   :0.02357   Min.   :0.01307   Min.   :0.001018   Min.   :0.01657  
# 1st Qu.: 82.03   1st Qu.:0.08114   1st Qu.:0.35959   1st Qu.:0.427388   1st Qu.:0.35767  
# Median : 95.30   Median :0.13165   Median :0.48645   Median :0.546923   Median :0.45543  
# Mean   : 99.02   Mean   :0.16417   Mean   :0.49489   Mean   :0.548926   Mean   :0.46164  
# 3rd Qu.:112.56   3rd Qu.:0.20588   3rd Qu.:0.63494   3rd Qu.:0.676206   3rd Qu.:0.55723  
# Max.   :189.48   Max.   :1.02473   Max.   :0.96448   Max.   :0.997138   Max.   :0.99970  
# phiL      
# Min.   :2.740  
# 1st Qu.:2.966  
# Median :3.084  
# Mean   :3.088  
# 3rd Qu.:3.214  
# Max.   :3.598  



#  comment: Ctrl+Shift+C


# 95% CI and CIW, CIW% ---------------------------------

#define functions for lower and upper 95% CI limits ----------
quantile_lowlim <- function(x) {quantile(x,  0.025)}
quantile_upperlim <- function(x) {quantile(x,  0.975)}

head(expB_MA5_Linfmin50_Sys.Time.res_1000)

#lower 95%CI limit
#quantile(res_MA3_Linf50to200$Linf, 0.025)# OK, 68.6 cm, OK
#quantile_lowlim(res_MA3_Linf50to200$Linf)# tested, OK, 68.6 cm OK

# apply functions by columns, for for lower and upper 95% CI limits -------------

lapply(expB_MA5_Linfmin50_Sys.Time.res_1000, quantile_lowlim)
lapply(expB_MA5_Linfmin50_Sys.Time.res_1000, quantile_upperlim)

# make a table with median and 95%CI --------------

vec1b <- unlist(lapply(expB_MA5_Linfmin50_Sys.Time.res_1000,median))
vec2b <- unlist(lapply(expB_MA5_Linfmin50_Sys.Time.res_1000, quantile_lowlim))
vec3b <- unlist(lapply(expB_MA5_Linfmin50_Sys.Time.res_1000, quantile_upperlim))

class(vec2b)

table.1cOK.MA5_Linfmin50.res_median_CI <- data.frame( median = vec1b,
                                                    lowCI = vec2b,
                                                    upperCI = vec3b)

table.1b.MA5_Linfmin50.res_median_CI <- table.1cOK.MA5_Linfmin50.res_median_CI
# View(table1a.MA3_Linfmin50.res_median_CI)

# additional columns (text95%CI, CIW, CIW%,text95%CI)

attach(table.1b.MA5_Linfmin50.res_median_CI)

table.1b.MA5_Linfmin50.res_median_CI$CIW <-  upperCI -lowCI      

attach(table.1b.MA5_Linfmin50.res_median_CI)

table.1b.MA5_Linfmin50.res_median_CI$CIW_perc <-   ( CIW /  median  ) * 100        

attach(table.1b.MA5_Linfmin50.res_median_CI)


# round values (3 digits)                                                       
table.1b.MA5_Linfmin50.res_median_CI <- round (table.1b.MA5_Linfmin50.res_median_CI[,1:5], 3)                                                    

#  create a new column `text95CI` ----------------                                                       

# columns to paste together for text
cols2 <- c("lowCI", "upperCI"  )

# create a new column `text95CI` with the three columns collapsed together (character string)
table.1b.MA5_Linfmin50.res_median_CI$text95CI <- apply( table.1b.MA5_Linfmin50.res_median_CI[ , cols2 ] , 1 , paste , collapse = "_to_" )

table.1b.MA5_Linfmin50.res_median_CI

# export (write to hard drive)   ---------------------

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

write.table(table.1b.MA5_Linfmin50.res_median_CI, file = "table.1b.MA5_Linfmin50.res_median_CI.csv")

write.table(expB_MA5_Linfmin50_Sys.Time.res_1000, file = "expB_MA5_Linfmin50_Sys.Time.res_1000.csv")

table.1b.MA5_Linfmin50.res_median_CI

# create  an empty fishboot "res" object ----------------------------

data(alba)

# ?ELEFAN_SA_boot

# settings (these settings may not be optimal - for demo only)
MA <- 7
init_par <- NULL
low_par <- list(Linf = 9, K = 0.6, t_anchor = 0)
up_par <- list(Linf = 12, K = 1, t_anchor = 1)
SA_time <- 1
SA_temp <- 1e5
nresamp <- 2


## parallel version
library(parallel)
t1 <- Sys.time()
res <- ELEFAN_SA_boot(lfq=alba, MA = MA, seasonalised = FALSE,
                      init_par = init_par, up_par = up_par, low_par = low_par,
                      SA_time = SA_time, SA_temp = SA_temp,
                      nresamp = nresamp, parallel = TRUE, no_cores = detectCores()-1,
                      seed = 1
)
t2 <- Sys.time()
t2 - t1
res

# insert  MA 5 Steenbras data into the res object ------------

res_tempMA5_Linfmin50 <- res

res_tempMA5_Linfmin50$bootRaw <- expB_MA5_Linfmin50_Sys.Time.res_1000

# plotting with fishboot (MA5)  -----------------------

median(res_tempMA5_Linfmin50$bootRaw$Linf)
median(res_tempMA5_Linfmin50$bootRaw$K)
median(res_tempMA5_Linfmin50$bootRaw$C)
median(res_tempMA5_Linfmin50$bootRaw$phiL)



opar <- par()

# univariate density plot of bootstrapped# 

res <- res_tempMA5_Linfmin50


# par(opar)



univariate_density(res)


# par(opar)



# Linf / K scatterhist GA
LinfK_scatterhist(res, Linf.breaks=50, K.breaks=50, 
                  phi.contour = TRUE, phi.contour.lty = 1, phi.contour.col = 8,
                  pt.col = rgb(0,0,0,0.1), pt.pch = "."
)


# par(opar)



# VBGF by time growth curve plot
CIinfo <- vbgfCI_time(
  res = res, 
  agemax = 12, CI = 95,
  perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
  ci.col = 1, ci.lty = 2, ci.lwd = 1, 
  maxd.col = 1, maxd.lty = 1, maxd.lwd = 2
)


# par(opar)

 par(opar)

 
 
 

###############
##########
#######
##
# I. C. MA5, Linf 30 to 200 cm --------------------------------
## objective: compare Linf_min 30 vs Linf_min 50 -------------------
# Linf_min 30, MA5


# I.C. Linf 30 to 200 cm (MA 5) -----------------------------------

# I. C 1:  set wd -------------------
 
path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
 
# C 2: read data; MA5, Linfmin 30,  Linf 30 to 200 cm (MA 5) -------------------

# Old (seed = 1)
res_EXP3_S_MA5_slowLinf30_200_n500_3b <- read.csv("res_EXP3_S_MA5_slowLinf30_200_n500_3b.csv", sep="")

summary(res_EXP3_S_MA5_slowLinf30_200_n500_3b)

# New (seed = Sys.time), MA5, Linfmin 30 cm

res_seed_Sys_time_expC10_MA5_n200_Linfmin30_10 <- read.csv("res_seed_Sys_time_expC10_MA5_n200_Linfmin30_10.csv", sep = ";")
res_seed_Sys_time_exp6_MA5_n300_Linfmin30_6b  <- read.csv("res_seed_Sys_time_exp6_MA5_n300_Linfmin30_6b.csv", sep = ";" )

# compare (identical results?)

head(res_EXP3_S_MA5_slowLinf30_200_n500_3b)
head( res_seed_Sys_time_exp6_MA5_n300_Linfmin30_6b)
head(res_seed_Sys_time_expC10_MA5_n200_Linfmin30_10)


summary(res_EXP3_S_MA5_slowLinf30_200_n500_3b)
summary( res_seed_Sys_time_exp6_MA5_n300_Linfmin30_6b)
summary(res_seed_Sys_time_expC10_MA5_n200_Linfmin30_10)



# C 3 merge de the data, MA5 Linfmin 30 c ------------------------

expC_MA5_Linfmin30_Sys.Time.res_1000 <- rbind(res_EXP3_S_MA5_slowLinf30_200_n500_3b,
                        res_seed_Sys_time_exp6_MA5_n300_Linfmin30_6b,
                        res_seed_Sys_time_expC10_MA5_n200_Linfmin30_10)

# View(expC_MA5_Linfmin30_Sys.Time.res_1000) # 1000 rows, OK

head(expC_MA5_Linfmin30_Sys.Time.res_1000)

summary(expC_MA5_Linfmin30_Sys.Time.res_1000)

length(expC_MA5_Linfmin30_Sys.Time.res_1000$Linf)# 1000 rows, OK

#  Save and Analyse the results, MA5 ----------------------------------

write.table(expC_MA5_Linfmin30_Sys.Time.res_1000, 
            "expC_MA5_Linfmin30_Sys.Time.res_1000.csv", sep = ";") 


summary(expC_MA5_Linfmin30_Sys.Time.res_1000)

# # NEW (MA5, Linfmin30cm)
# # 
# > summary(expC_MA5_Linfmin30_Sys.Time.res_1000)
# Linf              K             t_anchor             C                 ts         
# Min.   : 39.98   Min.   :0.0255   Min.   :0.02715   Min.   :0.07716   Min.   :0.01239  
# 1st Qu.: 44.69   1st Qu.:0.1446   1st Qu.:0.40795   1st Qu.:0.44803   1st Qu.:0.36798  
# Median : 51.34   Median :0.4880   Median :0.51203   Median :0.55343   Median :0.45129  
# Mean   : 69.23   Mean   :0.4610   Mean   :0.51110   Mean   :0.55920   Mean   :0.45492  
# 3rd Qu.: 90.39   3rd Qu.:0.7396   3rd Qu.:0.61911   3rd Qu.:0.68302   3rd Qu.:0.53802  
# Max.   :187.63   Max.   :1.3895   Max.   :0.98189   Max.   :0.99591   Max.   :0.99985  
# phiL      
# Min.   :2.748  
# 1st Qu.:3.044  
# Median :3.147  
# Mean   :3.118  
# 3rd Qu.:3.211  
# Max.   :3.382  


# Comment: Ctrl+Shift+C

# 95% CI and CIW, CIW% ---------------------------------

#define functions for lower and upper 95% CI limits ----------
quantile_lowlim <- function(x) {quantile(x,  0.025)}
quantile_upperlim <- function(x) {quantile(x,  0.975)}

head(expC_MA5_Linfmin30_Sys.Time.res_1000)

#lower 95%CI limit
quantile(expC_MA5_Linfmin30_Sys.Time.res_1000$Linf, 0.025)# OK, 41.8 cm, OK
quantile_lowlim(expC_MA5_Linfmin30_Sys.Time.res_1000$Linf)# tested, OK, 41.8 cm OK

# apply functions by columns, for  lower and upper 95% CI limits -------------

lapply(expC_MA5_Linfmin30_Sys.Time.res_1000, quantile_lowlim)

lapply(expC_MA5_Linfmin30_Sys.Time.res_1000, quantile_upperlim)

# make a table with median and 95%CI --------------

vec1c <- unlist(  lapply(expC_MA5_Linfmin30_Sys.Time.res_1000,median))
vec2c <- unlist(lapply(expC_MA5_Linfmin30_Sys.Time.res_1000, quantile_lowlim))
vec3c <- unlist(lapply(expC_MA5_Linfmin30_Sys.Time.res_1000, quantile_upperlim))

class(vec2c)

table.1c.MA5_Linfmin30.res_median_CI <- data.frame(median = vec1c,
                                                   lowCI = vec2c,
                                                   upperCI = vec3c)

table.1c.MA5_Linfmin30.res_median_CI
# View(table1a.MA3_Linfmin50.res_median_CI)

# additional columns (text95%CI, CIW, CIW%,text95%CI)

attach(table.1c.MA5_Linfmin30.res_median_CI)

table.1c.MA5_Linfmin30.res_median_CI$CIW <-  upperCI -lowCI      

attach(table.1c.MA5_Linfmin30.res_median_CI)

table.1c.MA5_Linfmin30.res_median_CI$CIW_perc <-   (CIW / median) * 100        

attach(table.1c.MA5_Linfmin30.res_median_CI)


# round values (3 digits)                                                       
table.1c.MA5_Linfmin30.res_median_CI <- round (table.1c.MA5_Linfmin30.res_median_CI[,1:5], 3)                                                    

#  create a new column `text95CI` ----------------                                                       

# columns to paste together for text
cols2 <- c("lowCI", "upperCI"  )

# create a new column `text95CI` with the three columns collapsed together (character string)
table.1c.MA5_Linfmin30.res_median_CI$text95CI <- apply( table.1c.MA5_Linfmin30.res_median_CI[ , cols2 ] , 1 , paste , collapse = "_to_" )


table.1c.MA5_Linfmin30.res_median_CI

# export (write to hard drive)   ---------------------

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

write.table(table.1c.MA5_Linfmin30.res_median_CI, file = "table.1c.MA5_Linfmin30.res_median_CI.csv")

write.table(expC_MA5_Linfmin30_Sys.Time.res_1000, file = "expC_MA5_Linfmin30_Sys.Time.res_1000.csv")



# create  an empty fishboot "res" object ----------------------------

data(alba)

# ?ELEFAN_SA_boot

# settings (these settings may not be optimal - for demo only)
MA <- 7
init_par <- NULL
low_par <- list(Linf = 9, K = 0.6, t_anchor = 0)
up_par <- list(Linf = 12, K = 1, t_anchor = 1)
SA_time <- 1
SA_temp <- 1e5
nresamp <- 2


## parallel version
library(parallel)
t1 <- Sys.time()
res <- ELEFAN_SA_boot(lfq=alba, MA = MA, seasonalised = FALSE,
                      init_par = init_par, up_par = up_par, low_par = low_par,
                      SA_time = SA_time, SA_temp = SA_temp,
                      nresamp = nresamp, parallel = TRUE, no_cores = detectCores()-1,
                      seed = 1
)
t2 <- Sys.time()
t2 - t1
res

# insert  MA  = 5, Linfmin = 30cm  Steenbras data into the res object ------------

res$bootRaw <- expC_MA5_Linfmin30_Sys.Time.res_1000


median (res$bootRaw$Linf)

median (res$bootRaw$K)


# plotting with fishboot (MA5)  -----------------------


# univariate density plot of bootstrapped# pars


# par(opar)




univariate_density(res)


# bimodal distribution! males and females!

# par(opar)



#Result (MA = 5 , Linfmin = 30 cm):
# huge uncertainty, 
# max density (most likely optimum) for Linf  at 44.40
# similar to tag-recap (most likely optimum) for Linf 
# (growth of males dominates the shape of the posterior distrib,)
# Kimura plot with two groups of results (males and females?) 


# Linf / K scatterhist GA
LinfK_scatterhist(res, Linf.breaks=50, K.breaks=50, 
                  phi.contour = TRUE, phi.contour.lty = 1, phi.contour.col = 8,
                  pt.col = rgb(0,0,0,0.1), pt.pch = "."
)


par(opar)



# VBGF by time growth curve plot
CIinfo <- vbgfCI_time(
  res = res, 
  agemax = 12, CI = 95,
  perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
  ci.col = 1, ci.lty = 2, ci.lwd = 1, 
  maxd.col = 1, maxd.lty = 1, maxd.lwd = 2
)


par(opar)


# plot raw and restructured LFQ data bs = 2 ----------------

# 1.2. Load data =====

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

lfq2 <- read.csv("stbrasLF-trunc.csv")
# View(lfq2)

lfq2
names(lfq2)
lfq2$DATE <- as.Date(lfq2$DATE, format = "%d/%m/%Y") #Date format needs to be input like that
lfq2

lfq2new <- lfqCreate(data = lfq2, Lname = "Length_CM", Dname = "DATE", Fname = "Frequency")
plot(lfq2new, Fname = "catch")
lfq2new

lfq <- lfq2new
lfq
# preliminary (!) mostly likely optimum (T = test)
 KT = 0.09
 LinfT = 44.40
 tanchorT = 0.51
 CT = 0.48
 tsT = 0.5

# adjust bin size (on the histogram) #Bin size 2!
synLFQ7a <- lfqModify(lfq, bin_size = 2) #mid-sizes are already at 0.5

lfqbin <- lfqRestructure(synLFQ7a, MA = 5, addl.sqrt = TRUE) #MA = moving average sample size 

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "catch", date.axis = "modern")
plot(lfqbin, Fname = "rcounts", date.axis = "modern")
par(opar)

# plot with most likely growth curve, MA = 5, Linf_min = 30cm --------
# (preliminary) ---------
plot(lfqbin, Fname = "rcounts", date.axis = "modern")
lt <- lfqFitCurves(lfqbin, 
                   par = list(Linf=LinfT, K=KT, 
                              t_anchor=tanchorT, 
                              C= CT, ts=tsT),
                   draw = TRUE, col = "grey", lty = 1, lwd=1.5)



# plot with most likely growth curve, MA = 5, Linf_min = 30cm --------
# (preliminary) ---------
plot(lfqbin, Fname = "catch", date.axis = "modern")
lt <- lfqFitCurves(lfqbin, 
                   par = list(Linf=LinfT, K=KT, 
                              t_anchor=tanchorT, 
                              C= CT, ts=tsT),
                   draw = TRUE, col = "grey", lty = 1, lwd=1.5)


par(opar)



#######
#####
###
# I. D. MA3, Linf 30 to 200 cm ---------------------------------
## objective: compare Linf_min 30 vs Linf_min 50 -------------------

# I.D. Linf 30 to 200 cm (MA 3) -----------------------------------

# I. D 1:  set wd -------------------

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# D 2: read data, Linf 30 to 200 cm (MA 3) -------------------

res_seed_Sys_time_exp5_MA3_Linfmin30_5a <- read.csv("res_seed_Sys_time_exp5_MA3_Linfmin30_5a.csv", sep = ";")# 100 runs
res_seed_Sys_time_exp5_MA3_Linfmin30_5b <- read.csv("res_seed_Sys_time_exp5_MA3_Linfmin30_5b.csv", sep = ";")# 100 runs
res_seed_Sys_time_exp5_MA3_n200_Linfmin30_5c <- read.csv("res_seed_Sys_time_exp5_MA3_n200_Linfmin30_5c.csv", sep = ";") # 200 runs
res_seed_Sys_time_exp5_MA3_n300_Linfmin30_5d <- read.csv("res_seed_Sys_time_exp5_MA3_n300_Linfmin30_5d.csv", sep = ";") # 300 runs

res_seed_Sys_time_expD11_MA3_n200_Linfmin30_11 <- read.csv("res_seed_Sys_time_expD11_MA3_n200_Linfmin30_11.csv", sep = ";") # 300 runs


# OLD analysis (seed = 1)
res_EXP4_S_MA3_slowLinf30_200_n1000_4a <- read.csv("res_EXP4_S_MA3_slowLinf30_200_n1000_4a.csv", sep = " ") # 300 runs
res_EXP4_S_MA3_slowLinf30_200_n100_4b <- res_EXP4_S_MA3_slowLinf30_200_n1000_4a[1:100,]


# check for identical results (MA = 3, Linfmin = 30) 
head(res_seed_Sys_time_exp5_MA3_Linfmin30_5a)
head(res_seed_Sys_time_exp5_MA3_Linfmin30_5b)
head(res_seed_Sys_time_exp5_MA3_n200_Linfmin30_5c)
head(res_seed_Sys_time_exp5_MA3_n300_Linfmin30_5d)

head(res_seed_Sys_time_expD11_MA3_n200_Linfmin30_11)

head(res_EXP4_S_MA3_slowLinf30_200_n100_4b)


summary(res_seed_Sys_time_exp5_MA3_Linfmin30_5a )
summary(res_seed_Sys_time_exp5_MA3_Linfmin30_5b)
summary(res_seed_Sys_time_exp5_MA3_n200_Linfmin30_5c)
summary(res_seed_Sys_time_exp5_MA3_n300_Linfmin30_5d)
            
summary(res_seed_Sys_time_expD11_MA3_n200_Linfmin30_11)
                    
summary(res_EXP4_S_MA3_slowLinf30_200_n100_4b)

                    
# D 3 merge de the data, (MA = 3, Linfmin = 30)  ------------------------

 expD_MA3_Linfmin30_Sys.Time.res_1000 <- rbind(res_seed_Sys_time_exp5_MA3_Linfmin30_5a,
                                      res_seed_Sys_time_exp5_MA3_Linfmin30_5b,
                                      res_seed_Sys_time_exp5_MA3_n200_Linfmin30_5c,
                                      res_seed_Sys_time_exp5_MA3_n300_Linfmin30_5d,
                                      res_seed_Sys_time_expD11_MA3_n200_Linfmin30_11,
                                      res_EXP4_S_MA3_slowLinf30_200_n100_4b)
                    
        
View(expD_MA3_Linfmin30_Sys.Time.res_1000) # 1000 rows, OK

# 4 Save and Analyse the results, MA5 ----------------------------------

write.table(expD_MA3_Linfmin30_Sys.Time.res_1000, 
             "expD_MA3_Linfmin30_Sys.Time.res_1000.csv", sep = ";") 


summary(expD_MA3_Linfmin30_Sys.Time.res_1000)

# MA3_Linfmin30 (NEW, OK) 
# > summary(expD_MA3_Linfmin30_Sys.Time.res_1000)
# Linf              K              t_anchor              C                 ts          
# Min.   : 52.87   Min.   :0.02498   Min.   :0.008982   Min.   :0.02772   Min.   :0.001587  
# 1st Qu.: 74.92   1st Qu.:0.07161   1st Qu.:0.381160   1st Qu.:0.39366   1st Qu.:0.346291  
# Median : 83.47   Median :0.11472   Median :0.503464   Median :0.52450   Median :0.479403  
# Mean   : 88.67   Mean   :0.13794   Mean   :0.501707   Mean   :0.52263   Mean   :0.474885  
# 3rd Qu.: 96.17   3rd Qu.:0.17944   3rd Qu.:0.624604   3rd Qu.:0.65029   3rd Qu.:0.589547  
# Max.   :181.57   Max.   :0.88826   Max.   :0.984738   Max.   :0.97983   Max.   :0.982956  
# phiL      
# Min.   :2.667  
# 1st Qu.:2.837  
# Median :2.907  
# Mean   :2.931  
# 3rd Qu.:3.018  
# Max.   :3.486 


#  comment: Ctrl+Shift+C

# 95% CI and CIW, CIW% ---------------------------------

#define functions for lower and upper 95% CI limits ----------
quantile_lowlim <- function(x) {quantile(x,  0.025)}
quantile_upperlim <- function(x) {quantile(x,  0.975)}

head(expD_MA3_Linfmin30_Sys.Time.res_1000)

res_MA3_Linf30_200_n1000 <- expD_MA3_Linfmin30_Sys.Time.res_1000

#lower 95%CI limit
quantile(res_MA3_Linf30_200_n1000$Linf, 0.025)# OK, 41.8 cm, OK
quantile_lowlim(res_MA3_Linf30_200_n1000$Linf)# tested, OK, 41.8 cm OK

# apply functions by columns, for  lower and upper 95% CI limits -------------

lapply(res_MA3_Linf30_200_n1000, quantile_lowlim)

lapply(res_MA3_Linf30_200_n1000, quantile_upperlim)

# make a table with median and 95%CI --------------

vec1d <- unlist(  lapply(res_MA3_Linf30_200_n1000,median))
vec2d <- unlist(lapply(res_MA3_Linf30_200_n1000, quantile_lowlim))
vec3d <- unlist(lapply(res_MA3_Linf30_200_n1000, quantile_upperlim))

class(vec2d)

table.1d.MA3_Linfmin30.res_median_CI <- data.frame( median = vec1d,
                                                    lowCI = vec2d,
                                                    upperCI = vec3d)

table.1d.MA3_Linfmin30.res_median_CI
# View(table1a.MA3_Linfmin50.res_median_CI)

# additional columns (text95%CI, CIW, CIW%,text95%CI)

attach(table.1d.MA3_Linfmin30.res_median_CI)

table.1d.MA3_Linfmin30.res_median_CI$CIW <- upperCI -lowCI      

attach(table.1d.MA3_Linfmin30.res_median_CI)

table.1d.MA3_Linfmin30.res_median_CI$CIW_perc <- (CIW / median) * 100        

attach(table.1d.MA3_Linfmin30.res_median_CI)


# round values (3 digits)                                                       
table.1d.MA3_Linfmin30.res_median_CI <- round (table.1d.MA3_Linfmin30.res_median_CI[,1:5], 3)                                                    

#  create a new column `text95CI` ----------------                                                       

# columns to paste together for text
cols2 <- c("lowCI", "upperCI"  )

# create a new column `text95CI` with the three columns collapsed together (character string)
table.1d.MA3_Linfmin30.res_median_CI$text95CI <- apply( table.1d.MA3_Linfmin30.res_median_CI[ , cols2 ] , 1 , paste , collapse = "_to_" )


table.1d.MA3_Linfmin30.res_median_CI

# export (write to hard drive)   ---------------------

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

write.table(table.1d.MA3_Linfmin30.res_median_CI, file = "table.1d.MA3_Linfmin30.res_median_CI.csv")

write.table(res_MA3_Linf30_200_n1000, file = "res_MA3_Linf30_200_n1000.csv")


# create  an empty fishboot "res" object ----------------------------

data(alba)

# ?ELEFAN_SA_boot

# settings (these settings may not be optimal - for demo only)
MA <- 7
init_par <- NULL
low_par <- list(Linf = 9, K = 0.6, t_anchor = 0)
up_par <- list(Linf = 12, K = 1, t_anchor = 1)
SA_time <- 1
SA_temp <- 1e5
nresamp <- 2


## parallel version
library(parallel)
t1 <- Sys.time()
res <- ELEFAN_SA_boot(lfq=alba, MA = MA, seasonalised = FALSE,
                      init_par = init_par, up_par = up_par, low_par = low_par,
                      SA_time = SA_time, SA_temp = SA_temp,
                      nresamp = nresamp, parallel = TRUE, no_cores = detectCores()-1,
                      seed = 1
)
t2 <- Sys.time()
t2 - t1
res

# insert  MA 5 Steenbras data into the res object ------------

res$bootRaw <- expD_MA3_Linfmin30_Sys.Time.res_1000

# plotting with fishboot (MA5)  -----------------------


par(opar)



# VBGF by time growth curve plot
CIinfo <- vbgfCI_time(
  res = res, 
  agemax = 12, CI = 95,
  perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
  ci.col = 1, ci.lty = 2, ci.lwd = 1, 
  maxd.col = 1, maxd.lty = 1, maxd.lwd = 2
)


par(opar)


# univariate density plot of bootstrapped# pars


# par(opar)

res$bootRaw

univariate_density(res)

# bimodal distribution! males and females?

# par(opar)

# Result (Linfmin = 30 cm):
# huge uncertainty, 
# max density (most likely optimum) for Linf  at 44.65
# similar to tag-recap (most likely optimum) for Linf 
# (growth of males dominates the shape of the posterior distrib,)
# Kimura plot with two groups of results (males and females?) 

par(opar)

# Linf / K scatterhist GA
LinfK_scatterhist(res, Linf.breaks=50, K.breaks=50, 
                  phi.contour = TRUE, phi.contour.lty = 1, phi.contour.col = 8,
                  pt.col = rgb(0,0,0,0.1), pt.pch = "."
)



 par(opar)

 
 
###########
######
##
# I.E  T&R ,  Tag-recapture (TR1b data set) ---------------------------
# TR1b, 2020-2022 data were added to th original dataset
# and removed two fish -13 and - 7.5 growth
# N = 80

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

#   read data, -------------------

res_E_1000StoutputB1_9_from_Margit <- read.csv("1000StoutputB1_9_from_Margit.csv", sep=",")

res_E_Tag_recap1000 <- res_E_1000StoutputB1_9_from_Margit


# 95% CI and CIW, CIW% ---------------------------------

#define functions for lower and upper 95% CI limits ----------
quantile_lowlim <- function(x) {quantile(x,  0.025)}
quantile_upperlim <- function(x) {quantile(x,  0.975)}

head(res_E_Tag_recap1000)

attach(res_E_Tag_recap1000)

res_E_Tag_recap1000$Phi_L <- log10(res_E_Tag_recap1000$K) + (2 * log10(res_E_Tag_recap1000$Linf))
#0.103  68.1 2.679131
log10(0.103) + (2 * log10(68.1)) # OK, 2.679131


head(res_E_Tag_recap1000)# OK

# reorganize, Linf first
res_E_Tag_recap1000 <- res_E_Tag_recap1000[,c(2,1,3)]

head(res_E_Tag_recap1000)# OK


#lower 95%CI limit
quantile(res_E_Tag_recap1000$Linf, 0.025)# OK, 43.14  cm, OK
quantile_lowlim(res_E_Tag_recap1000$Linf)# tested, OK, 43.14  cm OK

# apply functions by columns, for for lower and upper 95% CI limits -------------

lapply(res_E_Tag_recap1000, quantile_lowlim)

lapply(res_E_Tag_recap1000, quantile_upperlim)

# make a table with median and 95%CI --------------

vec1e <- unlist(lapply(res_E_Tag_recap1000,median))
vec2e <- unlist(lapply(res_E_Tag_recap1000, quantile_lowlim))
vec3e <- unlist(lapply(res_E_Tag_recap1000, quantile_upperlim))

class(vec2e)

table.1e.Tag_recap1000 <- data.frame(median = vec1e,
                                     lowCI = vec2e,
                                     upperCI = vec3e)

table.1e.Tag_recap1000
# View(table.1e.Tag_recap1000)

# additional columns (text95%CI, CIW, CIW%,text95%CI)

attach(table.1e.Tag_recap1000)

table.1e.Tag_recap1000$CIW <-  upperCI -lowCI      

attach(table.1e.Tag_recap1000)

table.1e.Tag_recap1000$CIW_perc <-   (CIW / median) * 100        

attach(table.1e.Tag_recap1000)


# round values (3 digits)                                                       
table.1e.Tag_recap1000 <- round (table.1e.Tag_recap1000[,1:5], 3)                                                    

#  create a new column `text95CI` ----------------                                                       

# columns to paste together for text
cols2 <- c("lowCI", "upperCI"  )

# create a new column `text95CI` with the three columns collapsed together (character string)
table.1e.Tag_recap1000$text95CI <- apply( table.1e.Tag_recap1000[ , cols2 ] , 1 , paste , collapse = "_to_" )


table.1e.Tag_recap1000

# export (write to hard drive)   ---------------------

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

write.table(table.1e.Tag_recap1000, file = "table.1e.Tag_recap1000.csv")

write.table(res_E_Tag_recap1000, file = "res_E_Tag_recap1000.csv")

res_E_Tag_recap1000




res$bootRaw <- res_1000StoutputB1_9_from_Margit

# plotting with fishboot (T&R)  -----------------------

# univariate density plot of bootstrapped# pars


# par(opar)

res$bootRaw
res$bootRaw$t_anchor <- rep(0, length (res$bootRaw$Linf))



univariate_density(res)

# bimodal distribution! males and females?

# par(opar)



# Linf / K scatterhist GA
LinfK_scatterhist(res, Linf.breaks=50, K.breaks=50, 
                  phi.contour = TRUE, phi.contour.lty = 1, phi.contour.col = 8,
                  pt.col = rgb(0,0,0,0.1), pt.pch = ".", 
                  xlim = c(0,200))



# par(opar)



# VBGF by time growth curve plot
CIinfo <- vbgfCI_time(
  res = res, 
  agemax = 12, CI = 95,
  perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
  ci.col = 1, ci.lty = 2, ci.lwd = 1, 
  maxd.col = 1, maxd.lty = 1, maxd.lwd = 2
)


# I.F.a. Plot LFDs with otolith-derived growth curve -------------------

######################
#####################
###################
# FOR PAPER 
# LFDs & VBGF from OTOLITHs
# FOR PAPER - plot LFDs with otolith-derived growth curve
# plot raw and restructured LFQ data bs = 2 ----------------

# 1.2. Load data =====

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

lfq2 <- read.csv("stbrasLF-trunc.csv")

# View(lfq2)

lfq2
names(lfq2)
lfq2$DATE <- as.Date(lfq2$DATE, format = "%d/%m/%Y") #Date format needs to be input like that
lfq2

lfq2new <- lfqCreate(data = lfq2, Lname = "Length_CM", Dname = "DATE", Fname = "Frequency")
plot(lfq2new, Fname = "catch")
lfq2new

lfq <- lfq2new
lfq

# preliminary (!) mostly likely optimum (T = test)
# MA 5 , Linf min = 50 cm

KT = 0.10
LinfT = 87.97
tanchorT = 0.45
CT = 0.51
tsT = 0.49


# adjust bin size (on the histogram) #Bin size 2!
synLFQ7a <- lfqModify(lfq, bin_size = 2) #mid-sizes are already at 0.5

lfqbin <- lfqRestructure(synLFQ7a, MA = 5, addl.sqrt = TRUE) #MA = moving average sample size 
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "catch", date.axis = "modern")
plot(lfqbin, Fname = "rcounts", date.axis = "modern")
par(opar)

# plot with most likely growth curve, MA = 5, Linf_min = 30cm --------
# (preliminary) ---------
plot(lfqbin, Fname = "rcounts", date.axis = "modern")
lt <- lfqFitCurves(lfqbin, 
                   par = list(Linf=LinfT, K=KT, 
                              t_anchor=tanchorT, 
                              C= CT, ts=tsT),
                   draw = TRUE, col = "grey", lty = 1, lwd=1.5)



# plot with most likely growth curve, MA = 5, Linf_min = 30cm --------
# (preliminary) ---------
plot(lfqbin, Fname = "catch", date.axis = "modern")
lt <- lfqFitCurves(lfqbin, 
                   par = list(Linf=LinfT, K=KT, 
                              t_anchor=tanchorT, 
                              C= CT, ts=tsT),
                   draw = TRUE, col = "grey", lty = 1, lwd=1.5)


# plot with otolith-derived growth curve, 
# 
# (preliminary) ---------

# otolith-derived paranters K and Linf,  
KT = 0.066 # # otolith-derived
LinfT = 77.29 # otolith-derived
tanchorT = 0.45 # from LFA, MA5, Linfmin = 50
CT = 0.51 # # from LFA, MA5, Linfmin = 50
tsT = 0.49# from LFA, MA5, Linfmin = 50


plot(lfqbin, Fname = "catch", date.axis = "modern")
lt <- lfqFitCurves(lfqbin, 
                   par = list(Linf=77.29, K=0.066, 
                              t_anchor=tanchorT, 
                              C= CT, ts=tsT),
                   draw = TRUE, col = "grey", lty = 1, lwd=1.5)


plot(lfqbin, Fname = "catch", date.axis = "modern")
lt <- lfqFitCurves(lfqbin, 
                   par = list(Linf=77.29, K=0.066, 
                              t_anchor=tanchorT, 
                              C= 0, ts=0),
                   draw = TRUE, col = "grey", lty = 1, lwd=1.5)


# I.F.b. Analyse otolith data S stock (Length-at age, LAA) S stock

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

#   read data, -------------------

res.oto.S.all.1000_ini <- read.csv("oto-S-all-1000results.csv")

res.oto.S.all.1000 <- res.oto.S.all.1000_ini

#   View(res.oto.S.all.1000)

summary(res.oto.S.all.1000)
# 
# Linf              K                  t0         
# Min.   : 57.12   Min.   :0.002288   Min.   :-15.346  
# 1st Qu.: 69.85   1st Qu.:0.043391   1st Qu.: -7.699  
# Median : 77.29   Median :0.066345   Median : -5.712  
# Mean   : 99.90   Mean   :0.065436   Mean   : -6.289  
# 3rd Qu.: 93.57   3rd Qu.:0.087917   3rd Qu.: -4.417  
# Max.   :915.22   Max.   :0.170270   Max.   : -1.719  

# Comment: Ctrl + Shift + C


# 95% CI and CIW, CIW% ---------------------------------

#define functions for lower and upper 95% CI limits ----------
quantile_lowlim <- function(x) {quantile(x,  0.025)}
quantile_upperlim <- function(x) {quantile(x,  0.975)}

head(res.oto.S.all.1000)

attach(res.oto.S.all.1000)

res.oto.S.all.1000$Phi_L <- log10(res.oto.S.all.1000$K) + (2 * log10(res.oto.S.all.1000$Linf))
#0.103  68.1 2.679131
log10(0.103) + (2 * log10(68.1)) # OK, 2.679131


head(res.oto.S.all.1000)# OK

# reorganize, Linf first
res.oto.S.all.1000 <- res.oto.S.all.1000[,c(1,2,4)]

head(res.oto.S.all.1000)# OK


#lower 95%CI limit
quantile(res.oto.S.all.1000$Linf, 0.025)# OK, 63.08  cm, OK
quantile_lowlim(res.oto.S.all.1000$Linf)# tested, OK, 63.08  cm OK

# apply functions by columns, for for lower and upper 95% CI limits -------------

lapply(res.oto.S.all.1000, quantile_lowlim)

lapply(res.oto.S.all.1000, quantile_upperlim)

# make a table with median and 95%CI --------------

vec1f <- unlist(  lapply(res.oto.S.all.1000,median))
vec2f <- unlist(lapply(res.oto.S.all.1000, quantile_lowlim))
vec3f <- unlist(lapply(res.oto.S.all.1000, quantile_upperlim))

class(vec2f)

table.1f.S_oto_1000 <- data.frame(median = vec1f,
                                  lowCI = vec2f,
                                  upperCI = vec3f)

table.1f.S_oto_1000
# View(table.1f.S_oto_1000)

# additional columns (text95%CI, CIW, CIW%,text95%CI)

attach(table.1f.S_oto_1000)

table.1f.S_oto_1000$CIW <-  table.1f.S_oto_1000$upperCI -table.1f.S_oto_1000$lowCI      

attach(table.1f.S_oto_1000)

table.1f.S_oto_1000$CIW_perc <- (table.1f.S_oto_1000$CIW / table.1f.S_oto_1000$median) * 100        

attach(table.1f.S_oto_1000)


# round values (3 digits)                                                       
table.1f.S_oto_1000 <- round (table.1f.S_oto_1000[,1:5], 3)                                                    

# create a new column `text95CI` ----------------                                                       

# columns to paste together for text
cols2 <- c("lowCI", "upperCI"  )

# create a new column `text95CI` with the three columns collapsed together (character string)
table.1f.S_oto_1000$text95CI <- apply( table.1f.S_oto_1000[ , cols2 ] , 1 , paste , collapse = "_to_" )


table.1f.S_oto_1000

# export (write to hard drive)   ---------------------

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

write.table(table.1f.S_oto_1000, file = "table.1f.S_oto_1000.csv")

table.1f.S_oto_1000



# plotting with fishboot (otoliths)  -----------------------


res$bootRaw <- res.oto.S.all.1000_ini


# univariate density plot of bootstrapped# pars


# par(opar)

res$bootRaw




univariate_density(res)


 par(opar)



# Linf / K scatterhist GA

  LinfK_scatterhist(res, Linf.breaks=50, K.breaks=50, 
                   phi.contour = TRUE, phi.contour.lty = 1, phi.contour.col = 8,
                   pt.col = rgb(0,0,0,0.1), pt.pch = ".",
                   xlim = c(0,400) )
#  
#  
# LinfK_scatterhist(res, Linf.breaks=50, K.breaks=50, 
#                   phi.contour = TRUE, phi.contour.lty = 1, phi.contour.col = 8,
#                   pt.col = rgb(0,0,0,0.1), pt.pch = "."
#                   , xlim = c(0,200))
# 


# par(opar)

  res$bootRaw

# VBGF by time growth curve plot
  
  par(opar)
  
  res$bootRaw$t_anchor <-  res.oto.S.all.1000_ini$t0
 
  
CIinfo <- vbgfCI_time(
  res = res, 
  agemax = 12, CI = 95,
  perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
  ci.col = 1, ci.lty = 2, ci.lwd = 1, 
  maxd.col = 1, maxd.lty = 1, maxd.lwd = 2
)



# I.F.c. Analyse otolith data N stock (Length-at age, LAA) N stock

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

#   read data, -------------------


res.oto.N.all.1000_ini <- read.csv("oto-N-all-1000results.csv")


res.oto.N.all.1000 <- res.oto.N.all.1000_ini

#   View(res.oto.N.all.1000)

#summary(res.oto.S.all.1000)
# S Stock (LAA)
# Linf              K                  t0         
# Min.   : 57.12   Min.   :0.002288   Min.   :-15.346  
# 1st Qu.: 69.85   1st Qu.:0.043391   1st Qu.: -7.699  
# Median : 77.29   Median :0.066345   Median : -5.712  
# Mean   : 99.90   Mean   :0.065436   Mean   : -6.289  
# 3rd Qu.: 93.57   3rd Qu.:0.087917   3rd Qu.: -4.417  
# Max.   :915.22   Max.   :0.170270   Max.   : -1.719  

summary(res.oto.N.all.1000)
# N Stock (LAA)
# > summary(res.oto.N.all.1000)
# Linf              K                 t0         
# Min.   : 67.42   Min.   :0.04048   Min.   :-5.8495  
# 1st Qu.: 81.27   1st Qu.:0.07034   1st Qu.:-3.5163  
# Median : 85.64   Median :0.08062   Median :-3.0273  
# Mean   : 85.85   Mean   :0.08268   Mean   :-3.0597  
# 3rd Qu.: 89.85   3rd Qu.:0.09324   3rd Qu.:-2.5227  
# Max.   :111.59   Max.   :0.16346   Max.   :-0.9613  

# Comment: Ctrl + Shift + C


# 95% CI and CIW, CIW% ---------------------------------

#define functions for lower and upper 95% CI limits ----------
quantile_lowlim <- function(x) {quantile(x,  0.025)}
quantile_upperlim <- function(x) {quantile(x,  0.975)}

head(res.oto.N.all.1000)

attach(res.oto.N.all.1000)

res.oto.N.all.1000$Phi_L <- log10(res.oto.N.all.1000$K) + (2 * log10(res.oto.N.all.1000$Linf))
#0.103  68.1 2.679131
log10(0.103) + (2 * log10(68.1)) # OK, 2.679131


head(res.oto.N.all.1000)# OK

# reorganize, Linf first
res.oto.N.all.1000 <- res.oto.N.all.1000[,c(1,2,4)]

head(res.oto.N.all.1000)# OK


#lower 95%CI limit
quantile(res.oto.N.all.1000$Linf, 0.025)# OK, 63.08  cm, OK
quantile_lowlim(res.oto.N.all.1000$Linf)# tested, OK, 63.08  cm OK

# apply functions by columns, for for lower and upper 95% CI limits -------------

lapply(res.oto.N.all.1000, quantile_lowlim)

lapply(res.oto.N.all.1000, quantile_upperlim)

# make a table with median and 95%CI --------------


vec1f <- unlist(lapply(res.oto.N.all.1000,median))
vec2f <- unlist(lapply(res.oto.N.all.1000, quantile_lowlim))
vec3f <- unlist(lapply(res.oto.N.all.1000, quantile_upperlim))

class(vec2f)

table.1g.N_oto_1000 <- data.frame( median = vec1f,
                                   lowCI = vec2f,
                                   upperCI = vec3f)

table.1g.N_oto_1000
# View(table.1g.N_oto_1000)

# additional columns (text95%CI, CIW, CIW%,text95%CI)

attach(table.1g.N_oto_1000)

table.1g.N_oto_1000$CIW <- table.1g.N_oto_1000$upperCI -table.1g.N_oto_1000$lowCI      

attach(table.1g.N_oto_1000)

table.1g.N_oto_1000$CIW_perc <- (table.1g.N_oto_1000$CIW / table.1g.N_oto_1000$median) * 100        

attach(table.1g.N_oto_1000)


# round values (3 digits)                                                       
table.1g.N_oto_1000 <- round (table.1g.N_oto_1000[,1:5], 3)                                                    

#  create a new column `text95CI` ----------------                                                       

# columns to paste together for text
cols2 <- c("lowCI", "upperCI"  )

# create a new column `text95CI` with the three columns collapsed together (character string)
table.1g.N_oto_1000$text95CI <- apply( table.1g.N_oto_1000[ , cols2 ] , 1 , paste , collapse = "_to_" )


table.1g.N_oto_1000

# export (write to hard drive)   ---------------------

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

write.table(table.1g.N_oto_1000, file = "table.1g.N_oto_1000.csv")

table.1g.N_oto_1000





# plotting with fishboot (otoliths)  -----------------------


res$bootRaw <- res.oto.N.all.1000_ini


# univariate density plot of bootstrapped# pars


# par(opar)

res$bootRaw




univariate_density(res)


par(opar)



# Linf / K scatterhist GA

LinfK_scatterhist(res, Linf.breaks=50, K.breaks=50, 
                  phi.contour = TRUE, phi.contour.lty = 1, phi.contour.col = 8,
                  pt.col = rgb(0,0,0,0.1), pt.pch = ".",
                  xlim = c(0,400) )
#  
#  
# LinfK_scatterhist(res, Linf.breaks=50, K.breaks=50, 
#                   phi.contour = TRUE, phi.contour.lty = 1, phi.contour.col = 8,
#                   pt.col = rgb(0,0,0,0.1), pt.pch = "."
#                   , xlim = c(0,200))
# 


# par(opar)

res$bootRaw

# VBGF by time growth curve plot

par(opar)

res$bootRaw$t_anchor <-  res.oto.S.all.1000_ini$t0


CIinfo <- vbgfCI_time(
  res = res, 
  agemax = 12, CI = 95,
  perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
  ci.col = 1, ci.lty = 2, ci.lwd = 1, 
  maxd.col = 1, maxd.lty = 1, maxd.lwd = 2
)




######### --------
######### --------
######### --------


# II. Chapter II: Comparisons and tests ------------------------------
# II. Comparisons ----------------


# table.1a.MA3_Linfmin50.res_median_CI
# table.1b.MA5_Linfmin50.res_median_CI
# table.1c.MA5_Linfmin30.res_median_CI
# table.1d.MA3_Linfmin30.res_median_CI
# table.1e.Tag_recap1000
# table.1f.S_oto_1000
# table.1g.N_oto_1000

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# read in data (posteriors)   ------------------

# New results (seed = Sys.time!)

#MA3_Linfmin50
res_MA3_Linf50to200 <- read.csv ("expA_MA3_Linfmin50_Sys.Time.res_1000.csv", sep = ";")
# old.res_MA3_Linf50to200 <- read.csv ("res_MA3_Linf50to200.csv", sep = "")
head(res_MA3_Linf50to200)

#MA5_Linfmin50
res_MA5_Linf50to200 <- read.csv ("expB_MA5_Linfmin50_Sys.Time.res_1000.csv", sep = ";")
#old.MA5_res_EXP_S_slow_exp1_ALLSizes <- read.csv ("MA5_res_EXP_S_slow_exp1_ALLSizes.csv", sep = "")
head(res_MA5_Linf50to200)

#MA3_Linfmin30
# expD_MA3_Linfmin30_Sys.Time.res_1000
res_MA3_Linf30_200_n1000  <- read.csv ("expD_MA3_Linfmin30_Sys.Time.res_1000.csv", sep = ";")
#old.res_MA3_Linf30_200_n1000  <- read.csv ("res_MA3_Linf30_200_n1000.csv", sep = "")
head(res_MA3_Linf30_200_n1000)

#MA5_Linfmin30
# expC_MA5_Linfmin30_Sys.Time.res_1000.csv
res_MA5_EXP3_S_slowLinf30_200_n1000 <- read.csv ("expC_MA5_Linfmin30_Sys.Time.res_1000.csv", sep = ";")

#old.res_MA5_EXP3_S_slowLinf30_200_n1000 <- read.csv ("res_MA5_EXP3_S_slowLinf30_200_n1000.csv", sep = "")
head(res_MA5_EXP3_S_slowLinf30_200_n1000)


# # MA3_Linf50to200
# res_MA3_Linf50to200 <- read.csv ("res_MA3_Linf50to200.csv", sep = "")
# head(res_MA3_Linf50to200)
# 
# #MA5_Linf50to200
# MA5_res_EXP_S_slow_exp1_ALLSizes <- read.csv ("MA5_res_EXP_S_slow_exp1_ALLSizes.csv", sep = "")
# head(MA5_res_EXP_S_slow_exp1_ALLSizes)
# res_MA5_Linf50to200 <- MA5_res_EXP_S_slow_exp1_ALLSizes
# head(res_MA5_Linf50to200)
# 
# #MA3_Linf30_200
# res_MA3_Linf30_200_n1000  <- read.csv ("res_MA3_Linf30_200_n1000.csv", sep = "")
# head(res_MA3_Linf30_200_n1000)
# 
# #MA5_Linfmin30
# res_MA5_EXP3_S_slowLinf30_200_n1000 <- read.csv ("res_MA5_EXP3_S_slowLinf30_200_n1000.csv", sep = "")
# head(res_MA5_EXP3_S_slowLinf30_200_n1000)


res_E_1000StoutputB1_9_from_Margit <- read.csv("res_E_Tag_recap1000.csv", sep="")
res_E_Tag_recap1000 <- res_E_1000StoutputB1_9_from_Margit
head(res_E_Tag_recap1000)


res.oto.S.all.1000 <- read.csv("oto-S-all-1000results.csv", sep=",")
head(res.oto.S.all.1000)#1000, OK

res_MA3_Linfmin50 <- res_MA3_Linf50to200
length(res_MA3_Linfmin50$Linf)#1000, OK

res_MA5_Linfmin50  <- MA5_res_EXP_S_slow_exp1_ALLSizes # Linfmin50cm
length(res_MA5_Linfmin50$Linf)#1000, OK

res_MA3_Linfmin30   <-  res_MA3_Linf30_200_n1000
length(res_MA3_Linfmin30$Linf)#1000, OK

res_MA5_Linfmin30   <-  res_MA5_EXP3_S_slowLinf30_200_n1000
length(res_MA5_Linfmin30$Linf)# 1000, OK
head(res_MA5_Linfmin30)

res_E_Tag_recap1000
length(res_E_Tag_recap1000$Linf)#
res_E_Tag_recap1000 <- res_E_Tag_recap1000[1:1000,]
head(res_E_Tag_recap1000)
length(res_E_Tag_recap1000$Linf)#1000, OK


res.oto.S.all.1000$Phi_L <-  log10(res.oto.S.all.1000$K) + (2 * log10(res.oto.S.all.1000$Linf))
res.oto.S.all.1000
head(res.oto.S.all.1000)
length(res_E_Tag_recap1000$Linf)#1000, OK

# Make Tables organized by variables (Linf, K, Phi') ------------------------

Linf_table <- data.frame( MA3_Linfmin50 = res_MA3_Linfmin50$Linf, 
                          MA5_Linfmin50 = res_MA5_Linfmin50$Linf,
                          MA3_Linfmin30 = res_MA3_Linfmin30$Linf,
                          MA5_Linfmin30 = res_MA5_Linfmin30$Linf,
                          Tag_recap_S =   res_E_Tag_recap1000$Linf,
                          Otoliths_S  =   res.oto.S.all.1000$Linf  )

head(Linf_table)


K_table  <- data.frame( MA3_Linfmin50 = res_MA3_Linfmin50$K, 
                            MA5_Linfmin50 = res_MA5_Linfmin50$K,
                            MA3_Linfmin30 = res_MA3_Linfmin30$K,
                            MA5_Linfmin30 = res_MA5_Linfmin30$K,
                            Tag_recap_S =   res_E_Tag_recap1000$K,
                            Otoliths_S  =   res.oto.S.all.1000$K  )

head(K_table)



Phi_table  <- data.frame( MA3_Linfmin50 = res_MA3_Linfmin50$phiL, 
                        MA5_Linfmin50 = res_MA5_Linfmin50$phiL,
                        MA3_Linfmin30 = res_MA3_Linfmin30$phiL,
                        MA5_Linfmin30 = res_MA5_Linfmin30$phiL,
                        Tag_recap_S =   res_E_Tag_recap1000$Phi_L,
                        Otoliths_S  =   res.oto.S.all.1000$Phi_L  )

head(Phi_table)

# II.1 compare medians (are medians different between methods?)   ----------------

# Medians - Test for difference in Medians, based on the bootstrap posteriors---------
# use the "bootstrapped.median.test"
# pairwise.bootstrapped.median.test --------------------

# #  The Bootstrapped two-sample test, based on posteriors 
# Ralf Schwamborn
# March 20, 2022
# Citation: 

# Schwamborn, R., 2022. The bootstrapped two-sample test. 
#〈https://rpubs.com/rschwamborn/880212 (accessed 22 March, 2022).

#Source: https://rpubs.com/rschwamborn/880212

#see also:
# Schwamborn, R., 2019. The interquantile range test. 
#〈http://rpubs.com/rschwamborn/515419〉 (accessed 25 July 2019).


# define the function "boot.p.two.sample.medians.diff.posteriors()" ------------------


boot.p.two.sample.medians.diff.posteriors <- function (boot.post.B, boot.post.A) {
  
  # difference between posteriors
  
  boot.statisticsA <- boot.post.A # posterior for sample A
  boot.statisticsB <- boot.post.B # posterior for sample B
  
  if ( median(boot.post.B) < median(boot.post.A) ) {
  # data set "B", has larger values
  
  boot.statisticsDiff <- boot.statisticsA - boot.statisticsB
  N.below.zero <- length(boot.statisticsDiff [boot.statisticsDiff< 0]) # 7 are below 0 
  
  N.total <- length(boot.statisticsDiff ) # all 
  
  (one.sided.p_value <- N.below.zero / N.total) # one-sided "p" value
  
  (two.sided.p_value <- 2* (N.below.zero / N.total)) # two-sided "p" value
  
  ans <-   data.frame( text = "two.sided.p_value", p.value = two.sided.p_value)  # output
  }
  
  
  
  else 
    {boot.statisticsDiff <- boot.statisticsB - boot.statisticsA
    
  # calculate "p" value  from bootstrap diff.  
  
  N.below.zero <- length(boot.statisticsDiff [boot.statisticsDiff< 0]) # 7 are below 0 
  
  N.total <- length(boot.statisticsDiff ) # all 
  
  (one.sided.p_value <- N.below.zero / N.total) # one-sided "p" value
  
  (two.sided.p_value <- 2* (N.below.zero / N.total)) # two-sided "p" value
  
  ans <-   data.frame(text = "two.sided.p_value", p.value = two.sided.p_value)  # output
    }
  
 paste ( as.list(ans) )   # output
  
}





boot.p.two.sample.medians.diff.posteriors_B <- function (boot.post.B, boot.post.A) {
  
  # difference between posteriors
  
  boot.statisticsA <- boot.post.A # posterior for sample A
  boot.statisticsB <- boot.post.B # posterior for sample B
  
  if ( median(boot.post.B) < median(boot.post.A) ) {
    # data set "B", has larger values
    
    boot.statisticsDiff <- boot.statisticsA - boot.statisticsB
    N.below.zero <- length(boot.statisticsDiff [boot.statisticsDiff<= 0]) # 7 are below 0 
    
    N.total <- length(boot.statisticsDiff ) # all 
    
    (one.sided.p_value <- N.below.zero / N.total) # one-sided "p" value
    
    (two.sided.p_value <- 2* (N.below.zero / N.total)) # two-sided "p" value
    
    ans <-   data.frame( text = "two.sided.p_value", p.value = two.sided.p_value)  # output
  }
  
  
  
  else 
  {boot.statisticsDiff <- boot.statisticsB - boot.statisticsA
  
  # calculate "p" value  from bootstrap diff.  
  
  N.below.zero <- length(boot.statisticsDiff [boot.statisticsDiff <= 0]) # n below 0 or equal 0 
  
  N.total <- length(boot.statisticsDiff ) # all 
  
  (one.sided.p_value <- N.below.zero / N.total) # one-sided "p" value
  
  (two.sided.p_value <- 2* (N.below.zero / N.total)) # two-sided "p" value
  
  ans <-   data.frame( text = "two.sided.p_value", p.value = two.sided.p_value)  # output
  }
  
  paste ( as.list(ans) )   # output
  
}



# test the function "boot.p.two.sample.medians.diff.posteriors()" ----------

boot.p.two.sample.medians.diff.posteriors (Linf_table$MA3_Linfmin30 ,Linf_table$MA3_Linfmin50 )

boot.p.two.sample.medians.diff.posteriors (Linf_table$MA3_Linfmin30 ,Linf_table$MA5_Linfmin50 )

boot.p.two.sample.medians.diff.posteriors (Linf_table$Otoliths_S ,Linf_table$MA5_Linfmin50 )

boot.p.two.sample.medians.diff.posteriors (Linf_table$Otoliths_S ,Linf_table$Tag_recap_S )

### LOOP ---------------
# Pairwise Bootstrapped two-sample test ---------------------
# make a loop for pairwise Bootstrapped two-sample test

# create an empty results matrix (filled with "NAs")

res_p_mat_Linf <- matrix(nrow = length(Linf_table) , ncol = length(Linf_table))
colnames(res_p_mat_Linf) <- names(Linf_table)
rownames(res_p_mat_Linf) <- names(Linf_table)

res_p_mat_Linf[1,2]

# fill in one row 
w = 2

p_res <- 1:length(Linf_table)

for (w in 1: length(Linf_table) ) {
p_res[w] <-  boot.p.two.sample.medians.diff.posteriors_B (K_table[,2], K_table[,w]            )[2] 
}
p_res

res_p_mat_Linf[1,] <- p_res
res_p_mat_Linf

# all rows (loop in loop) ----------------

# pairwise tests for Linf (difference in medians)    ------------
# all rows (loop in loop)----------------


for ( row in 1: length(Linf_table) ) {
  
  
  for (w in 1: length(Linf_table) ) {
    p_res[w] <-  boot.p.two.sample.medians.diff.posteriors_B (Linf_table[,row], Linf_table[,w]            )[2] 
  }
  p_res
  res_p_mat_Linf[row,] <- p_res
  
  
}

res_p_mat_Linf # "p" values are shown, for differences in medians of the bootstrap posteriorsLinf
      
# 
#                  MA3_Linfmin50 MA5_Linfmin50 MA3_Linfmin30 MA5_Linfmin30 Tag_recap_S Otoliths_S
# MA3_Linfmin50 "0"           "0.72"        "0.884"       "0.516"       "0.302"     "0.734"   
# MA5_Linfmin50 "0.72"        "0"           "0.632"       "0.342"       "0.278"     "0.62"    
# MA3_Linfmin30 "0.884"       "0.632"       "0"           "0.576"       "0.32"      "0.784"   
# MA5_Linfmin30 "0.516"       "0.342"       "0.576"       "0"           "0.994"     "0.636"   
# Tag_recap_S   "0.302"       "0.278"       "0.32"        "0.994"       "0"         "0.366"   
# Otoliths_S    "0.734"       "0.62"        "0.784"       "0.636"       "0.366"     "0"    


# pairwise tests for K (medians)    -------------
# all rows (loop in loop)----------------

res_p_mat_K <- matrix(nrow = length(K_table) , ncol = length(K_table))
colnames(res_p_mat_K) <- names(K_table)
rownames(res_p_mat_K) <- names(K_table)

p_res <- 1:length(K_table)

   

for ( row in 1: length(K_table) ) {
  
  
  for (w in 1: length(K_table) ) {
    p_res[w] <-  boot.p.two.sample.medians.diff.posteriors (K_table[,row], K_table[,w]            )[2] 
  }
  p_res
  res_p_mat_K[row,] <- p_res
  
  
}

res_p_mat_K # "p" values are shown, for differences in medians of the bootstrap posteriors  in K
# 
#                MA3_Linfmin50 MA5_Linfmin50 MA3_Linfmin30 MA5_Linfmin30 Tag_recap_S Otoliths_S
# MA3_Linfmin50 "0"           "0.808"       "0.868"       "0.368"       "0.644"     "0.498"   
# MA5_Linfmin50 "0.808"       "0"           "0.928"       "0.414"       "0.802"     "0.368"   
# MA3_Linfmin30 "0.868"       "0.928"       "0"           "0.412"       "0.684"     "0.436"   
# MA5_Linfmin30 "0.368"       "0.414"       "0.412"       "0"           "0.544"     "0.18"    
# Tag_recap_S   "0.644"       "0.802"       "0.684"       "0.544"       "0"         "0.38"    
# Otoliths_S    "0.498"       "0.368"       "0.436"       "0.18"        "0.38"      "0"  



# pairwise tests for Phi' (medians)    -------------
# all rows (loop in loop)----------------

res_p_mat_Phi <- matrix(nrow = length(Phi_table ) , ncol = length(Phi_table))
colnames(res_p_mat_Phi) <- names(Phi_table)
rownames(res_p_mat_Phi) <- names(Phi_table)

p_res <- 1:length(Phi_table)


for ( row in 1: length(Phi_table) ) {
  
  for (w in 1: length(Phi_table) ) {
    p_res[w] <-  boot.p.two.sample.medians.diff.posteriors (Phi_table[,row], Phi_table[,w]            )[2] 
  }
  p_res
  res_p_mat_Phi[row,] <- p_res
  
}

res_p_mat_Phi # "p" values are shown,  for differences in medians of the bootstrap posteriors in Phi'
# 
#                MA3_Linfmin50 MA5_Linfmin50 MA3_Linfmin30 MA5_Linfmin30 Tag_recap_S Otoliths_S
# MA3_Linfmin50 "0"           "0.35"        "0.926"       "0.27"        "0.266"     "0.052"   
# MA5_Linfmin50 "0.35"        "0"           "0.442"       "0.896"       "0.074"     "0.024"   
# MA3_Linfmin30 "0.926"       "0.442"       "0"           "0.304"       "0.262"     "0.044"   
# MA5_Linfmin30 "0.27"        "0.896"       "0.304"       "0"           "0.052"     "0.016"   
# Tag_recap_S   "0.266"       "0.074"       "0.262"       "0.052"       "0"         "0.56"    
# Otoliths_S    "0.052"       "0.024"       "0.044"       "0.016"       "0.56"      "0"   
# 

# Descriptive Text and conclusions (diff in medians) -----------

# Conclusions: no differences in medians for Linf and K (all medians were similar between methods), 
            # Explanation (intepretation): high uncertainty in Linf and K (very large CI widths in all methods )

#             but there are significnt differences in median Phi'  
#             for Otoliths_S vs MA5_Linfmin50 (p = 0.024), Otoliths_S vs  MA3_Linfmin30 ( p = 0.044),
#                and Otoliths_S vs  MA5_Linfmin30 (p = 0.016). 
#             Explanation (intepretation):
#             Otoliths analysis yielded  much lower median Phi' (median Phi'  = 2.617)  than LFA (median Phi'  values of 2.894,	2.903,	3.08, and	3.143, depending on MA and Linfmin)
#                             Tag-recapture (median Phi': 2.717) produced very similar (not significantly different) median Phi values as  Otoliths LAA (median Phi'  = 2.617).

# Explanation (interpretation): Our results indicate a very high uncertainty in Linf and K (very large CI widths in all methods ),
#                               However, confidence envelopes are well aligned along Phi' isolines (Kimura plots) , 
#                               thus much lower CI widths fpr Phi' 
#                                and much lower uncertainty  regarding  Phi' (also, Phi' is based on log-transformed values)



# Text for  Results chapter:

# There were no significant differences in medians (p > 0.05, Bootstrapped two-sample test), between methods (LFA, LAA and T&R),  for Linf and K (i.e., all median Linf and K posteriors were similar between methods).
# Also, median Phi' posteriors were not significantly different 
# between otoliths (LAA), Tag-recapture, and LFA with MA = 3 and Linfmin = 50. 
#        However,  we detected  significant differences in median Phi' posteriors 
#             for Otoliths vs LFA with MA = 5 and Linfmin = 50 (p = 0.024), Otoliths vs  MA = 3 and Linfmin = 30 ( p = 0.044),
#                and Otolith vs  MA = 5 and Linfmin = 30 (p = 0.016), where Otoliths analysis (median Phi'  = 2.617)  showed significantly lower median Phi' than these three LFA analyses . 

# Text for  Discussion chapter:

# The results of pairwise Bootstrapped two-sample tests indicate a very high uncertainty in Linf and K (very large CI widths in all methods).
# However, confidence envelopes were well aligned along Phi' isolines (Kimura plots, Fig XX) , 
# thus producing much lower CI widths for Phi', than for K and Linf, 
# and much lower uncertainty  regarding  Phi' (also, Phi' is based on log-transformed values). 
# This explains why we could detect significant differences in median Phi' posteriors 
# for Otoliths vs LFA for three different LFA analysis settings, with significant lower growth performance obtained from otoliths than from LFA. 
# Also, this highlights the need to conduct LFA always testing different search space limits and ELEFAN settings (REFs).  




# II.1.b compare CI widths (CI widths are different?)  MA3 vs MA5 -------------

# The interquantile range test, V01
# Ralf Schwamborn
# July 25, 2019
# The interquantile range test, Version 0.1.
# (Copyright: R. Schwamborn, 2018)
# This is a “R” script that describes the interquantile range test, 
# as implemented the function interquant_r.test, is a significance test 
# to verify significant differences in 95% interquantile range.
# More specifically, the range from the 0.025 to 0.975 quantile is 
# analysed, that encompassed 95% of the data. This test is intended to 
# compare the precision of statistical procedures (ie., to test for 
# diffences in 95% confidence intervals, based on bootstrap posteriors).

# https://rpubs.com/rschwamborn/515419

#install.packages("WRS2")
library(WRS2)
##

interquant_r.test <- function(dat.A, dat.B, n.boot) {
  
  
  dat.A.stand <-   dat.A - quantile(dat.A, 0.025 )                             
  
  dat.B.stand <-   dat.B - quantile(dat.B, 0.025 )                             
  
  
  dat.frameAB <- data.frame(dat.A.stand, dat.B.stand)
  dat.frameAB.stack <-  stack(dat.frameAB)
  
  
  res <- Qanova(values ~ ind, dat.frameAB.stack, q =  0.975, nboot = n.boot)
  
  res.p.value <-   res$p.value[1,2]
  
  
  if (res.p.value == 0)
    
    paste("p-value < 0.0001")
  
  else
    
    paste ("p-value = ", res.p.value )
  
}


# Compare Linf (CI width MA 3 vs MA5),  ---------------------
#test for differences in precision

interquant_r.test(res_E_Tag_recap1000$Linf, res.oto.S.all.1000$Linf , 1000)

res_p <- interquant_r.test(res_E_Tag_recap1000$Linf, res.oto.S.all.1000$Linf , 1000)



# par(opar)

two.posteriorsLinf <- data.frame(T_R = res_E_Tag_recap1000$Linf, LAA = res.oto.S.all.1000$Linf)

mydata <- stack(two.posteriorsLinf)
summary(mydata)
boxplot(mydata$values ~ mydata$ind, data = mydata,
        xlab= "Method", ylab = "Linf") 

library(vioplot)
vioplot(mydata$values ~ mydata$ind, data = mydata,
        xlab= "Method", ylab = "Linf") 


quantile(res.oto.S.all.1000$Linf, c(0.025, 0.5, 0.975))
quantile(res_E_Tag_recap1000$Linf, c(0.025, 0.5, 0.975))


# vioplot with quantiles (median and 95% CI)
vioplot(mydata$values ~ mydata$ind, data = mydata,
        xlab= "MA Value", ylab = "Linf") 

abline(col = "blue", h =  quantile(res.oto.S.all.1000$Linf, 
                                   c(0.025, 0.5, 0.975)))
abline(col = "red", h =  quantile(res_E_Tag_recap1000$Linf, 
                                  c(0.025, 0.5, 0.975)))

### LOOP ---------------
# Pairwise interquantile range test ---------------------
# make a loop for pairwise interquantile range tests

# create an empty results matrix (filled with "NAs")

res_p_mat_Linf <- matrix(nrow = length(Linf_table) , ncol = length(Linf_table))
colnames(res_p_mat_Linf) <- names(Linf_table)
rownames(res_p_mat_Linf) <- names(Linf_table)

res_p_mat_Linf[1,2]

# fill in one row 
w = 2

t1 <- Sys.time()

p_res <- 1:length(Linf_table)

for (w in 1: length(Linf_table) ) {
  p_res[w] <-  interquant_r.test (Linf_table[,1], Linf_table[,w]           ,1000 )
}
p_res


res_p_mat_Linf[1,] <- p_res
res_p_mat_Linf
t2 <- Sys.time()

(duration <- t2-t1) # duration: 38.4 secs per row


# all rows (loop in loop) ----------------

# pairwise tests for Linf (difference in medians)    ------------
# all rows (loop in loop)----------------


t1 <- Sys.time()


for ( row in 1: length(Linf_table) ) {
  
  
  for (w in 1: length(Linf_table) ) {
    p_res[w] <-  interquant_r.test (Linf_table[,row], Linf_table[,w], 1000 )
  }
  p_res
  res_p_mat_Linf[row,] <- p_res
  
  
}

t2 <- Sys.time()
(duration <- t2-t1) # duration: Time difference of 6.04 mins


res_p_mat_Linf # "p" values are shown, for Linf

res_p_mat_Linf_interquantile <- res_p_mat_Linf

write.table(  as.data.frame(res_p_mat_Linf_interquantile), file = "res_p_mat_Linf_interquantile.csv" )

# Bootstrapped pairwise interquantile range test ("p" for difference in 95% Confidence interval widths,
# i.e., differences in precision in Linf)
# 
# res_p_mat_Linf # "p" values are shown, for Linf

# MA3_Linfmin50      MA5_Linfmin50      MA3_Linfmin30      MA5_Linfmin30      Tag_recap_S       
# MA3_Linfmin50 "p-value =  0.999" "p-value < 0.0001" "p-value =  0.417" "p-value =  0.001" "p-value =  0.001"
# MA5_Linfmin50 "p-value < 0.0001" "p-value =  0.993" "p-value < 0.0001" "p-value =  0.041" "p-value =  0.017"
# MA3_Linfmin30 "p-value =  0.452" "p-value < 0.0001" "p-value =  0.965" "p-value < 0.0001" "p-value =  0.001"
# MA5_Linfmin30 "p-value =  0.001" "p-value =  0.037" "p-value < 0.0001" "p-value =  0.957" "p-value =  0.005"
# Tag_recap_S   "p-value =  0.001" "p-value =  0.013" "p-value =  0.001" "p-value =  0.005" "p-value =  0.992"
# Otoliths_S    "p-value < 0.0001" "p-value < 0.0001" "p-value < 0.0001" "p-value < 0.0001" "p-value =  0.056"
# Otoliths_S        
# MA3_Linfmin50 "p-value < 0.0001"
# MA5_Linfmin50 "p-value < 0.0001"
# MA3_Linfmin30 "p-value < 0.0001"
# MA5_Linfmin30 "p-value < 0.0001"
# Tag_recap_S   "p-value =  0.079"
# Otoliths_S    "p-value =  0.997"
# > 


# Results of Bootstrapped two-sample test ("p" for difference in medians of posteriors of Linf estimates), 
#  i.e., differences in accuracy)

#                  MA3_Linfmin50 MA5_Linfmin50 MA3_Linfmin30 MA5_Linfmin30 Tag_recap_S Otoliths_S
# MA3_Linfmin50 "0"           "0.72"        "0.884"       "0.516"       "0.302"     "0.734"   
# MA5_Linfmin50 "0.72"        "0"           "0.632"       "0.342"       "0.278"     "0.62"    
# MA3_Linfmin30 "0.884"       "0.632"       "0"           "0.576"       "0.32"      "0.784"   
# MA5_Linfmin30 "0.516"       "0.342"       "0.576"       "0"           "0.994"     "0.636"   
# Tag_recap_S   "0.302"       "0.278"       "0.32"        "0.994"       "0"         "0.366"   
# Otoliths_S    "0.734"       "0.62"        "0.784"       "0.636"       "0.366"     "0"    


# pairwise tests for K (medians) -------------
# all rows (loop in loop) ----------------

t1 <- Sys.time()

res_p_mat_K <- matrix(nrow = length(K_table) , ncol = length(K_table))
colnames(res_p_mat_K) <- names(K_table)
rownames(res_p_mat_K) <- names(K_table)

p_res <- 1:length(K_table)



for ( row in 1: length(K_table) ) {
  
  
  for (w in 1: length(K_table) ) {
    p_res[w] <-  interquant_r.test (K_table[,row], K_table[,w]  ,1000   )
  }
  p_res
  res_p_mat_K[row,] <- p_res
  
  
}


t2 <- Sys.time()
(duration <- t2-t1) # duration: 3.3 minutes



res_p_mat_K # "p" values are shown, for K

res_p_mat_K_interquantile <- res_p_mat_K

write.table(  as.data.frame(res_p_mat_K_interquantile), file = "res_p_mat_K_interquantile.csv" )

# Bootstrapped pairwise interquantile range test ("p" for difference in 95% Confidence interval widths,
# i.e., differences in precision in K)
# 
# 
# 
# res_p_mat_K # "p" values are shown, for K

#                MA3_Linfmin50      MA5_Linfmin50      MA3_Linfmin30      MA5_Linfmin30      Tag_recap_S       
# MA3_Linfmin50 "p-value =  0.972" "p-value =  0.002" "p-value =  0.482" "p-value < 0.0001" "p-value < 0.0001"
# MA5_Linfmin50 "p-value =  0.002" "p-value =  0.984" "p-value < 0.0001" "p-value < 0.0001" "p-value =  0.84" 
# MA3_Linfmin30 "p-value =  0.5"   "p-value < 0.0001" "p-value =  0.988" "p-value < 0.0001" "p-value < 0.0001"
# MA5_Linfmin30 "p-value < 0.0001" "p-value < 0.0001" "p-value < 0.0001" "p-value =  0.974" "p-value < 0.0001"
# Tag_recap_S   "p-value < 0.0001" "p-value =  0.88"  "p-value < 0.0001" "p-value < 0.0001" "p-value =  0.976"
# Otoliths_S    "p-value < 0.0001" "p-value < 0.0001" "p-value < 0.0001" "p-value < 0.0001" "p-value < 0.0001"
# Otoliths_S        
# MA3_Linfmin50 "p-value < 0.0001"
# MA5_Linfmin50 "p-value < 0.0001"
# MA3_Linfmin30 "p-value < 0.0001"
# MA5_Linfmin30 "p-value < 0.0001"
# Tag_recap_S   "p-value < 0.0001"
# Otoliths_S    "p-value =  0.979"
 



# Results of Bootstrapped two-sample test ("p" for difference in medians in K), 
#  i.e., differences in accuracy in K)
# 
#                MA3_Linfmin50 MA5_Linfmin50 MA3_Linfmin30 MA5_Linfmin30 Tag_recap_S Otoliths_S
# MA3_Linfmin50 "0"           "0.808"       "0.868"       "0.368"       "0.644"     "0.498"   
# MA5_Linfmin50 "0.808"       "0"           "0.928"       "0.414"       "0.802"     "0.368"   
# MA3_Linfmin30 "0.868"       "0.928"       "0"           "0.412"       "0.684"     "0.436"   
# MA5_Linfmin30 "0.368"       "0.414"       "0.412"       "0"           "0.544"     "0.18"    
# Tag_recap_S   "0.644"       "0.802"       "0.684"       "0.544"       "0"         "0.38"    
# Otoliths_S    "0.498"       "0.368"       "0.436"       "0.18"        "0.38"      "0"  



# pairwise tests for Phi' (medians)    -------------
# all rows (loop in loop)----------------

res_p_mat_Phi <- matrix(nrow = length(Phi_table ) , ncol = length(Phi_table))
colnames(res_p_mat_Phi) <- names(Phi_table)
rownames(res_p_mat_Phi) <- names(Phi_table)

p_res <- 1:length(Phi_table)

t1 <- Sys.time()

for ( row in 1: length(Phi_table) ) {
  
  for (w in 1: length(Phi_table) ) {
    p_res[w] <-  interquant_r.test (Phi_table[,row], Phi_table[,w]   ,1000         )
  }
  p_res
  res_p_mat_Phi[row,] <- p_res
  
}


t2 <- Sys.time()
(duration <- t2-t1) # duration: 3.9 minutes


res_p_mat_Phi_interquantile <- res_p_mat_Phi

write.table(  as.data.frame(res_p_mat_Phi_interquantile), file = "res_p_mat_Phi_interquantile.csv" )



res_p_mat_Phi # "p" values are shown, for Phi'

# Bootstrapped pairwise interquantile range test ("p" for difference in 95% Confidence interval widths,
# i.e., differences in precision in Phi')
# 
# 
# res_p_mat_Phi # "p" values are shown, for Phi'
# 
#   > res_p_mat_Phi # "p" values are shown, for Phi'
# MA3_Linfmin50      MA5_Linfmin50      MA3_Linfmin30      MA5_Linfmin30                  
# MA3_Linfmin50 "p-value =  0.986" "p-value =  0.809" "p-value =  0.215" "p-value < 0.0001"             
# MA5_Linfmin50 "p-value =  0.804" "p-value =  0.989" "p-value =  0.227" "p-value < 0.0001"             
# MA3_Linfmin30 "p-value =  0.198" "p-value =  0.274" "p-value =  0.994" "p-value < 0.0001"             
# MA5_Linfmin30 "p-value < 0.0001" "p-value < 0.0001" "p-value < 0.0001" "p-value =  0.957"             
# Tag_recap_S   "p-value < 0.0001" "p-value < 0.0001" "p-value < 0.0001" "p-value < 0.0001"             
# Otoliths_S    "p-value < 0.0001" "p-value < 0.0001" "p-value =  0.001" "p-value =  0.0600000000000001"
# Tag_recap_S        Otoliths_S                     
# MA3_Linfmin50 "p-value < 0.0001" "p-value =  0.001"             
# MA5_Linfmin50 "p-value < 0.0001" "p-value < 0.0001"             
# MA3_Linfmin30 "p-value < 0.0001" "p-value =  0.005"             
# MA5_Linfmin30 "p-value < 0.0001" "p-value =  0.0620000000000001"
# Tag_recap_S   "p-value =  1"     "p-value < 0.0001"             
# Otoliths_S    "p-value < 0.0001" "p-value =  0.958"             
# > 
# 
# 



# Results of Bootstrapped two-sample test ("p" for difference in medians in Phi'), 
#  i.e., differences in accuracy in Phi')
# 

#                MA3_Linfmin50 MA5_Linfmin50 MA3_Linfmin30 MA5_Linfmin30 Tag_recap_S Otoliths_S
# MA3_Linfmin50 "0"           "0.35"        "0.926"       "0.27"        "0.266"     "0.052"   
# MA5_Linfmin50 "0.35"        "0"           "0.442"       "0.896"       "0.074"     "0.024"   
# MA3_Linfmin30 "0.926"       "0.442"       "0"           "0.304"       "0.262"     "0.044"   
# MA5_Linfmin30 "0.27"        "0.896"       "0.304"       "0"           "0.052"     "0.016"   
# Tag_recap_S   "0.266"       "0.074"       "0.262"       "0.052"       "0"         "0.56"    
# Otoliths_S    "0.052"       "0.024"       "0.044"       "0.016"       "0.56"      "0"   
# 




# Descriptive Text and conclusions (differences in CIW, i.e. in precision) -----------




#########
########
######
####
# III. Chapter III: Plots, customized boxplots ------------------------
# III. 1. Plots  of 95% Confidence intervals (i.e., 95% quantiles of the posteriors)

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# read in data (posteriors)   ------------------
# New results (seed = Sys.time!)

#MA3_Linfmin50
res_MA3_Linf50to200 <- read.csv ("expA_MA3_Linfmin50_Sys.Time.res_1000.csv", sep = ";")

# old.res_MA3_Linf50to200 <- read.csv ("res_MA3_Linf50to200.csv", sep = "")
head(res_MA3_Linf50to200)

#MA5_Linfmin50
res_MA5_Linf50to200 <- read.csv ("expB_MA5_Linfmin50_Sys.Time.res_1000.csv", sep = ";")

#old.MA5_res_EXP_S_slow_exp1_ALLSizes <- read.csv ("MA5_res_EXP_S_slow_exp1_ALLSizes.csv", sep = "")
head(res_MA5_Linf50to200)

#MA3_Linfmin30
# expD_MA3_Linfmin30_Sys.Time.res_1000
res_MA3_Linf30_200_n1000  <- read.csv ("expD_MA3_Linfmin30_Sys.Time.res_1000.csv", sep = ";")

#old.res_MA3_Linf30_200_n1000  <- read.csv ("res_MA3_Linf30_200_n1000.csv", sep = "")
head(res_MA3_Linf30_200_n1000)

#MA5_Linfmin30
# expC_MA5_Linfmin30_Sys.Time.res_1000.csv
res_MA5_EXP3_S_slowLinf30_200_n1000 <- read.csv ("expC_MA5_Linfmin30_Sys.Time.res_1000.csv", sep = ";")

#old.res_MA5_EXP3_S_slowLinf30_200_n1000 <- read.csv ("res_MA5_EXP3_S_slowLinf30_200_n1000.csv", sep = "")
head(res_MA5_EXP3_S_slowLinf30_200_n1000)

res_E_1000StoutputB1_9_from_Margit <- read.csv("res_E_Tag_recap1000.csv", sep="")
res_E_Tag_recap1000 <- res_E_1000StoutputB1_9_from_Margit
head(res_E_Tag_recap1000)


res.oto.S.all.1000 <- read.csv("~/Papers/fishboot - Namibia - paper/LAA_SandN_Stock/oto-S-all-1000results.csv", sep=",")
head(res.oto.S.all.1000)#1000, OK


oto.N.all.1000results <- read.csv("~/Papers/fishboot - Namibia - paper/LAA_SandN_Stock/oto-N-all-1000results.csv")
res.oto.N.all.1000 <- oto.N.all.1000results
#   View(res.oto.S.all.1000)
head(res.oto.N.all.1000)#1000, OK



res_MA3_Linfmin50 <- res_MA3_Linf50to200
length(res_MA3_Linfmin50$Linf)#1000, OK

res_MA5_Linfmin50  <- res_MA5_Linf50to200 # Linfmin50cm
length(res_MA5_Linfmin50$Linf)#1000, OK

res_MA3_Linfmin30   <-  res_MA3_Linf30_200_n1000
length(res_MA3_Linfmin30$Linf)#1000, OK

res_MA5_Linfmin30   <-  res_MA5_EXP3_S_slowLinf30_200_n1000
length(res_MA5_Linfmin30$Linf)# 1000, OK
head(res_MA5_Linfmin30)

res_E_Tag_recap1000
length(res_E_Tag_recap1000$Linf)#
res_E_Tag_recap1000 <- res_E_Tag_recap1000[1:1000, ]
head(res_E_Tag_recap1000)
length(res_E_Tag_recap1000$Linf)#1000, OK


res.oto.S.all.1000$Phi_L <-  log10(res.oto.S.all.1000$K) + (2 * log10(res.oto.S.all.1000$Linf))

#res.oto.S.all.1000
head(res.oto.S.all.1000)
length(res.oto.S.all.1000$Linf)#1000, OK


res.oto.N.all.1000$Phi_L <-  log10(res.oto.N.all.1000$K) + (2 * log10(res.oto.N.all.1000$Linf))
#res.oto.N.all.1000
head(res.oto.N.all.1000)
length(res.oto.N.all.1000$Linf)#1000, OK



# Make Tables organized by variables (Linf, K, Phi') ------------------------

Linf_table <- data.frame( MA3_Linfmin50 = res_MA3_Linfmin50$Linf, 
                          MA5_Linfmin50 = res_MA5_Linfmin50$Linf,
                          MA3_Linfmin30 = res_MA3_Linfmin30$Linf,
                          MA5_Linfmin30 = res_MA5_Linfmin30$Linf,
                          Tag_recap_S =   res_E_Tag_recap1000$Linf,
                          Otoliths_S  =   res.oto.S.all.1000$Linf,
                          Otoliths_N  =   res.oto.N.all.1000$Linf)

head(Linf_table)


K_table  <- data.frame( MA3_Linfmin50 = res_MA3_Linfmin50$K, 
                        MA5_Linfmin50 = res_MA5_Linfmin50$K,
                        MA3_Linfmin30 = res_MA3_Linfmin30$K,
                        MA5_Linfmin30 = res_MA5_Linfmin30$K,
                        Tag_recap_S =   res_E_Tag_recap1000$K,
                        Otoliths_S  =   res.oto.S.all.1000$K,
                        Otoliths_N  =   res.oto.N.all.1000$K  )

head(K_table)



Phi_table  <- data.frame( MA3_Linfmin50 = res_MA3_Linfmin50$phiL, 
                          MA5_Linfmin50 = res_MA5_Linfmin50$phiL,
                          MA3_Linfmin30 = res_MA3_Linfmin30$phiL,
                          MA5_Linfmin30 = res_MA5_Linfmin30$phiL,
                          Tag_recap_S =   res_E_Tag_recap1000$Phi_L,
                          Otoliths_S  =   res.oto.S.all.1000$Phi_L ,
                          Otoliths_N  =   res.oto.N.all.1000$Phi_L )

head(Phi_table)



# code for plotting 95% quantiles:
# https://forum.posit.co/t/quantile-box-plot-which-is-not-an-outlier-box-plot/24460/5


# plots for Linf -------------------

# 95% quantile boxplot:
 

NewBoxPlot.ylim.Linf <- function(x,y, ylim = ylim){
  y <- as.character(y)
  nbox <- length(unique(y))
  f <- unique(y)
  lims <- c(0 - .5, nbox -.5)
  cents <- 0:(nbox-1)
  labs.y <- seq( 0 , max (x) , by = 50)# y labels

  #labs.y <- seq( 40 , quantile (x, 0.9999) , by = 20)# y labels
    plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '',xlim = lims, ylim = ylim)
  for (i in 1:nbox){
    xi <- x[y==f[i]]
    cent = cents[i]
    b <- quantile(xi, c(.025,.5,.975)) # quantiles
    mima <- range(xi)                        # max and min
    #rect(cent - .2, b[3], cent + .2,  border = 'firebrick')
   # points(cent + c(-.2,.2), c(b[4], b[4]), type = 'l', col = 'firebrick', lwd = 3)
  #  arrows(cent, b[5], cent, mima[2], col = 'firebrick', angle = 90)
    #arrows(cent, b[3], cent, mima[1], col = 'firebrick', angle = 90)
    points(rep(cent, 1), b[c(2)], pch = '-', cex = 1.5, col = 'navyblue')
        points(rep(cent, 2), b[c(1,3)], pch = '-', cex = 1.5, col = 'firebrick')
  }
  axis(side = 1, at = cents, labels = f)
  axis(side = 2, at = labs.y, labels = round(labs.y,1))
}


Linf_table_stacked <- stack(Linf_table)
head(Linf_table_stacked)
names(Linf_table_stacked) <- c("Linf", "Method")
head(Linf_table_stacked) # OK


# NewBoxPlot(Data$x,Data$y)

opar <- par() 

# NewBoxPlot.ylim.Linf( Linf_table_stacked$Linf,Linf_table_stacked$Method, 
#            ylim = c( quantile(Linf_table_stacked$Linf, 0.02), 
#                     quantile(Linf_table_stacked$Linf, 0.997 )))
# 

 NewBoxPlot.ylim.Linf( Linf_table_stacked$Linf,Linf_table_stacked$Method, 
            ylim = c( 0, 400 ))
 
library(vioplot)

vioplot::vioplot(Linf_table_stacked$Linf ~ Linf_table_stacked$Method, 
                 ylim = c( 0, 400 ), col = "lightblue",
                 cex.axis = 0.8, cex.names=0.7,  las=2)
                 
         
# plots for K ----------------------------- 

# 95% quantile boxplot:


NewBoxPlot.ylim.K <- function(x,y, ylim = ylim){
  y <- as.character(y)
  nbox <- length(unique(y))
  f <- unique(y)
  lims <- c(0 - .5, nbox -.5)
  cents <- 0:(nbox-1)
  labs.y <- seq( 0 , max (x) , by = 0.25)# y labels
  
  #labs.y <- seq( 40 , quantile (x, 0.9999) , by = 20)# y labels
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '',xlim = lims, ylim = ylim)
  for (i in 1:nbox){
    xi <- x[y==f[i]]
    cent = cents[i]
    b <- quantile(xi, c(.025,.5,.975)) # quantiles
    mima <- range(xi)                        # max and min
    #rect(cent - .2, b[3], cent + .2,  border = 'firebrick')
    # points(cent + c(-.2,.2), c(b[4], b[4]), type = 'l', col = 'firebrick', lwd = 3)
    #  arrows(cent, b[5], cent, mima[2], col = 'firebrick', angle = 90)
    #arrows(cent, b[3], cent, mima[1], col = 'firebrick', angle = 90)
    points(rep(cent, 1), b[c(2)], pch = '-', cex = 1.5, col = 'navyblue')
    points(rep(cent, 2), b[c(1,3)], pch = '-', cex = 1.5, col = 'firebrick')
  }
  axis(side = 1, at = cents, labels = f)
  axis(side = 2, at = labs.y, labels = round(labs.y,1))
}


K_table_stacked <- stack(K_table)
head(K_table_stacked)
names(K_table_stacked) <- c("K", "Method")
head(K_table_stacked) # OK


# NewBoxPlot(Data$x,Data$y)

opar <- par() 

# NewBoxPlot.ylim.Linf( Linf_table_stacked$Linf,Linf_table_stacked$Method, 
#            ylim = c( quantile(Linf_table_stacked$Linf, 0.02), 
#                     quantile(Linf_table_stacked$Linf, 0.997)))
# 

NewBoxPlot.ylim.K( K_table_stacked$K,K_table_stacked$Method, 
                      ylim = c( 0, 1.5 ))

library(vioplot)

vioplot::vioplot(K_table_stacked$K ~ K_table_stacked$Method, 
                 ylim = c( 0, 1.5 ),  col = "lightblue" ,
                 cex.axis = 0.8, cex.names=0.7,  las=2)


# plots for Phi' ----------------------------- 

# 95% quantile boxplot:


NewBoxPlot.ylim.Phi <- function(x,y, ylim = ylim){
  y <- as.character(y)
  nbox <- length(unique(y))
  f <- unique(y)
  lims <- c(0 - .5, nbox -.5)
  cents <- 0:(nbox-1)
  labs.y <- seq( 2.2 , 3.4 , by = 0.2)# y labels
  
  #labs.y <- seq( 40 , quantile (x, 0.9999) , by = 20)# y labels
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '',xlim = lims, ylim = ylim)
  for (i in 1:nbox){
    xi <- x[y==f[i]]
    cent = cents[i]
    b <- quantile(xi, c(.025,.5,.975)) # quantiles
    mima <- range(xi)                        # max and min
    #rect(cent - .2, b[3], cent + .2,  border = 'firebrick')
    # points(cent + c(-.2,.2), c(b[4], b[4]), type = 'l', col = 'firebrick', lwd = 3)
    #  arrows(cent, b[5], cent, mima[2], col = 'firebrick', angle = 90)
    #arrows(cent, b[3], cent, mima[1], col = 'firebrick', angle = 90)
    points(rep(cent, 1), b[c(2)], pch = '-', cex = 1.5, col = 'navyblue')
    points(rep(cent, 2), b[c(1,3)], pch = '-', cex = 1.5, col = 'firebrick')
  }
  axis(side = 1, at = cents, labels = f)
  axis(side = 2, at = labs.y, labels = round(labs.y,1))
}

Phi_table_stacked <- stack(Phi_table)
head(Phi_table_stacked)
names(Phi_table_stacked) <- c("Phi", "Method")
head(Phi_table_stacked) # OK

# NewBoxPlot(Data$x,Data$y)

opar <- par() 

# NewBoxPlot.ylim.Linf( Linf_table_stacked$Linf,Linf_table_stacked$Method, 
#            ylim = c( quantile(Linf_table_stacked$Linf, 0.02), 
#                     quantile(Linf_table_stacked$Linf, 0.997 )))
# 

NewBoxPlot.ylim.Phi( Phi_table_stacked$Phi ,Phi_table_stacked$Method, 
                   ylim = c( 2.2 , 3.4))

library(vioplot)

vioplot::vioplot(Phi_table_stacked$Phi ~ Phi_table_stacked$Method, 
                 ylim = c(2.2 , 3.4),  col = "lightblue" ,
                 cex.axis = 0.8, cex.names=0.7,  las=2)







#########
########
#######
###
# III. 2. Barplots of 95% Confidence interval widths 
# (i.e., 95% quantiles of the posteriors, 
#  or graph comparing the precision)

# read data -----------
path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

table.1a.MA3_Linfmin50.res_median_CI <- read.table( file = "table.1a.MA3_Linfmin50.res_median_CI.csv")
head(table.1a.MA3_Linfmin50.res_median_CI)

table.1b.MA5_Linfmin50.res_median_CI <- read.table( file = "table.1b.MA5_Linfmin50.res_median_CI.csv")
head(table.1b.MA5_Linfmin50.res_median_CI)

table.1c.MA5_Linfmin30.res_median_CI <- read.table( file = "table.1c.MA5_Linfmin30.res_median_CI.csv")
head(table.1c.MA5_Linfmin30.res_median_CI)

table.1d.MA3_Linfmin30.res_median_CI <- read.table( file = "table.1d.MA3_Linfmin30.res_median_CI.csv")
head(table.1d.MA3_Linfmin30.res_median_CI)

table.1e.Tag_recap1000 <- read.table( file = "table.1e.Tag_recap1000.csv")
head(table.1e.Tag_recap1000)

table.1f.S_oto_1000 <- read.table( file = "table.1f.S_oto_1000.csv")
head(table.1f.S_oto_1000)

table.1g.N_oto_1000 <- read.table( file = "table.1g.N_oto_1000.csv")
head(table.1g.N_oto_1000)




Linf_CIWs_table <- data.frame (table.1a.MA3_Linfmin50.res_median_CI$CIW[1],
                              table.1b.MA5_Linfmin50.res_median_CI$CIW[1],
                              table.1d.MA3_Linfmin30.res_median_CI$CIW[1] ,          
                              table.1c.MA5_Linfmin30.res_median_CI$CIW[1],
                              table.1e.Tag_recap1000$CIW[1],
                              table.1f.S_oto_1000$CIW[1]     ,
                              table.1g.N_oto_1000$CIW[1]) 


K_CIWs_table <- data.frame (table.1a.MA3_Linfmin50.res_median_CI$CIW[2],
                            table.1b.MA5_Linfmin50.res_median_CI$CIW[2],
                            table.1d.MA3_Linfmin30.res_median_CI$CIW[2] ,          
                            table.1c.MA5_Linfmin30.res_median_CI$CIW[2],
                            table.1e.Tag_recap1000$CIW[2],
                            table.1f.S_oto_1000$CIW[2]   ,
                            table.1g.N_oto_1000$CIW[2] ) 


Phi_CIWs_table <- data.frame (table.1a.MA3_Linfmin50.res_median_CI$CIW[3],
                            table.1b.MA5_Linfmin50.res_median_CI$CIW[3],
                            table.1d.MA3_Linfmin30.res_median_CI$CIW[3] ,          
                            table.1c.MA5_Linfmin30.res_median_CI$CIW[3],table.1e.Tag_recap1000$CIW[3],
                            table.1f.S_oto_1000$CIW[3],
                            table.1g.N_oto_1000$CIW[3]) 


names( Linf_CIWs_table) <- names(Linf_table) 
names( K_CIWs_table) <- names(Linf_table) 
names( Phi_CIWs_table) <- names(Linf_table) 

Linf_CIWs_table_stacked <- stack(Linf_CIWs_table)
barplot( Linf_CIWs_table_stacked$values ~ Linf_CIWs_table_stacked$ind,
         col = "lightblue",  main = "Precison (CIW) of methods for Linf ")
barplot( Linf_CIWs_table_stacked$values ~ Linf_CIWs_table_stacked$ind,
         col = "lightblue",  main = "Precison (CIW) of methods for Linf ",
         ylab = "95%Confid. Interval width for Linf (cm)",
         cex.axis = 0.8, cex.names=0.7,  las=2)

K_CIWs_table_stacked <- stack(K_CIWs_table)
barplot( K_CIWs_table_stacked$values ~ K_CIWs_table_stacked$ind,
         col = "lightblue",  main = "Precison (CIW) of methods for K ")
barplot( K_CIWs_table_stacked$values ~ K_CIWs_table_stacked$ind,
         col = "lightblue",  main = "Precison (CIW) of methods for K ",
ylab = "95%Confid. Interval width for K (y-1)",
cex.axis = 0.8, cex.names=0.7,  las=2)




Phi_CIWs_table_stacked <- stack(Phi_CIWs_table)
barplot( Phi_CIWs_table_stacked$values ~ Phi_CIWs_table_stacked$ind,
         col = "lightblue", main = "Precison (CIW) of methods for Phi' ")

barplot( Phi_CIWs_table_stacked$values ~ Phi_CIWs_table_stacked$ind,
         col = "lightblue", main = "Precison (CIW) of methods for Phi'",
         ylab = "95%Confid. Interval width for Phi (log cm + log y-1)",
         cex.axis = 0.8, cex.names=0.7,  las=2)

##########
#########
# FOR PAPER -----------------
# New complete barplots ( j = 7 classes) -----------
# vioplots with LFA (j = 4), T&R (j = 1), S and N Stock LAA (j = 2)


######
#####
####
##
# IV. Chapter IV: S vs N stock, males vs females (otoliths only) ------------------------
# S vs N stock (otoliths only)
# males vs females (otoliths only)

# IV.1 S vs N stock (otoliths only) ----------------------------------

# read data ----------------

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

oto.N.all.1000results <- read.csv("oto-N-all-1000results.csv")
#   View(oto.N.all.1000results)

oto.S.all.1000results <- read.csv("oto-S-all-1000results (1).csv")
#   View(`oto.S.all.1000results

# IV.1 S vs N stock compare medians  -----------


# Medians - Test for difference in Medians, based on the bootstrap posteriors---------
# use the "bootstrapped.median.test"
# pairwise.bootstrapped.median.test --------------------

# #  The Bootstrapped two-sample test, based on posteriors 
# Ralf Schwamborn
# March 20, 2022
# Citation: 

# Schwamborn, R., 2022. The bootstrapped two-sample test. 
#〈https://rpubs.com/rschwamborn/880212〉 (accessed 22 March, 2022).

#Source:
# https://rpubs.com/rschwamborn/880212

#see also:
# Schwamborn, R., 2019. The interquantile range test. 
#〈http://rpubs.com/rschwamborn/515419〉 (accessed 25 July 2019).


# define the function "boot.p.two.sample.medians.diff.posteriors_B()" ------------------


boot.p.two.sample.medians.diff.posteriors_B <- function (boot.post.B, boot.post.A) {
  
  # difference between posteriors
  
  boot.statisticsA <- boot.post.A # posterior for sample A
  boot.statisticsB <- boot.post.B # posterior for sample B
  
  if ( median(boot.post.B) < median(boot.post.A) ) {
    # data set "B", has larger values
    
    boot.statisticsDiff <- boot.statisticsA - boot.statisticsB
    N.below.zero <- length(boot.statisticsDiff [boot.statisticsDiff<= 0]) # 7 are below 0 
    
    N.total <- length(boot.statisticsDiff ) # all 
    
    (one.sided.p_value <- N.below.zero / N.total) # one-sided "p" value
    
    (two.sided.p_value <- 2* (N.below.zero / N.total)) # two-sided "p" value
    
    ans <-   data.frame( text = "two.sided.p_value", p.value = two.sided.p_value)  # output
  }
  
  
  
  else 
  {boot.statisticsDiff <- boot.statisticsB - boot.statisticsA
  
  # calculate "p" value  from bootstrap diff.  
  
  N.below.zero <- length(boot.statisticsDiff [boot.statisticsDiff <= 0]) # n below 0 or equal 0 
  
  N.total <- length(boot.statisticsDiff ) # all 
  
  ( one.sided.p_value <- N.below.zero / N.total) # one-sided "p" value
  
  (two.sided.p_value <- 2* (N.below.zero / N.total)) # two-sided "p" value
  
  ans <-   data.frame( text = "two.sided.p_value", p.value = two.sided.p_value)  # output
  }
  
  paste ( as.list(ans) )   # output
  
}



# test the function "boot.p.two.sample.medians.diff.posteriors()" ----------

# Linf
boot.p.two.sample.medians.diff.posteriors_B (oto.S.all.1000results$Linf ,oto.N.all.1000results$Linf )

# median Linf values are not significantly different between S and N stocks,
# , two.sided.p_value = 0.71


# K
boot.p.two.sample.medians.diff.posteriors_B (oto.S.all.1000results$K ,oto.N.all.1000results$K )

# median K values  are not significantly different between S and N stocks,
# , two.sided.p_value = 0.668


# to
boot.p.two.sample.medians.diff.posteriors_B (oto.S.all.1000results$t0 ,oto.N.all.1000results$t0 )

# median t0 values  are not significantly different between S and N stocks,
# , two.sided.p_value = 0.128


# Phi'
oto.S.all.1000results$Phi_L <-  log10(oto.S.all.1000results$K) + (2 * log10(oto.S.all.1000results$Linf))
oto.N.all.1000results$Phi_L <-  log10(oto.N.all.1000results$K) + (2 * log10(oto.N.all.1000results$Linf))

boot.p.two.sample.medians.diff.posteriors_B (oto.S.all.1000results$Phi_L ,oto.N.all.1000results$Phi_L )

# median Phi'values are not significantly different between S and N stocks,
# , two.sided.p_value = 0.08

# plots S vs N stocks (otoliths)

library (vioplot)

vioplot (oto.S.all.1000results$Linf ,oto.N.all.1000results$Linf , 
         main = "Linf")

vioplot (oto.S.all.1000results$Linf ,oto.N.all.1000results$Linf , 
         main = "Linf", ylim = c(40, 160))


vioplot (oto.S.all.1000results$K,oto.N.all.1000results$K , 
         main = "K")

vioplot (oto.S.all.1000results$Phi_L ,oto.N.all.1000results$Phi_L , 
                  main = "Phi'")
   
 
  
#######
#####
###
# IV.2 Test for differences between sexes ---------------------
  
# IV.2 Males vs Females vs Hermaphrodites, compare medians (otolith LAA) -----------
  
  
  
# read data ----------------
  
path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
  
  # oto.N.all.1000results <- read.csv("oto-N-all-1000results.csv")
  # #   View(oto.N.all.1000results)
  # 
  # oto.S.all.1000results <- read.csv("oto-S-all-1000results (1).csv")
  # #   View(oto.S.all.1000results)
  
  
    oto.S.male<- read.csv("oto-S-male-1000results.csv")
  # >   View(oto.S.male.1000results)
    oto.N.male <- read.csv("oto-N-male-1000results.csv")
  # >   View(oto.N.male.1000results)
   oto.N.female <- read.csv("oto-N-female-1000results.csv")
  # >   View(oto.N.female.1000results)
    oto.S.female <- read.csv("oto-S-female-1000results.csv")
  # >   View(oto.S.female.1000results)
   
   oto.N.herm <- read.csv("oto-N-herm-1000results.csv")
   #   View(oto.N.herm.1000results)
    oto.S.herm <- read.csv("oto-S-herm-1000results.csv")
  #    View(oto.S.herm.1000results)
         
    
  
  # IV.2 Females vs Males vs Hermaphr.,  compare medians (S and N stock)  -----------
  
  
  # Medians - Test for difference in Medians, based on the bootstrap posteriors---------
  # use the "bootstrapped.median.test"
  # pairwise.bootstrapped.median.test --------------------
  
  # #  The Bootstrapped two-sample test, based on posteriors 
  # Ralf Schwamborn
  # March 20, 2022
  # Citation: 
  
  # Schwamborn, R., 2022. The bootstrapped two-sample test. 
  #〈https://rpubs.com/rschwamborn/880212〉 (accessed 22 March, 2022).
  
  #Source:
  # https://rpubs.com/rschwamborn/880212
  
  #see also:
  # Schwamborn, R., 2019. The interquantile range test. 
  #〈http://rpubs.com/rschwamborn/515419〉 (accessed 25 July 2019).
  
  
  # define the function "boot.p.two.sample.medians.diff.posteriors_B()" ------------------
  
  
  boot.p.two.sample.medians.diff.posteriors_B <- function (boot.post.B, boot.post.A) {
    
    # difference between posteriors
    
    boot.statisticsA <- boot.post.A # posterior for sample A
    boot.statisticsB <- boot.post.B # posterior for sample B
    
    if ( median(boot.post.B) < median(boot.post.A) ) {
      # data set "B", has larger values
      
      boot.statisticsDiff <- boot.statisticsA - boot.statisticsB
      N.below.zero <- length(boot.statisticsDiff [boot.statisticsDiff<= 0]) # 7 are below 0 
      
      N.total <- length(boot.statisticsDiff ) # all 
      
      ( one.sided.p_value <- N.below.zero / N.total) # one-sided "p" value
      
      (two.sided.p_value <- 2* (N.below.zero / N.total)) # two-sided "p" value
      
      ans <-   data.frame( text = "two.sided.p_value", p.value = two.sided.p_value)  # output
    }
    
    
    
    else 
    {boot.statisticsDiff <- boot.statisticsB - boot.statisticsA
    
    # calculate "p" value  from bootstrap diff.  
    
    N.below.zero <- length(boot.statisticsDiff [boot.statisticsDiff <= 0]) # n below 0 or equal 0 
    
    N.total <- length(boot.statisticsDiff ) # all 
    
    ( one.sided.p_value <- N.below.zero / N.total) # one-sided "p" value
    
    (two.sided.p_value <- 2* (N.below.zero / N.total)) # two-sided "p" value
    
    ans <-   data.frame( text = "two.sided.p_value", p.value = two.sided.p_value)  # output
    }
    
    paste ( as.list(ans) )   # output
    
  }
  
  
  
  # test the function "boot.p.two.sample.medians.diff.posteriors()" ----------
  
    # rename function for convenience (short name: boot.test.medians)
    boot.test.medians <- boot.p.two.sample.medians.diff.posteriors_B
    
# IV 2.a Males vs Females,  South and N stock -------------------------------------------
    # test for differences in medians
    
  # Linf
    boot.test.medians (oto.S.male$Linf ,oto.S.female$Linf )
  # medians are not significantly different, "two.sided.p_value" "0.718"
   
    boot.test.medians (oto.N.male$Linf ,oto.N.female$Linf )
  # medians are not significantly different, "two.sided.p_value" "0.636"
    
    
    
    # K
    boot.test.medians (oto.S.male$K ,oto.S.female$K )
    # medians are not significantly different, "two.sided.p_value" "0.772"
   
    boot.test.medians (oto.N.male$K ,oto.N.female$K )
    # medians are not significantly different, "two.sided.p_value" "0.546"
    
     
    # t0
    boot.test.medians (oto.S.male$t0 ,oto.S.female$t0 )
    # medians are not significantly different, "two.sided.p_value" "0.77"
 
    # t0
    boot.test.medians (oto.N.male$t0 ,oto.N.female$t0 )
    # medians are not significantly different, "two.sided.p_value" "0.432"
    
    summary(oto.N.female$t0)
    summary(oto.N.male$t0)
    
       
    # Phi'
    
    oto.S.male$Phi_L <-  log10(oto.S.male$K) + (2 * log10(oto.S.male$Linf))
    oto.S.female$Phi_L <-  log10(oto.S.female$K) + (2 * log10(oto.S.female$Linf))
   
    oto.N.male$Phi_L <-  log10(oto.N.male$K) + (2 * log10(oto.N.male$Linf))
    oto.N.female$Phi_L <-  log10(oto.N.female$K) + (2 * log10(oto.N.female$Linf))
    
    oto.S.herm$Phi_L <-  log10(oto.S.herm$K) + (2 * log10(oto.S.herm$Linf))
    oto.N.herm$Phi_L <-  log10(oto.N.herm$K) + (2 * log10(oto.N.herm$Linf))
    
    
    
    boot.test.medians (oto.S.male$Phi_L ,oto.S.female$Phi_L )
    # medians are not significantly different, "two.sided.p_value" "0.824"
    
    boot.test.medians (oto.N.male$Phi_L ,oto.N.female$Phi_L )
    # medians are not significantly different, "two.sided.p_value" "0.326"

    
    # IV 2.b Males vs Hermaphrodites --------------
    
    # Males vs Hermaphr., South and N stock -------------------------------------------
    # test for differences in medians
    
    # Linf
    boot.test.medians (oto.S.male$Linf ,oto.S.herm$Linf )
    # medians are not significantly different, p > 0.05
    
    boot.test.medians (oto.N.male$Linf ,oto.N.herm$Linf )
    # medians are not significantly different, p > 0.05
    
    
    
    # K
    boot.test.medians (oto.S.male$K ,oto.S.herm$K )
    # medians are not significantly different, ", p > 0.05
    
    boot.test.medians (oto.N.male$K ,oto.N.herm$K )
    # medians are not significantly different, , p > 0.05
    
    
    # t0
    boot.test.medians (oto.S.male$t0 ,oto.S.herm$t0 )
    # medians are not significantly different, , p > 0.05
    
    # t0
    boot.test.medians (oto.N.male$t0 ,oto.N.herm$t0 )
    # medians are not significantly different,  p > 0.05
    
    summary(oto.N.herm$t0)
    summary(oto.N.male$t0)
    
    
    # Phi'
    
    oto.S.male$Phi_L <-  log10(oto.S.male$K) + (2 * log10(oto.S.male$Linf))
    oto.S.herm$Phi_L <-  log10(oto.S.herm$K) + (2 * log10(oto.S.herm$Linf))
    
    oto.N.male$Phi_L <-  log10(oto.N.male$K) + (2 * log10(oto.N.male$Linf))
    oto.N.herm$Phi_L <-  log10(oto.N.herm$K) + (2 * log10(oto.N.herm$Linf))
    
    oto.S.herm$Phi_L <-  log10(oto.S.herm$K) + (2 * log10(oto.S.herm$Linf))
    oto.N.herm$Phi_L <-  log10(oto.N.herm$K) + (2 * log10(oto.N.herm$Linf))
    
    
    
    boot.test.medians (oto.S.male$Phi_L ,oto.S.herm$Phi_L )
    # medians are not significantly different, p > 0.05
    
    boot.test.medians (oto.N.male$Phi_L ,oto.N.herm$Phi_L )
    # medians are not significantly different,  p > 0.05
    
    
    
    # IV 2.c Females vs Hermaphrodites --------------
    
    # Females vs Hermaphr., South and N stock -------------------------------------------
    # test for differences in medians
    
    # Linf
    boot.test.medians (oto.S.female$Linf ,oto.S.herm$Linf )
    # medians are not significantly different, p > 0.05
    
    boot.test.medians (oto.N.female$Linf ,oto.N.herm$Linf )
    # medians are not significantly different, p > 0.05
    
    # K
    boot.test.medians (oto.S.female$K ,oto.S.herm$K )
    # medians are not significantly different, ", p > 0.05
    
    boot.test.medians (oto.N.female$K ,oto.N.herm$K )
    # medians are not significantly different, , p > 0.05
    
    
    # t0
    boot.test.medians (oto.S.female$t0 ,oto.S.herm$t0 )
    # medians are not significantly different, , p > 0.05
    
    # t0
    boot.test.medians (oto.N.female$t0 ,oto.N.herm$t0 )
    # medians are not significantly different,  p > 0.05
    
    summary(oto.N.herm$t0)
    summary(oto.N.female$t0)
    
    
    # Phi'
    
    boot.test.medians (oto.S.female$Phi_L ,oto.S.herm$Phi_L )
    # medians are not significantly different, p > 0.05
    
    boot.test.medians (oto.N.female$Phi_L ,oto.N.herm$Phi_L )
    # medians are not significantly different,  p > 0.05
    
    
    
    
    # V. Chapter V: simple LFD plots for paper, with "optimum curve" ---------------------------
    
    # LFD plots (S tock) --------------------------

    
    # Data: 
    # Fish lengths (fork length) from 2004-08-25 to 2009-12-30 
    # (approx. five and a half years) 
    # Length range: 12 cm to 79 cm
    # N = 1399 indiv.
    
    
    #V. 1. Initialization ----------------------------------------------------------
    
    # V.1.1. required packages =====
    library(parallel)
    library(TropFishR)
    library(fishboot)
    library(ks)
    
    # clean memory
    gc()  
    gc(reset=T)
    rm(list = ls())
    
    
    # V.1.2. Load data =====
    
    
    path = dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(path)
    
    lfq2 <- read.csv("stbrasLF-trunc.csv")
    
    # View(lfq2)
    
    
    
    lfq2
    names(lfq2)
    lfq2$DATE <- as.Date(lfq2$DATE, format = "%d/%m/%Y") #Date format needs to be input like that
    lfq2
    
    lfq2new <- lfqCreate(data = lfq2, Lname = "Length_CM", Dname = "DATE", Fname = "Frequency")
    plot(lfq2new, Fname = "catch")
    lfq2new
    
    summary(lfq2$DATE)
    
    summary(lfq2$Length_CM)
    
    length(lfq2$Length_CM)
    
    
    # ---
    lfq <- lfq2new
    lfq
    
    ## 
    
    
    # 1.3 plot raw and restructured LFQ data --------------------------
    
    # plot raw and restructured LFQ data bs = 2 --------------------
    
    # adjust bin size (on the histogram) #Bin size 2!
    synLFQ7a <- lfqModify(lfq, bin_size = 2) #mid-sizes are already at 0.5
    
    lfqbin <- lfqRestructure(synLFQ7a, MA = 3, addl.sqrt = TRUE) #MA = moving average sample size 
    opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
    plot(lfqbin, Fname = "catch", date.axis = "modern")
    plot(lfqbin, Fname = "rcounts", date.axis = "modern")
    par(opar)
    
    
    # plot raw and restructured LFQ data bs = 2, Ma = 3 ----------------
    
    # adjust bin size (on the histogram) #Bin size 2!
    synLFQ7a <- lfqModify(lfq, bin_size = 2) #mid-sizes are already at 0.5
    
    lfqbin <- lfqRestructure(synLFQ7a, MA = 3, addl.sqrt = TRUE) #MA = moving average sample size 
    opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
    plot(lfqbin, Fname = "catch", date.axis = "modern")
    plot(lfqbin, Fname = "rcounts", date.axis = "modern", main = "MA3")
    par(opar)
    
    # plot raw and restructured LFQ data bs = 2, Ma = 5 ----------------
    
    # adjust bin size (on the histogram) #Bin size 2!
    synLFQ7a <- lfqModify(lfq, bin_size = 2) #mid-sizes are already at 0.5
    
    lfqbin <- lfqRestructure(synLFQ7a, MA = 5, addl.sqrt = TRUE) #MA = moving average sample size 
    opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
    plot(lfqbin, Fname = "catch", date.axis = "modern")
    plot(lfqbin, Fname = "rcounts", date.axis = "modern", main = "MA5")
    par(opar)
    
    
    
    
    # V.2. ELEFAN_GA (empty res object) ---------------------------------------------------------------
    # V.2.1. search settings GA =====
    # 
    
    # SUPER-SUPERFAST parameters (for creating an empty "res" objetc only)
    MA <- 5 #5 
    low_par <- list(Linf = 50, K = 0.01, t_anchor = 0, C = 0, ts = 0)
    up_par <- list(Linf = 160, K = 2, t_anchor = 1, C = 1, ts = 1)
    seasonalised <- TRUE
    popSize <- 12
    maxiter <- 5 #30 but start with 10
    run <- 4
    pmutation <- 0.2
    #nresamp <- 4 #Number of resamplings; 20 but start with 5
   
     
    
    # V. 2.2. test settings GA =====
    # (check for sufficiently long plateau of Rn scores)
    
    t1 <- Sys.time()
    set.seed(1)
    res0 <- ELEFAN_GA(synLFQ7a, MA = MA, seasonalised = seasonalised, #used to be x=synLFQ7a
                      up_par = up_par, low_par = low_par,
                      popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation
    )
    res0$Rn_max
    unlist(res0$par)
    
    
    
    
    plot(res0) # plot resulting growth curve (test object)
    plot(res0, Fname = "rcounts") # plot resulting growth curve against restructured scores
    res0
    t2 <- Sys.time()
    t2 - t1
    
    
    
    #V.3 FINAL LFD plots for paper -----------------
   
    #set wd
    
    path = dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(path)
    
    
     
    # V.3a:  MA3, bs 2, Linfmin30 ------------------
  
    MA <- 3 #3
    
    
    # MA3, Linfmin30
    # Linf              K              t_anchor             C          
    # Min.   : 61.14   Min.   :0.01743   Min.   :0.06839   Min.   :0.01911  
    # 1st Qu.: 75.43   1st Qu.:0.07327   1st Qu.:0.38110   1st Qu.:0.40199  
    # Median : 83.89   Median :0.11483   Median :0.50126   Median :0.53389  
    # Mean   : 88.25   Mean   :0.13473   Mean   :0.50223   Mean   :0.53257  
    # 3rd Qu.: 95.34   3rd Qu.:0.17200   3rd Qu.:0.61981   3rd Qu.:0.66038  
    # Max.   :180.69   Max.   :0.75513   Max.   :0.96651   Max.   :0.98990  
    
    # ts                phiL      
    # Min.   :0.001587   Min.   :2.679  
    # 1st Qu.:0.312684   1st Qu.:2.839  
    # Median :0.468051   Median :2.903  
    # Mean   :0.466835   Mean   :2.929  
    # 3rd Qu.:0.597131   3rd Qu.:3.018  
    # Max.   :0.990533   Max.   :3.525 
    
    
  
      
    t1 <- Sys.time()
    set.seed(1)
    resMA3 <- ELEFAN_GA(synLFQ7a, MA = MA, seasonalised = seasonalised, #used to be x=synLFQ7a
                      up_par = up_par, low_par = low_par,
                      popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation
    )
    resMA3$Rn_max
    unlist(resMA3$par)
    
    
    # insert "optimum values into res object" --------------------- 
    
    
# insert Median of 1000 bootstrap runs, bs2, MA3, Linf 30 to 200 cm ---------
# length(res_MA3_Linf30to200$Linf) # 1000 lines
   
    
    
    
    resMA3.OK <- resMA3
    
        
    # Median Linf: 83.89
    resMA3.OK$par$Linf <- 83.89
  
      
    # Median K: 0.11483
    resMA3.OK$par$K <- 0.11483
      
       
    # Mediant_anchor :0.50126
    resMA3.OK$par$t_anchor <- 0.50126
        
       
    # Median C:0.53389  
    resMA3.OK$par$C <- 0.53389 
        
        
    # Median ts:0.468051 
    resMA3.OK$par$ts <- 0.468051 
     
        
    # Median Phi:2.903 
    resMA3.OK$par$phiL <- 2.903 
        
    
      
    # adjust bin size (on the histogram) #Bin size 2!
    synLFQ7a <- lfqModify(lfq, bin_size = 2) #mid-sizes are already at 0.5
    
    lfqbin <- lfqRestructure(synLFQ7a, MA = 3, addl.sqrt = TRUE) #MA = moving average sample size 
    opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
    plot(lfqbin, Fname = "catch", date.axis = "modern")
    plot(resMA3.OK, Fname = "rcounts", date.axis = "modern", 
         main = "bs = 2cm, MA = 3, Linf_min = 30cm")
    par(opar)
    
     # filename: "LFDs_MA3_wCurve.svg"
    
    
    
    
    # V.3b: MA5, bs 2 ------------------
    
        
    
    MA <- 5 #3
    
    
    
    t1 <- Sys.time()
    set.seed(1)
    resMA5 <- ELEFAN_GA(synLFQ7a, MA = MA, seasonalised = seasonalised, #used to be x=synLFQ7a
                        up_par = up_par, low_par = low_par,
                        popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation
    )
    resMA5$Rn_max
    unlist(resMA5$par)
    
    
    # insert "optimum values into res object" --------------------- 
    
    
    
    # insert Median of 1000 bootstrap runs, bs2, MA5, Linf 30 to 200 cm ---------
    # length(res_MA5_Linf30to200$Linf) # 1000 lines
    
    # MA5, Linf 30 to 200 cm
    # res_EXP3_S_slowLinf30_200_n1000
    
    # Linf              K              t_anchor             C          
    # Min.   : 39.98   Min.   :0.02673   Min.   :0.02715   Min.   :0.08297  
    # 1st Qu.: 44.70   1st Qu.:0.14485   1st Qu.:0.40597   1st Qu.:0.43408  
    # Median : 51.30   Median :0.49404   Median :0.51023   Median :0.53562  
    # Mean   : 69.02   Mean   :0.46342   Mean   :0.50638   Mean   :0.54777  
    # 3rd Qu.: 89.78   3rd Qu.:0.74217   3rd Qu.:0.61605   3rd Qu.:0.66748  
    # Max.   :187.63   Max.   :1.38949   Max.   :0.95265   Max.   :0.99591  
    # ts               phiL      
    # Min.   :0.01239   Min.   :2.748  
    # 1st Qu.:0.37517   1st Qu.:3.042  
    # Median :0.45977   Median :3.143  
    # Mean   :0.45915   Mean   :3.116  
    # 3rd Qu.:0.54166   3rd Qu.:3.214  
    # Max.   :0.98729   Max.   :3.379 
    
    
    

    
    resMA5.OK <- resMA5
    
    # Median : 51.30
    # Median Linf: 
    resMA5.OK$par$Linf <-   51.30
    
    #Median K: 
    # Median :0.49404
    resMA5.OK$par$K <- 0.49404
    
    #Median t_anchor: 
    #Median :0.51023
    resMA5.OK$par$t_anchor <-0.51023
    
    # Median C: 
    #Median :0.53562  
    resMA5.OK$par$C <-    0.53562  
    
    # Median ts : 
    # Median :0.45977
    resMA5.OK$par$ts <-    0.45977
    
    
    # Median Phi':
    # Median :3.143  
      
    resMA5.OK$par$phiL <-   3.143 
     
    
    
    
    
    
    # adjust bin size (on the histogram) #Bin size 2!
    synLFQ7a <- lfqModify(lfq, bin_size = 2) #mid-sizes are already at 0.5
    
    lfqbin <- lfqRestructure(synLFQ7a, MA = 5, addl.sqrt = TRUE) #MA = moving average sample size 
    opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
    plot(lfqbin, Fname = "catch", date.axis = "modern")
    plot(resMA5.OK, Fname = "rcounts", date.axis = "modern", 
         main = "bs = 2cm, MA = 5, Linfin = 30")
    par(opar)
    
    # filename: "LFDs_MA5_wCurve.svg"
    
    
    #############
    ############
    ###########
    # VI. Chapter VI: simple T&R plots for paper,  --------------------------------------
    #      with "curve swarms", increments  and "optimum curve" --------------------
    
    
    # clean memory
    gc()  
    gc(reset=T)
    rm(list = ls()) 
    
    # install.packages("fishmethods")
    # install.packages("devtools")
    # library(devtools)
    # install_github("rschwamborn/fishboot")
    # install.packages("devtools")
    
     
    
    # required packages =====
    library(parallel)
    library(TropFishR)
    library(fishboot)
#   library(ks)
    library(fishmethods)
    
    opar <- par()   # save default graphical parameters (for plotting) 
    
    
    
    # VI.1 : read in the T&R data (posterior distributions)
  
    #  T&R ,  Tag-recapture (TR1b data set) ---------------------------
    # Dataset "TR1b",    2020-2022 data were added to th original dataset
    # and removed two fish -13 and - 7.5 growth
    # N = 80
    
    path = dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(path)
    
    #  VI.1a read data -------------------
    
    res_E_1000StoutputB1_9_from_Margit <- read.csv("1000StoutputB1_9_from_Margit.csv", sep=",")
    
    res_E_Tag_recap1000 <- res_E_1000StoutputB1_9_from_Margit
    
    # Comment: Ctrl + Shift + C
    
    
    # 95% CI and CIW, CIW% ---------------------------------
    
    #define functions for lower and upper 95% CI limits ----------
    quantile_lowlim <- function(x) {quantile(x,  0.025)}
    quantile_upperlim <- function(x) {quantile(x,  0.975)}
    
    head(res_E_Tag_recap1000)
    
    attach(res_E_Tag_recap1000)
    
    res_E_Tag_recap1000$Phi_L <- log10(res_E_Tag_recap1000$K) + (2 * log10(res_E_Tag_recap1000$Linf))
    #0.103  68.1 2.679131
    log10(0.103) + (2 * log10(68.1)) # OK, 2.679131
    
    
    head(res_E_Tag_recap1000)# OK
    
    # reorganize, Linf first
    res_E_Tag_recap1000 <- res_E_Tag_recap1000[,c(2,1,3)]
    
    head(res_E_Tag_recap1000)# OK
    
    
    
    # #  VI.1b Create "res0" object to be filed with T& R posteriors -------------------
    # Near-EMPTY res object ELEFAN_GA (creates  "res0" object with 3 runs) ---------------------------------------------------------------
   
    # search settings GA =====
    # SUPER-SUPERFAST parameters (for creating a near-empty "res0" object only)
 
    #  GA with test settings (creates near-empty "res0" object with 3 runs)  =====
    
    # SUPER-Fast settings
    MA <- 5 #5 is best!
    low_par <- list(Linf = 50, K = 0.01, t_anchor = 0, C = 0, ts = 0)
    up_par <- list(Linf = 200, K = 2, t_anchor = 1, C = 1, ts = 1)
    seasonalised <- FALSE
    popSize <- 12
    maxiter <- 6 #30 but start with 10
    run <- 3
    pmutation <- 0.2
    nresamp <- 3 #Number of resamplings; 20 but start with 5
    
    lfq2 <- read.csv("stbrasLF-trunc.csv")    #View(steenbras1trunc)
    
    lfq2
    names(lfq2)
    lfq2$DATE <- as.Date(lfq2$DATE, format = "%d/%m/%Y") #Date format needs to be input like that
    lfq2
    
    lfq2new <- lfqCreate(data = lfq2, Lname = "Length_CM", Dname = "DATE", Fname = "Frequency")
    plot(lfq2new, Fname = "catch")
    lfq2new
    
    # ---
    lfq <- lfq2new
    lfq
    
    # adjust bin size (on the histogram) #Bin size 2!
    synLFQ7a <- lfqModify(lfq, bin_size = 2) #mid-sizes are already at 0.5
    
    
    
    #  test settings GA =====
    # (check for sufficiently long plateau of Rn scores)
    set.seed(1)
    res0 <- ELEFAN_GA(synLFQ7a, MA = MA, seasonalised = seasonalised, #used to be x=synLFQ7a
                      up_par = up_par, low_par = low_par,
                      popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation
    )
    res0$Rn_max
    unlist(res0$par)
    plot(res0) # plot resulting growth curve
    plot(res0, Fname = "rcounts") # plot resulting growth curve against restructured scores
    res0
    
    
    # full bootstrap GA=====
    t1 <- Sys.time()
    res0 <- ELEFAN_GA_boot(lfq=synLFQ7a, MA = MA, seasonalised = seasonalised,
                           up_par = up_par, low_par = low_par,
                           popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
                           nresamp = nresamp, parallel = TRUE, no_cores = detectCores(),
                           seed = 1, resample = TRUE
    )
    t2 <- Sys.time()
    t2 - t1
    res0$boot
    head(res0$seed)
    
   
    
    
    
    #  VI.1c insert the T & R posteriors into the empty "res0" object ---------------
   
    # create new "res.T_R" object  for T & R posteriors --------------
    #  T&R ,  Tag-recapture (TR1b data set) ---------------------------
    # Dataset "TR1b",    2020-2022 data were added to th original dataset
    # and removed two fish -13 and - 7.5 growth
    # N = 80
    
    head(res0$bootRaw)# OK
    
    res.T_R <- res0
   
    head(res_E_Tag_recap1000)# OK
    
    
    res.T_R$bootRaw <- res_E_Tag_recap1000
    
    head(res.T_R$bootRaw)# OK
    
    # insert t_anchor (zero values ) for plotting
    
    t_anchor.Z <- rep (0, length(res.T_R$bootRaw$Linf) )
    
    res.T_R$bootRaw$t_anchor <-  t_anchor.Z
    
    head(res.T_R$bootRaw)# OK
    
    
    
    
    # VI.2 : plot the optimum curve with increments
    
    # read in T & R data (80 Increments)
    # T&R ,  Tag-recapture (TR1b data set) ---------------------------
    # Dataset "TR1b",    2020-2022 data were added to the original dataset
    # and removed two fish -13 and - 7.5 growth
    # N = 80
    
    
    steenTR1b <- read.csv("steenTR1b.csv")
    #   View(steenTR1b)
    
    length(steenTR1b$T1) # 80 growth increments, OK
   
    summary(steenTR1b) # Max. Length  :61.50 
    
     
    # calculate "optimum" (median of posteriors) growth curve
    
    summary(res.T_R$bootRaw)
    
    # Linf               K              Phi_L          t_anchor
    # Min.   :  38.80   Min.   :0.0010   Min.   :1.994   Min.   :0  
    # 1st Qu.:  45.50   1st Qu.:0.0900   1st Qu.:2.615   1st Qu.:0  
    # Median :  49.70   Median :0.1850   Median :2.717   Median :0  
    # Mean   :  74.53   Mean   :0.1934   Mean   :2.709   Mean   :0  
    # 3rd Qu.:  67.40   3rd Qu.:0.2790   3rd Qu.:2.806   3rd Qu.:0  
    # Max.   :3710.70   Max.   :0.6300   Max.   :4.139   Max.   :0 
    
    # median values:
    # Median Linf:  49.70   Median K:0.1850     t_anchor:0   
    
    
    
    
    opar4 <- par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(0,0,0,0), mgp = c(2,0.5,0), 
                tcl = -0.25,
                cex = 1)
    
    # problem Linf below Lmax ... (!!!)
    # growthTraject(0.185, 49.70,
    #               lentag=steenTR1b$L1, lenrec=steenTR1b$L2,
    #               timelib=c(steenTR1b$T2-steenTR1b$T1), 
    #               ylim = c(20,65), xlim = c(2,25), 
    #               main = "Growth increments & one fitted VBGF curve", ylab = "Fork length (cm)", xlab = "Relative age (years)" ) #K, Linf
 
    # three plots:
    

    # 1. Linf = 49.70 cm, no increments at all, "empty plot"
   
    library(TropFishR)
    t <- seq(0,24,0.1)
    Lt <- VBGF(list(Linf=49.70, K=0.185, t0= 0),t=t)
    plot(t, Lt, t="l", ylim = c(20,65), xlim = c(2,25), col = "red")
   
    
    # 2. Linf = 49.70 cm, small ones only (not plotting ALL increments, just the ones below Linf)   
    steenTR1b.small <- subset(steenTR1b, L1 < 49.70  ) 
    
    
      growthTraject(0.185, 49.70,
                  lentag=steenTR1b.small$L1, lenrec=steenTR1b.small$L2,
                  timelib=c(steenTR1b.small$T2-steenTR1b.small$T1), 
                  ylim = c(20,65), xlim = c(2,25), 
                  main = "Growth increments & one fitted VBGF curve", ylab = "Fork length (cm)", xlab = "Relative age (years)" ) #K, Linf
      
      lines(t, Lt, t="l", ylim = c(20,65), xlim = c(2,25), col = "red")
      
   
     # 3. Linf = 62 cm, ALL increments (just for plotting ALL increments)   
    growthTraject(0.185, 62,
                  lentag=steenTR1b$L1, lenrec=steenTR1b$L2,
                  timelib=c(steenTR1b$T2-steenTR1b$T1), 
                  ylim = c(20,65), xlim = c(2,25), 
                  main = "Growth increments & one fitted VBGF curve", ylab = "Fork length (cm)", xlab = "Relative age (years)" ) #K, Linf
    
    
    # VI.2 : plot the  curve swarm 
    
    # help(package = "fishboot")
    # ? vbgfCI_time
    # vbgfCI_time(res, CI = 95, agemax = NULL, add_legend = TRUE,
    #             add_max_dens_legend = TRUE, xlab = "Relative time",
    #             ylab = "Length", perm.col = adjustcolor("grey50", 0.1),
    #             perm.lwd = 1, ci.col = 1, ci.lty = 2, ci.lwd = 1, maxd.col = 1,
    #             maxd.lty = 1, maxd.lwd = 2)
    # 
   
    
    library (ks)
    
    res.T_Rexp <- res.T_R
    res.T_Rexp$bootRaw$t_anchor <- runif(1017, -0.01, 0.01) 
    summary(res.T_Rexp$bootRaw)
    
    
    CIinfo <- vbgfCI_time(
      res = res.T_Rexp,
      agemax = 24, CI = 95,
      perm.col = adjustcolor("grey50",0.2)
    )
  
#    insert median optimum curve (median of posteriors), red line
      
    lines(t, Lt, t="l", col = "red")
    
    
     # VI.2 : plot the  curve swarm with optimum curve, with increments
    
    growthTraject(0.185, 49.70,
                  lentag=steenTR1b.small$L1, lenrec=steenTR1b.small$L2,
                  timelib=c(steenTR1b.small$T2-steenTR1b.small$T1), 
                  ylim = c(0,90), xlim = c(0,24), 
                  main = "Growth increments & one fitted VBGF curve", ylab = "Fork length (cm)", xlab = "Relative age (years)" ) #K, Linf
    
    lines(t, Lt, t="l",  col = "red")
    
          
    # 3. Linf = 62 cm, ALL increments (just for plotting ALL increments)   
    growthTraject(0.185, 62,
                  lentag=steenTR1b$L1, lenrec=steenTR1b$L2,
                  timelib=c(steenTR1b$T2-steenTR1b$T1), 
                  ylim = c(0,90), xlim = c(0,24), 
                  main = "Growth increments & one fitted VBGF curve", ylab = "Fork length (cm)", xlab = "Relative age (years)" ) #K, Linf
    
    
    
# VII. Chapter VII: simple LAA (otolith reading data) plots for paper,  --------------------------------------
#      with "curve swarms", ages by sex and "optimum curve" --------------------
  
      
    # clean memory
    gc()  
    gc(reset=T)
    rm(list = ls()) 
    
    # install.packages("fishmethods")
    # install.packages("devtools")
    # library(devtools)
    # install_github("rschwamborn/fishboot")
    # install.packages("devtools")
    
    
    
    # required packages =====
    library(parallel)
    library(TropFishR)
    library(fishboot)
    # library(ks)
    library(fishmethods)
    
    opar <- par()   # save default graphical parameters (for plotting) 
    
    # VII.1 Load data --------------------------

    path = dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(path)
    
    # S stock
    age.length.rawVS3 <- read.csv("age-length-rawVS.csv")
    
    #   View(age.length.rawVS3) 
    summary(age.length.rawVS3)
    hist(age.length.rawVS3$length_cm)
   
    age.length.rawVS3$Sex <- as.factor(age.length.rawVS3$Sex)
    
    summary(age.length.rawVS3)
    length(age.length.rawVS3$age_yr) # N = 155 data, S stock
    
    LAAS <-age.length.rawVS3
    
    # N stock
    age.length.rawVN3 <- read.csv("age-length-rawVN.csv")
    #   View(age.length.rawVN3)     
    
    age.length.rawVN3$Sex <- as.factor(age.length.rawVN3$Sex)
    summary(age.length.rawVN3)
    
    length(age.length.rawVN3$age_yr) # N = 104 data, N stock
    
    LAAN <-  age.length.rawVN3
    
    
    # VII.2 first plots ----------------------------
    
    plot(LAAS$length_cm ~ LAAS$age_yr, main = "S stock")
    plot(LAAN$length_cm ~ LAAN$age_yr, main = "N stock")
    
    summary(LAAS)
    # Fork length:length_cm ,  L RANGE = 26.3 TO 62.4 cm Fork length
    # Southern stock: narrow size range !
    hist(LAAS$length_cm, xlim = c(0, 100))
    
    summary(LAAN) # smaller and larger fish are only found in the N stock dataset !
    # Fork length:length_cm ,  L RANGE =  11.60  to  89.90 cm Fork length
    hist(LAAN$length_cm, xlim = c(0, 100))
    
    
    
    
    # Scatter plots with colours by sex ----------------- 
    
    # S stock, Scatter plot --------------------
    plot(LAAS$length_cm ~ LAAS$age_yr, main = "S stock",
         pch = 19,
         col = factor(LAAS$Sex))
    
    # Legend
    legend("topleft",
           legend = levels(factor(LAAS$Sex)),
           pch = 19,
           col = factor(levels(factor(LAAS$Sex))))
    
    
    
    # N stock, Scatter plot --------------------
    plot(LAAN$length_cm ~ LAAN$age_yr, main = "N stock",
         pch = 19,
         col = factor(LAAN$Sex),
         xlim = c(0, max(LAAN$age_yr) ) ,
         ylim = c(0, max(LAAN$length_cm)) )
    
    # Legend
    legend("topleft",
           legend = levels(factor(LAAN$Sex)),
           pch = 19,
           col = factor(levels(factor(LAAN$Sex))))
    
    
    
    # VII.3 first curve fits  (TropFishR::growth_length_age()) --------------
    
    library(TropFishR)
    
   #  ??TropFishR
    
   # ?growth_length_age
    
    # S Stock ------------------
    
    datS <- list(age = LAAS$age_yr, length = LAAS$length_cm)
    
    # outputS <- growth_length_age(param = datS, method = "LSM",
    #                             Linf_init = 70 , CI = TRUE, age_plot=NULL)
    # ERROR, no fit possible with LSM method....
    
    # growth_length_age(datS, method = "GullandHolt") # looks horrible...
    # 
    #  summary(outputS$mod)
    
    
    growth_length_age(datS, method = "BertalanffyPlot", Linf_est = 83.89 )
    # Linf  = 83.89 , Linf estimate from LFA, with bs = 2 MA3, Linfmin30
    # MA3, Linfmin30
    # Linf              K              t_anchor             C          
    # Min.   : 61.14   Min.   :0.01743   Min.   :0.06839   Min.   :0.01911  
    # 1st Qu.: 75.43   1st Qu.:0.07327   1st Qu.:0.38110   1st Qu.:0.40199  
    # Median : 83.89 
    
    
    # result from growth_length_age (S Stock)
    # Linf fixed !
    # method = "BertalanffyPlot", Linf_est = 83.89 )
    # "LSM" method in "growth_length_age" does not work for S stock!!!
    # $K
    # [1] 0.05534864
    # 
    # $t0
    # [1] -6.534579
     
    # N Stock ----------------------
     
    # N stock works OK
     
    datN <- list(age = LAAN$age_yr, length = LAAN$length_cm)
     
    output.N <- growth_length_age(param = datN, method = "LSM",
                                  Linf_init = 70 , CI = TRUE, age_plot=NULL)
      summary(output.N$mod)
      # 
      # Parameters:
      #   Estimate Std. Error t value Pr(>|t|)    
      # Linf 85.57933    4.06351  21.060  < 2e-16 ***
      #   K     0.07981    0.01120   7.125 1.57e-10 ***
      #   t0   -3.07853    0.55945  -5.503 2.84e-07 ***
      
      # VII.4 Nice plot, N Stock -------------------
      
      # N stock, Scatter plot --------------------
   
      output.N <- growth_length_age(param = datN, method = "LSM",
                                    Linf_init = 83.9 , CI = TRUE, age_plot=NULL,
                                    )
      
      
         points(LAAN$length_cm ~ LAAN$age_yr,
           pch = 19,
           col = factor(LAAN$Sex),
           )
      
      # Legend
      legend("topleft",
             legend = levels(factor(LAAN$Sex)),
             pch = 19,
             col = factor(levels(factor(LAAN$Sex))))
      
      
      
    
    
    # VII.5 Define and apply function lenage_boot (fishboot) ------------
    
    datN <-  as.data.frame (datN)
      
    # Define function lenage_boot (fishboot) ------------------
      
    lenage_boot <- function(input.data, nboot = 200) {
      
      data.len.age <- input.data
      
      B = nboot ## number of bootstraps
      res = data.frame(Linf = numeric(B) , K = numeric(B), t0 = numeric(B)) ## vector to hold results
      n = length(data.len.age$length) # number of data pairs
      
      
      for(b in 1:B){
        
        
      #  tryCatch({   
          
          seed <- round(as.numeric(Sys.time())+runif(1, min = 0, max = 1e4)) # maximize stochasticity
          set.seed(seed)
          
          i = sample(x = 1:n, size = n, replace = TRUE) ## sample indices
          
          bootsam.dataxy <- data.len.age[i,] ## get data
          bootsam.dataxy <- list(age=  round(bootsam.dataxy$age,0) , length= bootsam.dataxy$length)
          
          
          output <- growth_length_age(param = bootsam.dataxy, method = "LSM",
                                      Linf_init = (max(bootsam.dataxy$length)), CI = FALSE, age_plot=NULL)
          
          output$Linf
          
          ## store results
          res$Linf[b]  <-   output$Linf 
          res$K[b]  <- output$K  
          res$t0[b]<- output$t0
          
       # }, error=function(e){})
        
      }
      
      ret <- list()
      
      ret$bootRaw <- res
      
      ret$seed <- seed
      
      class(ret) <- "lfqBoot"
      
      
      
      return(ret)
      
      
    }
     
          
    names(datN) <- c("age", "length"  )
   
      
    # apply function lenage_boot (fishboot) ---------
      
    res.LAA.N.OK <- lenage_boot(datN, 1000)
    
      
      
    summary (res.LAA.N.OK$bootRaw) 
    # Linf              K                 t0        
    # Min.   : 66.59   Min.   :0.03966   Min.   :-6.051  
    # 1st Qu.: 81.69   1st Qu.:0.07081   1st Qu.:-3.539  
    # Median : 85.81   Median :0.08051   Median :-3.034  
    # Mean   : 85.73   Mean   :0.08258   Mean   :-3.073  
    # 3rd Qu.: 89.58   3rd Qu.:0.09141   3rd Qu.:-2.578  
    # Max.   :110.09   Max.   :0.15636   Max.   :-1.051   
    
    # save bootstrap posteriors results (1000 bootstrap posteriors), N stock -----------
    
    write.csv (res.LAA.N.OK$bootRaw,
                      file = "res_LAA_N_OKRS_1000.csv")
    
   
     
    # VII.6  curve swarms  (fishboot) and other FINAL plots for paper ------------
    
    # VII.6.1. simple scatter plot with median VBGF curve 
    # Linf = 85.81 cm, just the curve and the points
    # Medians from lenage_boot (N Stock)
    # Median Linf : 85.81, # median K = 0.08051, Median t_zero:-3.034  
  
    max(datN$age)
    #  38, plot to age = 40 years! 
    max(datN$length)
    #  89.9 cm, plot to Fork length = 90 cm! 
    
    library(TropFishR)
    
    t <- seq(-3,40,0.1)
    Lt <- VBGF(list(Linf=85.81, K=0.08051, t0= -3.034  ),t=t)
    plot(t, Lt, t="l", ylim = c(0,90), xlim = c(-5,40), col = "red")
    
    #    insert median optimum curve (median of posteriors), red line
    
    lines(t, Lt, t="l", col = "red")
    
    points(LAAN$length_cm ~ LAAN$age_yr,
           pch = 19,
           col = factor(LAAN$Sex),
    )
    
    # Legend
    legend("topleft",
           legend = levels(factor(LAAN$Sex)),
           pch = 19,
           col = factor(levels(factor(LAAN$Sex))))
    
    # VII.6.1. curve swarm (fishboot), with median VBGF curve and points
    
    # VII.6.1.a LAA N stock figure, for PAPER ------------------------------------
    
    
    # curve swarm (fishboot) --------------
    # with median VBGF curve  (red line) ---------
    # with color points coded by sex
    
    library(fishboot)
    
    
    opar6 <- par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(0,0,0,0), mgp = c(2,0.5,0), 
                 tcl = -0.25,
                 cex = 1)
    
    par.current <- par()
    
    # rename t0 -> t_anchor, needed for plotting function
    
    names(res.LAA.N.OK$bootRaw ) <- c ("Linf", "K", "t_anchor")  
    
    CIinfo <- vbgfCI_time(
      res = res.LAA.N.OK,
      agemax = 40, CI = 95,
      perm.col = adjustcolor("grey50",0.2)
    )
   
    lines(t, Lt, t="l", col = "red")
    
    points(LAAN$length_cm ~ LAAN$age_yr,
           pch = 19,
           col = factor(LAAN$Sex),
    )
    
    legend("topleft",
           legend = levels(factor(LAAN$Sex)),
           pch = 19,
           col = factor(levels(factor(LAAN$Sex))))
    
    
    # # VII.6.1.b LAA S stock figure, for PAPER ------------------------------------
    
    # load data for plot (1000 posteriors) -------------------
    
    path = dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(path)
    
    #   read data -------------------
    
    oto.S.all.1000results_b <- read.csv("oto-S-all-1000results.csv")
    
    #  View(oto.S.all.1000results_b)
   
    library(TropFishR) 
    library(fishboot)
    
    datN <- read.csv("age-length-rawVN.csv")
    
    res.LAA.N <- lenage_boot(datN, 1000)
    
    res.LAA.S.OK <- res.LAA.N

  
    res.LAA.S.OK$bootRaw <- oto.S.all.1000results_b
    
    # Comment: Ctrl + Shift + C
    
  
    head(res.LAA.S.OK$bootRaw)
    
    attach(res.LAA.S.OK$bootRaw)

    # insert a Phi' column    
    res.LAA.S.OK$bootRaw$Phi_L <- log10(res.LAA.S.OK$bootRaw$K) + (2 * log10(res.LAA.S.OK$bootRaw$Linf))

    
    head(res.LAA.S.OK$bootRaw)# OK
    
    
    # Median values from posteriors (S Stock):
    summary(res.LAA.S.OK$bootRaw)
    
    # Median Linf:  77.29 ,  Median :0.066345,   Median t0 :-5.712,  Median Phi': 2.617  
   
    # Linf              K                  t0              Phi_L      
    # Min.   : 57.12   Min.   :0.002288   Min.   :-15.346   Min.   :2.511  
    # 1st Qu.: 69.85   1st Qu.:0.043391   1st Qu.: -7.699   1st Qu.:2.590  
    # Median : 77.29   Median :0.066345   Median : -5.712   Median :2.617  
    # Mean   : 99.90   Mean   :0.065436   Mean   : -6.289   Mean   :2.633  
    # 3rd Qu.: 93.57   3rd Qu.:0.087917   3rd Qu.: -4.417   3rd Qu.:2.648  
    # Max.   :915.22   Max.   :0.170270   Max.   : -1.719   Max.   :3.304  
    
    
    
    
    ### ### 
    # load raw data (for plotting color points) --------------------
    # S stock
    age.length.rawVS3 <- read.csv("age-length-rawVS.csv")
   

    head(age.length.rawVS3)
    
    # View(age.length.rawVS3) 
    summary(age.length.rawVS3)
    hist(age.length.rawVS3$length_cm)
    
    
    age.length.rawVS3$Sex <- as.factor(age.length.rawVS3$Sex)
    
    summary(age.length.rawVS3)
    length(age.length.rawVS3$age_yr) # N = 155 data, S stock
    
    LAAS <-age.length.rawVS3
    
    LAAS$age_yr <- as.numeric( LAAS$age_yr)
    
    
    plot(LAAS$length_cm ~ LAAS$age_yr, main = "S stock")
     
    summary(LAAS)
    # Fork length:length_cm ,  L RANGE = 26.3 TO 62.4 cm Fork length
    # Southern stock: narrow size range !
    
    hist(LAAS$length_cm, xlim = c(0, 100))
    
    # Scatter plots with colors by sex ----------------- 
  
    # N stock, Scatter plot --------------------
    plot(datN$length_cm ~ datN$age_yr, main = "N & S stock, S: orange",
         pch = 19,
         col = factor(datN$Sex),
         ylim = c(0,90), xlim = c(0, 40) )
    
    # Legend
    legend("topleft",
           legend = levels(factor(datN$Sex)),
           pch = 19,
           col = factor(levels(factor(datN$Sex))))
    
    
      
    # S stock, Scatter plot --------------------
    points(LAAS$length_cm ~ LAAS$age_yr, main = "S stock",
         pch = 6, # pch = 19,
         col = "darkorange", #  col = factor(LAAS$Sex),
         ylim = c(0,90), xlim = c(0, 40) )
    
    # Legend
    legend("bottomright",
           legend = levels(factor(LAAS$Sex)),
            pch = 6,  #pch = 19,
             col = "darkorange") #   col = factor(levels(factor(LAAS$Sex))))
    
    
    
    # curve swarm (fishboot) FOR PAPER --------------
    # with median VBGF curve  (red line) ---------
    # with color points coded by sex
    
    library(fishboot)
    
    
    opar6 <- par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(0,0,0,0), mgp = c(2,0.5,0), 
                 tcl = -0.25,
                 cex = 1)
    
    par.current <- par()
    
    # rename t0 -> t_anchor, needed for plotting function
    
    names(res.LAA.S.OK$bootRaw ) <- c ("Linf", "K"  ,  "t_anchor", "Phi_L")  
    
    head(res.LAA.N.OK$bootRaw) 
    
    summary(res.LAA.S.OK$bootRaw)
    
    class(res.LAA.S.OK$bootRaw)
    class(res.LAA.N.OK)
    
    CIinfo <- vbgfCI_time(
      res = res.LAA.S.OK,
      agemax = 40, CI = 95,
      perm.col = adjustcolor("grey50",0.2)
    )
    
     
  
    
    # S stock, Scatter plot --------------------
    points(LAAS$length_cm ~ LAAS$age_yr,#  main = "S stock",
         pch = 19,
         col = factor(LAAS$Sex))
    
    # Legend
    legend("topleft",
           legend = levels(factor(LAAS$Sex)),
           pch = 19,
           col = factor(levels(factor(LAAS$Sex))))
    
  
    # Median values from posteriors:
    # Median Linf:  77.29 ,  Median K :0.066345,   Median t0 :-5.712,  
    # Median Phi': 2.617  
    
    library(TropFishR)

    t <- seq(-3,40,0.1)
    Lt <- VBGF(list(Linf= 77.29, K=0.066345, t0= -5.712 ),t=t)
    #plot(t, Lt, t="l", ylim = c(0,90), xlim = c(-5,40), col = "red")
    
    #    insert median optimum curve (median of posteriors), red line
    
    lines(t, Lt, t="l", col = "red")
   
    ########
    ######
    ####
    ###
    
    # VIII. Chapter VIII: simple LAA (otolith reading data) plots for paper,  --------------------------------------
    #      without "curve swarms", ages by sex and "optimum curve" --------------------
    #       Northern  vs Southern stock 
    #       Direct testing of medion length,  Northern  vs Southern stock
    #       Permutation test Northern  vs Southern stock, median length at age
    #       Simple boxplots
    #       Simple  plot (symbols and cols by region)
    
    
    # clean memory
    gc()  
    gc(reset=T)
    rm(list = ls()) 
    
    # install.packages("fishmethods")
    # install.packages("devtools")
    # library(devtools)
    # install_github("rschwamborn/fishboot")
    # install.packages("devtools")
    
    
    # required packages =====
    library(parallel)
    library(TropFishR)
    library(fishboot)
    # library(ks)
    library(fishmethods)
    
    opar <- par()   # save default graphical parameters (for plotting) 
    
    
    # VIII.1 Load data --------------------------
    
    path = dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(path)
    
    # S stock
    age.length.rawVS3 <- read.csv("age-length-rawVS.csv")
    
    #   View(age.length.rawVS3) 
    summary(age.length.rawVS3)
    hist(age.length.rawVS3$length_cm)

    age.length.rawVS3$Sex <- as.factor(age.length.rawVS3$Sex)
    
    summary(age.length.rawVS3)
    length(age.length.rawVS3$age_yr) # N = 155 data, S stock
    
    LAAS <-age.length.rawVS3
    
    # N stock
    age.length.rawVN3 <- read.csv("age-length-rawVN.csv")
    #   View(age.length.rawVN3)     
    
    age.length.rawVN3$Sex <- as.factor(age.length.rawVN3$Sex)
    summary(age.length.rawVN3)
    
    length(age.length.rawVN3$age_yr) # N = 104 data, N stock
    
    LAAN <-  age.length.rawVN3
    
    
    # VIII.2 simple plots ----------------------------
    
    plot(LAAS$length_cm ~ LAAS$age_yr, main = "S stock")
    plot(LAAN$length_cm ~ LAAN$age_yr, main = "N stock")
    
    summary(LAAS)
    # Fork length:length_cm ,  L RANGE = 26.3 TO 62.4 cm Fork length
    # Southern stock: narrow size range !
    hist(LAAS$length_cm, xlim = c(0, 100))
    
    summary(LAAN) # smaller and larger fish are only found in the N stock dataset !
    # Fork length:length_cm ,  L RANGE =  11.60  to  89.90 cm Fork length
    hist(LAAN$length_cm, xlim = c(0, 100))
    
    
    # Define one age  classe with at least 30 fish  ------------
    # "medium age" AGE CLASS: 4 to 12 years
    
    hist(LAAS$age_yr, xlim = c(0, 40))
    hist(LAAN$age_yr, xlim = c(0, 40))

  # S stock   
   length (LAAS$age_yr)# All fish: 155 ind. (S stock)
   medium.aged.S <- subset (LAAS, age_yr<13 &  age_yr > 3) # Medium-aged fish: 155
   length(medium.aged.S$age_yr)
   
   plot(LAAS$length_cm ~ LAAS$age_yr, main = "S stock",
        pch = 19,
        col = factor(LAAS$Sex), 
        xlim = c(0, 40 ) , 
        ylim = c(0, 90) )
   abline (v = c(3, 12))
   
   
   # N stock   
   length (LAAN$age_yr)# All fish: 104 ind. (N stock)
   medium.aged.N <- subset (LAAN, age_yr<13 &  age_yr > 3) # Medium-aged fish: 155
   length(medium.aged.N$age_yr)
   
   summary(medium.aged.N) # age: 4 to 12 years
  
    plot(LAAN$length_cm ~ LAAN$age_yr, main = "N stock",
        pch = 19,
        col = factor(LAAN$Sex), 
        xlim = c(0, 40 ) , 
        ylim = c(0, 90) )
    abline (v = c(3.5, 12.5))
    
    
   # plots for PAPER ------------------------------
   
   # S stock (orange)  ----------
   length (LAAS$age_yr)# All fish: 155 ind. (S stock)
   medium.aged.S <- subset (LAAS, age_yr<13 &  age_yr > 3) # Medium-aged fish: 155
   length(medium.aged.S$age_yr)
   
   plot(LAAS$length_cm ~ LAAS$age_yr,
        pch = 6,
        col = "darkorange", 
        xlim = c(0, 40 ) , 
        ylim = c(0, 90) ,
        xlab = "Age (y)",
        ylab = "Fork length (cm)")
   
   
   abline (v = c(3.5, 12.5))
   
   # N stock (blue)  ---------
   
   library (scales)
   
 
   length (LAAN$age_yr)# All fish: 104 ind. (N stock)
   medium.aged.N <- subset (LAAN, age_yr<13 &  age_yr > 3) # Medium-aged fish: 155
   length(medium.aged.N$age_yr)
   
   summary(medium.aged.N) # age: 4 to 12 years
   
   points(LAAN$length_cm ~ LAAN$age_yr,
        pch = 19,
        col =  alpha("navyblue", 0.3) 
       )
 
   df.S.N.test <- data.frame ( Stock = c( rep("N", length(medium.aged.N$length_cm)), rep("S", length(medium.aged.S$length_cm))  )  ,
                               length = c(medium.aged.N$length_cm , medium.aged.S$length_cm) )
   
   df.S.N.test$Stock <- as.factor(df.S.N.test$Stock)
   
   
   # Legend
   legend("bottomright",
          legend = levels(factor(df.S.N.test$Stock)),
          pch = 19,
          col = c(alpha("navyblue", 0.3)  , "darkorange"))
   
   
   
   
   
   
  # boxplot N vs S (medium-aged fish only) -----------
   
   df.S.N.test <- data.frame ( Stock = c( rep("N", length(medium.aged.N$length_cm)), rep("S", length(medium.aged.S$length_cm))  )  ,
                               length = c(medium.aged.N$length_cm , medium.aged.S$length_cm) )
   
   df.S.N.test$Stock <- as.factor(df.S.N.test$Stock)
   
   # View(df.S.N.test)
   
   boxplot(df.S.N.test$length ~df.S.N.test$Stock,
           ylab = "Fork length (cm)", xlab = "Stock",
           col = c(alpha("navyblue", 0.3)  , "darkorange"))
   
   
   median (medium.aged.S$length_cm) # 44.9
   median (medium.aged.N$length_cm)# 44.5 
   
   
   # VIII.3. PERMUTATION TESTS  --------------------------------------

   # N vs S stock, median length at age , 
   # for "medium-aged" (4 to 12 y) fish
   # (mostly young adults, mostly males and hermaphr.)
   
   df.S.N.test <- data.frame ( Stock = c( rep("N", length(medium.aged.N$length_cm)), rep("S", length(medium.aged.S$length_cm))  )  ,
                                length = c(medium.aged.N$length_cm , medium.aged.S$length_cm) )
    
   df.S.N.test$Stock <- as.factor(df.S.N.test$Stock)
     
  # View(df.S.N.test)
   
   boxplot(df.S.N.test$length ~df.S.N.test$Stock)
   
   # three different approaches are tested: package "coin", 
   # package "vegan", 
   # and  using the basic sample() function 
   
   # coin::independence_test( ) --------------------
   library(coin)

   independence_test(df.S.N.test$length ~ df.S.N.test$Stock)   
# p-value = 0.4516, n.sign.
    
# vegan::adonis2() ------------------
library(vegan)  

#?adonis2

adonis2(df.S.N.test$length ~ df.S.N.test$Stock, method = "euclidean" ) 
# p = 0.434, n.sign.

# using the basic sample () function, two scripts are tested ------------------
# 1. Script taken from  https://stats.stackexchange.com/questions/176691/running-a-permutation-test-with-different-sample-sizes-in-r

# 1.a diff. in means
t_stat=mean(medium.aged.N$length_cm)-mean(medium.aged.S$length_cm) 
combine = c("N", "S") 
difference=rep(0,924) 
for (i in c(1:924)) { 
  y=sample(combine) 
  t_stat1=mean(y[1:6])-mean(y[7:12])
  difference[i]=t_stat1 } 
p_value=sum(abs(difference)>abs(t_stat))/924 
p_value # 0.69 no signif. for difference in means

# 1.b diff. in medians
t_stat=median(medium.aged.N$length_cm)-median(medium.aged.S$length_cm) 
combine=c(N,S) 
difference=rep(0,924) 
for (i in c(1:924)) { 
  y=sample(combine) 
  t_stat1=median(y[1:6])-median(y[7:12])
  difference[i]=t_stat1 } 
p_value=sum(abs(difference)>abs(t_stat))/924 
p_value # p = 0.88095 no signif. for difference in means ---------


# 2. Script taken from   https://rpubs.com/rschwamborn/288625

# https://rpubs.com/rschwamborn/288625  
# Permanova and Bootstrap - simple tutorial
# Ralf Schwamborn
# 2 de julho de 2017
# PERMANOVA and Bootstrap
# Cheat sheet - Simple tutorial
# How to do a Permutation test (PERMANOVA) and Bootstrap
# (R. Schwamborn, 2017)
# I. PERMANOVA

# diff in means ---------------- 
am.A <- medium.aged.N$length_cm

am.B <- medium.aged.S$length_cm

# am.A <- 1:5 # sample A
# 
# am.B <- 11:15 # sample B

# Difference between the means of samples A and B
diff.AB.means.orig <- mean(am.B)-mean(am.A)
diff.AB.means.orig   # -1.113362


diff.AB.medians <- median(am.B)-median(am.A)
diff.AB.medians  #0.4

#-1.1 cm diff in means, raw data (medium aged fish) 
# 0.4 cm diff in medians, raw data (medium aged  fish) 

## I.0. Basic Question: Is the observed difference significant?

# standard test: Student's t-test 
#Student's t-test  ---------------
t.test(am.A, am.B) # p-value = 0.5687, not signif.

# Mann-whitney U test 
stats::wilcox.test(am.A, am.B) # p-value = 0.8675, not signif.




Data.3 <- df.S.N.test


#View(Data.3)

names(Data.3) # "Stock" "length"   

names(Data.3) <- c( "grupo", "values" )   
names(Data.3) # "grupo" "values"  

summary(Data.3)
# grupo     
# N: 31     
# S:147   


## [1] "values" "grupo"
t.test(Data.3$values ~Data.3$grupo) # p-value = 0.5687, not signif.
       
coin::wilcox_test (Data.3$values ~Data.3$grupo) # p-value = 0.866, not signif.



                     
## I.2. Take a random sample with values all mixed up,
## but factors ("A", "B") unchanged
sample (Data.3$values,  replace = FALSE)

## I.3. Take such a mixed up random sample, 
## then create a new data frame with this sample
Data.sam1 <-Data.3 
Data.sam1$values <- sample (Data.3$values, replace = FALSE)
##View(Data.sam1)
sam.A <- Data.sam1$values[1:31]
sam.B <- Data.sam1$values[32:138]
t.test(sam.A, sam.B) # p-value = n.s.
## 
##  Welch Two Sample t-test
## 
## data:  sam.A and sam.B
## t = 0.67135, df = 7.6582, p-value = 0.5217
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##-5.908252 10.708252
## sample estimates:
## mean of x mean of y 
##       9.2       6.8
## I.4. calculate the difference between "A" and "B" in the "mixed up" random sample 
diff.AB <- mean(sam.B)-mean(sam.A)
diff.AB   # difference between "A" and "B" in the a "mixed up" random sample 
## [1] 2.16 cm diff in means in the a "mixed up" random sample 
diff.AB.medians <- median(sam.B)-median(sam.A)
diff.AB.medians   # difference between "A" and "B" in the a "mixed up" random sample 
# 1.7 cm diff in means in the a "mixed up" random sample 

## I.5. Create a function that makes all this. 
## It takes a "mixed up" random sample  
## and then calculates the difference between "A" and "B"
## new function "sam_diff.fun1()"

# As Function

sam_diff.fun1 <- function (DF) {
  
  #take a sample
  Data.sam1 <-DF
  colnames(Data.sam1)[1] <- "values"
  Data.sam1$values <- sample (Data.3$values,(length(Data.3$values)),  
                              replace = FALSE)
  ##View(Data.sam1)
  sam.A <- Data.sam1$values[1:31]
  sam.B <- Data.sam1$values[32:138]
  
  diff.AB <- mean(sam.A)-mean(sam.B)
  ; diff.AB
  
  
} # END OF FUNCTION

sam_diff.fun1(Data.3)
## [1] -1.29

## I.6. Apply this function several times (10000 times)
# Each time, you get a diferent sample 
# and a different value for the diffence between means in the sample 
# By this, you will get 10000 values of diffence between means 
# the distrib. of these 10000 values is the posterior distib., 
# an approximation of the null distrib.. 
# Analyse: Is your original observed mean within the null distribution?

# Loop to obtain the null distribution

Nperm <- 10000 # Number of permutations
diff.vec <-  rep(0,Nperm)

for (i in 1:Nperm) {  
  
  diff.sam <-  sam_diff.fun1(Data.3)
  
  diff.vec[i] <- diff.sam
  
} # END OF LOOP

## I.7. Analyse the distribution obtained by permutation
# i.e., analyse the shape of the distribution obtained by repeated mixing up and sampling. 
# and determine where your original observed mean is located within the null distribution.

# Analyze Loop results

summary(diff.vec)
## Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## -5.883120 -1.016493 -0.005324  0.016234  1.058493  5.424508 
hist.graph1 <- hist(diff.vec, breaks = seq(-20, 20, length.out = 201))


hist.graph1$counts

hist.graph1 <- hist(diff.vec, breaks = seq(-20, 20, length.out = 201),
                    prob = TRUE,
                    main = "Null Distribution of the diff. in means")
abline(v = diff.AB.means.orig, col = "green")
xydens<- density(diff.vec, from = -20, to= 20,n = 201)
xydens$y <- xydens$y*2
lines(xydens, col = grey(0.3), lty = 2)              
text(-11.4, 0.3, "obs.", col = "green", cex = 0.8)

## now calculate "p" and 95% Conf.Interval
# calculate95% Conf.Interval
quantile(diff.vec, c(0.025, 0.975)) # 95% Conf.Interval:-6.8 to 6.8

CI_lower <- quantile(diff.vec, c(0.025)) # 95% Conf.Interval:-6.8 to 6.8
CI_upper <- quantile(diff.vec, c(0.975)) # 95% Conf.Interval:-6.8 to 6.8


##  2.5%        97.5% 
##  -2.971161  3.023621 
abline(v = CI_lower,col = "grey"); abline(v = CI_upper,col = "grey");
text(3, 0.32, "95% Conf. int.", col = "grey", cex = 0.8)

#  calculate "p"  based on counts (counting number of samples)
df.hist1 <- data.frame(hist.graph1$density,round(hist.graph1$mids,2), 
                       hist.graph1$counts) 
#View(df.hist1)
names(df.hist1)<-  c("density", "mids", "counts")
#View(df.hist1)

df.hist1.10 <- subset(df.hist1, df.hist1$mids < diff.AB.means.orig) 
#View(df.hist1.10)

sum(df.hist1.10$counts) # 2137 permutation runs out of 10,000 gave "diff = -10"
## [1] 2137
p_val <- 2 *  sum(df.hist1.10$counts)/Nperm # p = 0.8511
p_val ## [1] 0.4274, not signif.
#  calculate "p"  based on the kernel density



### TEST for PAPER ----------------
# test for difference  in medians -------------

am.A <- medium.aged.N$length_cm

am.B <- medium.aged.S$length_cm

# am.A <- 1:5 # sample A
# 
# am.B <- 11:15 # sample B

length(am.A) # 31 medium-aged fish (N stock, 4 to 12 Y) 
length(am.B) # 147 medium-aged fish (S stock, 4 to 12 Y) 
+31+147    # total: 178 medium-aged fish

# Difference between the means of samples A and B
diff.AB.medians.orig <- median(am.B)-median(am.A)
diff.AB.medians.orig   # 0.4 cm difference in medians of original data 
#                         medium aged fish (31, 147) 

diff.AB.medians <- median(am.B)-median(am.A)
diff.AB.medians  

# 0.4 cm diff in medians, raw data (medium aged  fish) 

## I.0. Basic Question: Is the observed difference significant?

# standard test: Student's t-test 
#Student's t-test  ---------------
t.test(am.A, am.B) # p-value = 0.5687, not signif.

# Mann-whitney U test 
stats::wilcox.test(am.A, am.B) # p-value = 0.8675, not signif.


Data.3 <- df.S.N.test


#View(Data.3)

names(Data.3) # "Stock" "Length"   

names(Data.3) <- c( "grupo", "values" )   
names(Data.3) # "grupo" "values"  

summary(Data.3)
# grupo     
# N: 31     
# S:147   


## [1] "values" "grupo"
t.test(Data.3$values ~Data.3$grupo) # p-value = 0.5687, not signif.

coin::wilcox_test (Data.3$values ~Data.3$grupo) # p-value = 0.866, not signif.




## I.5. Create a function that makes all this. 
## It takes a "mixed up" random sample  
## and then calculates the difference between "A" and "B"
## new function "sam_diff.fun1()"

# As Function

sam_diff.fun1 <- function (DF) {
  
  #take a sample
  Data.sam1 <-DF
  colnames(Data.sam1)[1] <- "values"
  Data.sam1$values <- sample (Data.3$values,(length(Data.3$values)),  
                              replace = FALSE)
  ##View(Data.sam1)
  sam.A <- Data.sam1$values[1:31]
  sam.B <- Data.sam1$values[32:138]
  
  diff.AB <- median(sam.A)-median(sam.B)
  ; diff.AB
  
  
} # END OF FUNCTION

sam_diff.fun1(Data.3)
## [1] 0.6
## I.6. Apply this function several times (10000 times)
# Each time, you get a diferent sample 
# and a different value for the diffence between medians in the sample 
# By this, you will get 10000 values of diffence between medians 
# the distrib. of these 10000 values is the posterior distib., 
# an approximation of the null distrib.. 
# Analyse: Is your original observed median within the null distribution?

# Loop to obtain the null distribution

Nperm <- 10000 # Number of permutations
diff.vec <-  rep(0,Nperm)

for (i in 1:Nperm) {  
  
  diff.sam <-  sam_diff.fun1(Data.3)
  
  diff.vec[i] <- diff.sam
  
} # END OF LOOP

## I.7. Analyse the distribution obtained by permutation
# i.e., analyse the shape of the  distribution obtained by repeated mixing up and sampling. 
# and determine where your original observed median is located within the null distribution.

# Analyze Loop results

summary(diff.vec)
##      Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##  -7.20000 -1.50000  0.00000  0.02234  1.40000  9.20000 
hist.graph1 <- hist(diff.vec, breaks = seq(-20, 20, length.out = 201))


hist.graph1$counts

hist.graph1 <- hist(diff.vec, breaks = seq(-20, 20, length.out = 201),
                    prob = TRUE,
                    main = "Null Distribution of the diff. in medians")
abline(v = diff.AB.medians.orig, col = "green")
xydens<- density(diff.vec, from = -20, to= 20,n = 201)
xydens$y <- xydens$y*2
lines(xydens, col = grey(0.3), lty = 2)              
text(-11.4, 0.3, "obs.", col = "green", cex = 0.8)

## now calculate "p" and 95% Conf.Interval
# calculate95% Conf.Interval
quantile(diff.vec, c(0.025, 0.975)) # 95% Conf.Interval: -3.9 to 4.7

CI_lower <- quantile(diff.vec, c(0.025)) 
CI_upper <- quantile(diff.vec, c(0.975)) 


##  2.5% 97.5% 
##  -3.9   4.7
abline(v = CI_lower,col = "grey"); abline(v = CI_upper,col = "grey");
text(3, 0.32, "95% Conf. int.", col = "grey", cex = 0.8)

#  calculate "p"  based on counts (counting number of samples)
df.hist1 <- data.frame(hist.graph1$density,round(hist.graph1$mids,2), 
                       hist.graph1$counts) 
#View(df.hist1)
names(df.hist1)<-  c("density", "mids", "counts")
#View(df.hist1)

df.hist1.10 <- subset(df.hist1, df.hist1$mids < diff.AB.medians.orig) 
#View(df.hist1.10)

sum(df.hist1.10$counts) # 36 permutation runs out of 10,000 gave "diff = -10"
## [1] 5585
p_val <- sum(df.hist1.10$counts)/Nperm 
p_val ## [1]  p = 0.5585, not signif.




    
    # Scatter plots with colours by sex ----------------- 
    
    # S stock, Scatter plot --------------------
    plot(LAAS$length_cm ~ LAAS$age_yr, main = "S stock",
         pch = 19,
         col = factor(LAAS$Sex), 
         xlim = c(0, 40 ) , 
         ylim = c(0, 90) )
    abline (h = c(30, 50))
    
    # Legend
    legend("topleft",
           legend = levels(factor(LAAS$Sex)),
           pch = 19,
           col = factor(levels(factor(LAAS$Sex))))
    
    
    
    # N stock, Scatter plot --------------------
    plot(LAAN$length_cm ~ LAAN$age_yr, main = "N stock",
         pch = 19,
         col = factor(LAAN$Sex),
          xlim = c(0, 40 ) , 
          ylim = c(0, 90) )
    abline (h = c(30, 50))
    
         max(LAAN$age_yr) # 38 yrs
         max (LAAN$length_cm) # 89.9 cm
    # Legend
    legend("topleft",
           legend = levels(factor(LAAN$Sex)),
           pch = 19,
           col = factor(levels(factor(LAAN$Sex))))
    
    
      
    
    
    
    
    # IX. Chapter IX. Recruitment analysis (numbers of small indiv.) -------------
    #         Numbers of small fish, and description of "first peak"
    
    
    ## IX. RECRUITMENT #######
    # LFD plots (S tock) --------------------------
    
    
    # Data: 
    # Fish lengths (fork length) from 2004-08-25 to 2009-12-30 
    # (approx. five and a half years) 
    # Length range: 12 cm to 79 cm
    # N = 1399 indiv.
    
    
    #Initialization ----------------------------------------------------------
    
    #  required packages =====
    library(parallel)
    library(TropFishR)
    library(fishboot)
    library(ks)
    
    # clean memory
    gc()  
    gc(reset=T)
    rm(list = ls())
    
    
    #  Load data =====
    
    lfq2 <- read.csv("stbrasLF-trunc.csv", header = T)
    
    # View(lfq2)
    
    barplot(lfq2$Length_CM)
    barplot(lfq2$Frequency)
    
    summary(lfq2)
    
    # catch curve (All years)
    plot(lfq2$Frequency ~lfq2$Length_CM)
    plot(lfq2$Frequency ~lfq2$Length_CM)
    
    plot(log10(lfq2$Frequency) ~lfq2$Length_CM)
    
    lfq2
    names(lfq2)
    lfq2$DATE <- as.Date(lfq2$DATE, format = "%d/%m/%Y") #Date format needs to be input like that
    lfq2
    
    summary(lfq2)# 12 to 79 cm
    length(lfq2$Length_CM) # 1399 fish cm
    
    
    lfq2new <- lfqCreate(data = lfq2, Lname = "Length_CM", Dname = "DATE", Fname = "Frequency")
   
     plot(lfq2new, Fname = "catch")
  
     # plot catch curve (All years) -----------
     
     Ncatch <-  rowSums(lfq2new$catch) 
     Ncatch.NA  <-  replace (Ncatch , 0, NA) 
    
    plot(Ncatch ~lfq2new$midLengths)
    plot(Ncatch.NA ~lfq2new$midLengths , main = "no zeros shown" )
   
    # RECRUITMENT -----------
    # N OF SMALL-SIZED FISH
    # make a plot with date (x axis) and N of small indiv. --------
    
    # Date vector (X) ------------
    
    Date <-  lfq2new$dates  
    
    lfq2new$midLengths# # 81 length classes, 
    #                    from 0.5 to 80.5 cm

    length(lfq2new$midLengths)    
    
    lfq2new$catch
    
    # sm21, 22, 23, 24, 25 ... : small ones <23, 24, 25, cm -------------------
   
    # sm21 : small ones <21 cm --------------
    
    lfq2sm21 <- subset ( lfq2 , Length_CM < 21)# 151 fish
    
     
    # sm22 : small ones <22 cm --------------
    
    lfq2sm22 <- subset ( lfq2 , Length_CM < 22)# 189 fish
    
    
     # sm23 : small ones <23 --------------
    
    lfq2sm23 <- subset ( lfq2 , Length_CM < 23)# 209 fish
    
    
    # sm24 : small ones <24 --------------
    
    lfq2sm24 <- subset ( lfq2 , Length_CM < 24)# 238 fish
    
    
    # sm25 : small ones <25 --------------
    
    lfq2sm25 <- subset ( lfq2 , Length_CM < 25)# 267 fish
   
    # sm26 : small ones <26 ---------------
    
    lfq2sm26 <- subset ( lfq2 , Length_CM < 26)# 296 fish
    
    # sm27 : small ones <27 ---------------
    
    lfq2sm27 <- subset ( lfq2 , Length_CM < 27)# 325 fish
    
    # sm28 : small ones <28 ---------------
    
    lfq2sm28 <- subset ( lfq2 , Length_CM < 28)# 354 fish
    
    # sm29 : small ones <29 ---------------
    
    lfq2sm29 <- subset ( lfq2 , Length_CM < 29)# 383 fish
    
    # FOR PAPER Catch curve --------------- 
   
    lfq2new
    
     barplot(   Ncatch.NA ~lfq2new$midLengths , 
               main = "catch curve, S stock, N = 1399",
               col = "lightblue")
    
    
    lfq2new <- lfqCreate(data = lfq2, Lname = "Length_CM", Dname = "DATE", Fname = "Frequency")
    
    plot(lfq2new, Fname = "catch")
    
    # plot catch curve (All years) -----------
    
    Ncatch <-  rowSums(lfq2new$catch) 
    Ncatch.NA  <-  replace (Ncatch , 0, NA) 
    
    plot( Ncatch ~lfq2new$midLengths  )
    plot( Ncatch.NA ~lfq2new$midLengths , main = "no zeros shown" )
    
    
    
    
    
    
    # FOR PAPER Catch curve --------------- 
    barplot(Ncatch.NA ~lfq2new$midLengths , 
               main = "catch curve, S stock, N = 1399",
               col = "lightblue")
    
     
    hist(lfq2$Length_CM)
    summary(lfq2$Length_CM)
    summary(lfq2$DATE)
    
     
     # Catch curve by month (as in TropfishR, but manually month by month)
        
     
     Ncatch1 <-  lfq2new$catch [,1]
     Ncatch.NA1  <-  replace(Ncatch1 , Ncatch1 == 0  , NA) 
     Ncatch.NA1
     
     
     
     plot(Ncatch1 ~lfq2new$midLengths  )
    
     plot(log10(Ncatch.NA1) ~lfq2new$midLengths ,
             main = paste(lfq2new$dates[1], " month 1, log10, no zeros shown") )
     abline(v = c(23, 25, 27, 29))
     text ("23, 25, 27, 29 cm", x = 15, y =1.25 ) 
     
     lfq2new$midLengths[25]
     lfq2new$midLengths[26]
     lfq2new$midLengths[27]
     lfq2new$midLengths[28]
     lfq2new$midLengths[29]
     lfq2new$midLengths[30]
     
     
     plot(   (Ncatch.NA1) ~lfq2new$midLengths ,
             main = paste(lfq2new$dates[1], " month 1, no zeros shown") )
     abline(v = c(23, 25, 27, 29))
     text ("23, 25, 27, 29 cm", x = 15, y =25 ) 
     text (paste ( "N = ",(sum(Ncatch1) )  ), x = 65, y =25 ) 
     text (paste ( "N (<29cm) = ",(sum(Ncatch1[1:29]  ) )  ), x = 65, y =21 ) 
     text (paste ( "N (<27cm) = ",(sum(Ncatch1[1:27] ) )  ), x = 65, y =17 ) 
     text (paste ( "N (<25cm) = ",(sum(Ncatch1[1:25] ) )  ), x = 65, y =13 ) 
     text (paste ( "N (<23cm) = ",(sum(Ncatch1[1:23] ) )  ), x = 65, y =9 ) 
     
     
     # MONTH-By-MONTH -------
     # (LOOP)
     
     # month 1
     
     Month = 1
     i = Month
     
     Ncatch.i  <-  lfq2new$catch [,i]
     Ncatch.NAi  <-  replace(Ncatch.i , Ncatch.i == 0  , NA) 
     Ncatch.NAi
     
     plot(   (Ncatch.NAi) ~lfq2new$midLengths ,
             main = paste(lfq2new$dates[Month], "no zeros shown") )
     abline(v = c(25, 27, 29, 31))
     text ("25, 27, 29, 31 cm", x = 15, y =max(Ncatch.i)*0.9) 
     text (paste ( "N = ",(sum(Ncatch.i) )  ), x = 65,  y =max(Ncatch.i)*0.9 ) 
     text (paste ( "N (<31cm) = ",(sum(Ncatch.i[1:31]  ) )  ), x = 65, y =max(Ncatch.i)*0.7  ) 
     text (paste ( "N (<29cm) = ",(sum(Ncatch.i[1:29] ) )  ), x = 65, y =max(Ncatch.i)*0.5 ) 
     text (paste ( "N (<27cm) = ",(sum(Ncatch.i[1:27] ) )  ), x = 65, y =max(Ncatch.i)*0.3 ) 
     text (paste ( "N (<25cm) = ",(sum(Ncatch.i[1:25] ) )  ), x = 65, y =max(Ncatch.i)*0.1 ) 
     
     N_total <- sum(Ncatch.i)
     N_below31 <- sum(Ncatch.i[1:31]  )
     N_below29 <- sum(Ncatch.i[1:29]  )
     N_below27 <- sum(Ncatch.i[1:27]  )
     N_below25 <- sum(Ncatch.i[1:25]  )
     N_below23 <- sum(Ncatch.i[1:23]  )
     
     
     
     # month 2
   
      
     Month = 2
     i = Month
     
     Ncatch.i  <-  lfq2new$catch [,i]
     Ncatch.NAi  <-  replace(Ncatch.i , Ncatch.i == 0  , NA) 
     Ncatch.NAi
     
     plot(   (Ncatch.NAi) ~lfq2new$midLengths ,
             main = paste(lfq2new$dates[Month], "no zeros shown", 
                          "plot " , i) )
     abline(v = c(25, 27, 29, 31))
     text ("25, 27, 29, 31 cm", x = 15, y =max(Ncatch.i)*0.9) 
     text (paste ( "N = ",(sum(Ncatch.i) )  ), x = 65,  y =max(Ncatch.i)*0.9 ) 
     text (paste ( "N (<31cm) = ",(sum(Ncatch.i[1:31]  ) )  ), x = 65, y =max(Ncatch.i)*0.7  ) 
     text (paste ( "N (<29cm) = ",(sum(Ncatch.i[1:29] ) )  ), x = 65, y =max(Ncatch.i)*0.5 ) 
     text (paste ( "N (<27cm) = ",(sum(Ncatch.i[1:27] ) )  ), x = 65, y =max(Ncatch.i)*0.3 ) 
     text (paste ( "N (<25cm) = ",(sum(Ncatch.i[1:25] ) )  ), x = 65, y =max(Ncatch.i)*0.1 ) 
     
     N_total <- sum(Ncatch.i)
     N_below31 <- sum(Ncatch.i[1:31]  )
     N_below29 <- sum(Ncatch.i[1:29]  )
     N_below27 <- sum(Ncatch.i[1:27]  )
     N_below25 <- sum(Ncatch.i[1:25]  )
     N_below23 <- sum(Ncatch.i[1:23]  )
       
     # START LOOP -----------------------
      
      ivec <- rep(NA, 29)
      monthvec <- rep(NA, 29)
       N_totalvec <- rep(NA, 29)
       N_below31vec <- rep(NA, 29)
       N_below29vec <- rep(NA, 29)
       N_below27vec <- rep(NA, 29)
       N_below25vec <- rep(NA, 29)
       N_below23vec <- rep(NA, 29)
       
       
       outp.table <- data.frame( i_vec = ivec,
                                month_vec = monthvec,
                                 N_total = N_totalvec,
                                N_below31= N_below31vec,
                                N_below29= N_below29vec,
                                N_below27= N_below27vec,
                                N_below25= N_below25vec,
                                N_below23= N_below23vec)
                                
       
        
       for (i in 1:29) {
     
       Ncatch.i  <-  lfq2new$catch [,i]
       Ncatch.NAi  <-  replace(Ncatch.i , Ncatch.i == 0  , NA) 
       Ncatch.NAi
       
       plot(   (Ncatch.NAi) ~lfq2new$midLengths ,
               main = paste(lfq2new$dates[i], 
                            "no zeros shown","plot " , i  ) )
       abline(v = c(25, 27, 29, 31))
       text ("25, 27, 29, 31 cm", x = 15, y =max(Ncatch.i)*0.9) 
       text (paste ( "N = ",(sum(Ncatch.i) )  ), x = 65,  y =max(Ncatch.i)*0.9 ) 
       text (paste ( "N (<31cm) = ",(sum(Ncatch.i[1:31]  ) )  ), x = 65, y =max(Ncatch.i)*0.7  ) 
       text (paste ( "N (<29cm) = ",(sum(Ncatch.i[1:29] ) )  ), x = 65, y =max(Ncatch.i)*0.5 ) 
       text (paste ( "N (<27cm) = ",(sum(Ncatch.i[1:27] ) )  ), x = 65, y =max(Ncatch.i)*0.3 ) 
       text (paste ( "N (<25cm) = ",(sum(Ncatch.i[1:25] ) )  ), x = 65, y =max(Ncatch.i)*0.1 ) 
       
       N_total <- sum(Ncatch.i)
       N_below31 <- sum(Ncatch.i[1:31]  )
       N_below29 <- sum(Ncatch.i[1:29]  )
       N_below27 <- sum(Ncatch.i[1:27]  )
       N_below25 <- sum(Ncatch.i[1:25]  )
       N_below23 <- sum(Ncatch.i[1:23]  )
      
       
       outp.table[i,1] <- i
       outp.table[i,2] <- paste (lfq2new$dates[i])
       outp.table[i,3] <-         N_total 
       outp.table[i,4] <-        N_below31 
       outp.table[i,5] <-        N_below29 
       outp.table[i,6] <-        N_below27 
       outp.table[i,7] <-       N_below25 
       outp.table[i,8] <-        N_below23 
       
        
       
       }
       
       outp.table
      
       
 # PERCENT - calculate percent N per month
     
       
 outp.table$N_below23_perc <- (outp.table$N_below23/outp.table$N_total)*100
 outp.table$N_below25_perc <- (outp.table$N_below25/outp.table$N_total)*100
 outp.table$N_below27_perc <- (outp.table$N_below27/outp.table$N_total)*100
 outp.table$N_below29_perc <- (outp.table$N_below29/outp.table$N_total)*100
 
       
       
       # plot outputs, N by month -------------
       
       outp.table$DATE <- as.Date(outp.table$month_vec, format = "%Y-%m-%d")   #Date format needs to be input like that
       
       plot(outp.table$i_vec ~ as.numeric(outp.table$DATE) )
       plot(outp.table$i_vec ~ (outp.table$DATE) )
       
       outp.table
       
     #barplots
      barplot(outp.table$i_vec)  
      barplot(outp.table$N_total, main = "N by month")
      barplot(outp.table$N_below25, main = "N<25 cm by month")
      
      # x-y plots
      lines(outp.table$N_total ~outp.table$DATE , main = "N by month")
      lines(outp.table$N_below25~outp.table$DATE, main = "N<25 cm by month")
     
      plot(outp.table$N_below25~outp.table$DATE, main = "N<25 cm by month")
        lines(outp.table$N_below25~outp.table$DATE, main = "N<25 cm by month")
      
    
        plot(outp.table$N_below25~outp.table$DATE, main = "N<25 cm by month")
        lines(outp.table$N_below25~outp.table$DATE, main = "N<25 cm by month")
      
        
        lines(outp.table$N_below23~outp.table$DATE, main = "N<25 cm by month")
        
        # FOR PAPER , N < 23 cm per month---------- 
        plot(outp.table$N_below23~outp.table$DATE, 
               main = "N < 23 cm per month",
               type = "b", pch = 19, 
               col = "red", xlab = "x", ylab = "y")
     
          
          
          
          # FOR PAPER , N < 23 cm per month---------- 
          
          plot(outp.table$N_below23~outp.table$DATE, 
               main = "N < 23 cm per month",
               type = "b", pch = 19, 
               col = "red", xlab = "x", ylab = "y")
          
          
          # FOR PAPER , N < 25 cm per month---------- 
          
          plot(outp.table$N_below25~outp.table$DATE, 
               main = "N < 25 cm per month",
               type = "b", pch = 19, 
               col = "red", xlab = "x", ylab = "y")
          text(outp.table$DATE[2],
               x = 300+outp.table$DATE[2],
               y = 63)
          text(outp.table$DATE[28],
               x = (outp.table$DATE[28]-80),
               y = 55)
          text(outp.table$DATE[21],
               x = outp.table$DATE[21],
               y = 68)
          
          text(outp.table$DATE[19],
               x = outp.table$DATE[19],
               y = 73)
          
          text(outp.table$DATE[16],
               x = outp.table$DATE[16],
               y = 57)
          
          # FOR PAPER , N < 27 cm per month---------- 
          
          plot(outp.table$N_below27~outp.table$DATE, 
               main = "N < 27 cm per month",
               type = "b", pch = 19, 
               col = "red", xlab = "x", ylab = "y")
          text(outp.table$DATE[2],
               x = 300+outp.table$DATE[2],
               y = 63)
          text(outp.table$DATE[28],
               x = (outp.table$DATE[28]-80),
               y = 55)
          text(outp.table$DATE[21],
               x = outp.table$DATE[21],
               y = 68)
          
          text(outp.table$DATE[19],
               x = outp.table$DATE[19],
               y = 73)
          
          text(outp.table$DATE[16],
               x = outp.table$DATE[16],
               y = 57)
          
          # FOR PAPER , N < 29 cm per month---------- 
          
          plot(outp.table$N_below29~outp.table$DATE, 
               main = "N < 29 cm per month",
               type = "b", pch = 19, 
               col = "red", xlab = "x", ylab = "y")
          # text(outp.table$DATE[2],
          #      x = 300+outp.table$DATE[2],
          #      y = 63)
          # text(outp.table$DATE[28],
          #      x = (outp.table$DATE[28]-80),
          #      y = 55)
          # text(outp.table$DATE[21],
          #      x = outp.table$DATE[21],
          #      y = 68)
          # 
          # text(outp.table$DATE[19],
          #      x = outp.table$DATE[19],
          #      y = 73)
          # 
          # text(outp.table$DATE[16],
          #      x = outp.table$DATE[16],
          #      y = 57)
          
          
          
          
          # FOR PAPER several lines 27 25 23 cm----------
         
          plot(outp.table$N_below27~outp.table$DATE, 
               main = "N < 27, 25, 23 cm per month",
               type = "b", pch = 19, 
               col = "red", xlab = "x", ylab = "y")
          points(outp.table$N_below25~outp.table$DATE, 
                 pch = 6,  type = "b",  
                 col = "darkgreen"   )  
                 points (  outp.table$N_below23~outp.table$DATE, 
                 pch = 4,  type = "b",  
                 col = "navy"   )  
        
    # FOR PAPER lines only, N < 27, 25, 23 cm per month ------------
                  
                  plot(outp.table$N_below27~outp.table$DATE, 
                       main = "N < 27, 25, 23 cm per month",
                       type = "l", 
                       xlab = "Year",
                       ylab = "N, indiv. per month",
                       col = "red",
                       ylim = c(0, 200))
                  
                  lines (  outp.table$N_below25~outp.table$DATE, 
                            col = "darkgreen"   )  
                  lines (  outp.table$N_below23~outp.table$DATE, 
                            col = "blue"   )  
                  lines (  outp.table$N_below27~outp.table$DATE, 
                           col = "red"   )  
                  legend("topleft", 
                         legend=c("N < 23cm", "N < 25cm", "N < 27cm"),
                         col=c( "blue","darkgreen" , "red"), 
                         lty=1, cex=0.8)
        
                  
        # CONCLUSION   (N small-sized indiv.  per month)       
        # Weak or negligible recruitment in 2005, 2006 and 2007!      
        # strong recruitment in November 2004 and from march 2008!      
        # Seasonal recruitment peaks mostly in October-November, 
        # (except  in 2008, where there were peaks throughout the year),
        # and in 2006 (negligible recruitment throughout the year).
          
                  
# PERCENT small ones
# FOR PAPER lines only, PERCENT < 27, 25, 23 cm per month ------------
                  
                  plot(outp.table$N_below27_perc ~outp.table$DATE, 
                       main = "% N per month, < 27, < 25, and < 23 cm ",
                       type = "l", 
                       xlab = "Year",
                       ylab = "Pecent of total N, per month (%)",
                       col = "red",
                       ylim = c(0, 60))
                  
                  lines (  outp.table$N_below25_perc~outp.table$DATE, 
                           col = "darkgreen"   )  
                  lines (  outp.table$N_below23_perc~outp.table$DATE, 
                           col = "blue"   )  
                  lines (  outp.table$N_below27_perc~outp.table$DATE, 
                           col = "red"   )  
                  legend("topleft", 
                         legend=c("% N < 23cm", "% N < 25cm", "% N < 27cm"),
                         col=c( "blue","darkgreen" , "red"), 
                         lty=1, cex=0.8)
                  
                  
# CONCLUSION (% small ones)         
# negligible recruitment in 2006
# strong recruitment  from march 2008!      
# Seasonal recruitment peaks mostly in October-November, 
 # (except  in 2007 and 2008, where there were peaks throughout the year,
 # and in 2006, when there was negligible recruitment throughout the year)
                  
                  
             
 # FOR PAPER Catch curve by year --------------- 
                  
# Catch curve by year --------------- 
                  
                  # ALL YEARS ---------
               barplot(Ncatch.NA ~lfq2new$midLengths , 
               main = "catch curve, S stock, N = 1399",
               col = "lightblue")
                  
                  
                  plot(   log10(Ncatch) ~lfq2new$midLengths  , main = "log scale")
                  plot(   Ncatch ~lfq2new$midLengths  )
                  
                  
                  
                  lfq2new
                  
                  summary(lfq2$DATE)
                  
                  summary(lfq2$Length_CM)
                  
                  length(lfq2$Length_CM)
                  
                  # extract years-------------------
                  
                  
                  
                  
                  
                  # FOR PAPER LFDs ----------------
                  # plot raw and restructured LFQ data --------------------------
                  
                  # plot raw and restructured LFQ data bs = 2 --------------------
                  
                  abline(h = 24)
                  abline(h = 26)
                  abline(h = 30)
                  
                  # plot small ones only 
                  
                  lfq2$Frequency
                  lfq2$DATE
                  lfq2new$Length_CM
                  
                  lfq2$Length_CM[1:46] # 46 length classes , fom 16 to 61 cm
                  
                  
                  
                  barplot(   Ncatch.NA ~lfq2new$midLengths , 
                main = "Catch Curve, all years, S stock, N = 20,090 indiv.",
                col = "lightblue")
     
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
  # FOR PAPER LFDs  ---------------
                  
                library(TropFishR)
                
                # load data ----------------   
                lfq2 <- read.csv("stbrasLF-trunc.csv")
                # View(lfq2)
                  
                summary(lfq2)
                  
                  plot(log10(lfq2$Frequency) ~lfq2$Length_CM)
                  
                  names(lfq2)
                  lfq2$DATE <- as.Date(lfq2$DATE, format = "%d/%m/%Y") #Date format needs to be input like that
                  
                  summary(lfq2)# 12 to 79 cm
                  
                  
lfq2new <- lfqCreate(data = lfq2, Lname = "Length_CM", Dname = "DATE", 
                     Fname = "Frequency")
                  
                  plot(lfq2new, Fname = "catch")
                  
                  # adjust bin size (on the histogram) #Bin size 2!
                  synLFQ7a <- lfqModify(lfq2new, bin_size = 2) #mid-sizes are already at 0.5
    # MA3
                                
 lfqbinMA3 <- lfqRestructure(synLFQ7a, MA = 3, addl.sqrt = TRUE) #MA = moving average sample size 
  opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
opar2 <- par()
     
 plot(lfqbinMA3, Fname = "catch", date.axis = "modern", main = "MA3")
plot(lfqbinMA3, Fname = "rcounts", date.axis = "modern")
  
# MA5

lfqbinMA5 <- lfqRestructure(synLFQ7a, MA = 5, addl.sqrt = TRUE) #MA = moving average sample size 
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
opar2 <- par()

plot(lfqbinMA5, Fname = "catch", date.axis = "modern", main = "MA5 ")
plot(lfqbinMA5, Fname = "rcounts", date.axis = "modern")

                  par(opar)
                  par(opar2)
                  
dev.off()     
              
                  
# ALL YEARS (Catch Curve) ---------
                  barplot(Ncatch.NA ~lfq2new$midLengths , 
                             main = "catch curve, S stock, N = 1399",
                             col = "lightblue")
                  
           