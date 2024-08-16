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
# Steenbras71
# Seed = Sys.time! 
# Exps, 1,2,3,4
# MA 3 ,  MA5, Linfmin 50 cm ("naive" wide search space)
# Linf: 50 to 200 cm, seasonal, bs  = 2
# MA5, and MA3
#Linfmin 50 and 30 cm (widened search space. after tag-recap results)
#  bs  = 2


# Lithognatus aureti growth
# Southern stock
# v3d 
# March 2024
# Ralf Schwamborn: ralf.schwamborn@ufpe.br

# Data: 
# Fish lengths (fork length) from 2004-08-25 to 2009-12-30 
# (approx. five and a half years) 
# Length range: 12 cm to 79 cm
# N = 1399 indiv.


# 1. Initialization ----------------------------------------------------------

# 1.1. required packages =====
# clean memory
gc()  
gc(reset=T)
rm(list = ls())

library(parallel)
library(TropFishR)
library(fishboot)
library(ks)

help(package= "TropFishR")





# BS = 2 -----------------------------------------------------------
# adjust bin size (on the histo

# 1.2. Load data =====

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

lfq2 <- read.csv("stbrasLF-trunc.csv", header = T)

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


# 1.3 plot raw and restructured LFQ data --------------------------

# 1.3a plot raw and restructured LFQ data bs = 2 --------------------

#Bin size 2! -----------------
synLFQ7a <- lfqModify(lfq, bin_size = 2) #mid-sizes are already at 0.5


# MA 3
lfqbin <- lfqRestructure(synLFQ7a, MA = 3, addl.sqrt = TRUE) #MA = moving average sample size 

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "catch", date.axis = "modern")
plot(lfqbin, Fname = "rcounts", date.axis = "modern", main = "MA = 3")
par(opar)


# MA 5
lfqbin <- lfqRestructure(synLFQ7a, MA = 5, addl.sqrt = TRUE) #MA = moving average sample size 
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "catch", date.axis = "modern")
plot(lfqbin, Fname = "rcounts", date.axis = "modern", main = "MA = 5")
par(opar)



 #Bin size 2! -----------------
synLFQ7a <- lfqModify(lfq, bin_size = 2) #mid-sizes are already at 0.5


# plot with different MA values, MA 3,5,7,9,11

lfqbin <- lfqRestructure(synLFQ7a, MA = 11, addl.sqrt = TRUE) #MA = moving average sample size 

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "rcounts", date.axis = "modern", main = "Restructured Data , bs = 2, MA = 11")
plot(lfqbin, Fname = "catch", date.axis = "modern", main = "Raw Data, bs = 2")

par(opar)


lfqbin <- lfqRestructure(synLFQ7a, MA = 9, addl.sqrt = TRUE) #MA = moving average sample size 

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "rcounts", date.axis = "modern", main = "Restructured Data , bs = 2, MA = 9")
plot(lfqbin, Fname = "catch", date.axis = "modern", main = "Raw Data, bs = 2")

par(opar)

lfqbin <- lfqRestructure(synLFQ7a, MA = 7, addl.sqrt = TRUE) #MA = moving average sample size 

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "rcounts", date.axis = "modern", main = "Restructured Data , bs = 2, MA = 7")
plot(lfqbin, Fname = "catch", date.axis = "modern", main = "Raw Data, bs = 2")

par(opar)

lfqbin <- lfqRestructure(synLFQ7a, MA = 5, addl.sqrt = TRUE) #MA = moving average sample size 

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "rcounts", date.axis = "modern", main = "Restructured Data , bs = 2, MA = 5")
plot(lfqbin, Fname = "catch", date.axis = "modern", main = "Raw Data, bs = 2")

par(opar)

lfqbin <- lfqRestructure(synLFQ7a, MA = 3, addl.sqrt = TRUE) #MA = moving average sample size 
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "rcounts", date.axis = "modern", main = "Restructured Data , bs = 2, MA = 3")
plot(lfqbin, Fname = "catch", date.axis = "modern", main = "Raw Data, bs = 2")

par(opar)

# Conclusion form looking at the pĺots: MA 3 is messy, with no consistent modal progression... best are MA 5 to 9

# MA = 5 gives more uncertainty... 


# USED for Analysis = SLOW, very wide Linf, search parameters: ----------------------------

# Exp 1 and 2 
## MA <- 3
# MA <- 5 #5 is best? MA = 3 is more realistic, MA = 5 gives less uncertainty...
# low_par <- list(Linf = 50, K = 0.01, t_anchor = 0, C = 0, ts = 0)
# up_par <- list(Linf = 200, K = 2, t_anchor = 1, C = 1, s = 1)
# seasonalised <- TRUE
# popSize <- 48
# maxiter <- 100 #30 but start with 10
# run <- 40
# pmutation <- 0.2
# nresamp <- 200 #Number of resamplings; 20 but start with 5


# 
# # Exp 3 ----------
# # USED for Analysis = SLOW, very wide Linf, search parameters: ----------------------------
# MA <- 5 # MA = 3 is more realistic, MA = 5 gives less uncertainty...
# low_par <- list(Linf = 30, K = 0.01, t_anchor = 0, C = 0, ts = 0)
# up_par <- list(Linf = 200, K = 2, t_anchor = 1, C = 1, s = 1)
# seasonalised <- TRUE
# popSize <- 48
# maxiter <- 100 #30 but start with 10
# run <- 40
# pmutation <- 0.2
# nresamp <- 200 #Number of resamplings; 20 but start with 5
# 
# 

# Exp 4 ----------
# 
# #SLOW, very wide Linf, search parameters: ----------------------------
# MA <- 5 # MA = 3 is more realistic, MA = 5 gives less uncertainty...
# low_par <- list(Linf = 30, K = 0.01, t_anchor = 0, C = 0, ts = 0)
# up_par <- list(Linf = 200, K = 2, t_anchor = 1, C = 1, s = 1)
# seasonalised <- TRUE
# popSize <- 48
# maxiter <- 100 #30 but start with 10
# run <- 40
# pmutation <- 0.2
# nresamp <- 1000 #Number of resamplings; 20 but start with 5
# 
# 



# Used for actual analysis ----------

# USED for Analysis = SLOW, very wide Linf, search parameters: ----------------------------
MA <- 3 # MA = 3 is more realistic, MA = 5 gives less uncertainty...
low_par <- list(Linf = 50, K = 0.01, t_anchor = 0, C = 0, ts = 0)
up_par <- list(Linf = 200, K = 2, t_anchor = 1, C = 1, s = 1)
seasonalised <- TRUE
popSize <- 48
maxiter <- 100 #30 but start with 10
run <- 40
pmutation <- 0.2
nresamp <- 300 #Number of resamplings; 200 but start with 5



#  full bootstrap GA, runs I =====

set.seed(Sys.time())
t1 <- Sys.time()
res_seed_Sys_time_exp8A_MA3_n300_Linfmin50_8c <- ELEFAN_GA_boot(lfq=synLFQ7a, MA = 3,
                                                               seasonalised = seasonalised,
                                                               up_par = up_par, low_par = low_par,
                                                               popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
                                                               nresamp = 300, parallel = TRUE, no_cores = detectCores(),
                                                               seed = as.numeric(Sys.time()), resample = TRUE
)
t2 <- Sys.time()
t2 - t1
# 
# apprx 1 minute per run, 300 min from 14:40 = 
# 300 /60 = 5 hs
# 5 hrs from 14:30 = 19:30 hrs estimated end time

head(res_seed_Sys_time_exp8A_MA3_n300_Linfmin50_8c$bootRaw)

# Time difference of 2.6 hours for 300 runs!

# 
# t2 <- Sys.time()
# > t2 - t1
# Time difference of 59.44136 mins
# > 
#   > 
#   > 
#   > head(res_seed_Sys_time_exp6_MA5_n300_Linfmin30_6b$bootRaw)
# Linf          K  t_anchor         C        ts     phiL
# 1  46.34145 0.68975344 0.4881445 0.5776292 0.3859077 3.170633
# 2 103.28955 0.07412741 0.8015209 0.6193316 0.4179676 2.898092
# 3  46.86904 0.65887986 0.4620354 0.3332298 0.6192563 3.160578
# 4  44.09055 0.91736495 0.7592372 0.4745363 0.4944776 3.251233
# 5  44.12041 0.58028093 0.3823612 0.5649906 0.5883433 3.052917


# define working directory

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Save the results on the hard drive --------------------------

#  t2 - t1
# Time difference of 51 mins for 100 runs
#  1.5 hs for 200 runs (MA)


### SAVE results ---------------
write.table(res_seed_Sys_time_exp8A_MA3_n300_Linfmin50_8c$bootRaw  ,  
            row.names = FALSE, sep = ";",
            file = "res_seed_Sys_time_exp8A_MA3_n300_Linfmin50_8c.csv"  )



#  full bootstrap GA, runs II =====

low_par <- list(Linf = 50, K = 0.01, t_anchor = 0, C = 0, ts = 0)

set.seed(Sys.time())
t1 <- Sys.time()
res_seed_Sys_time_expB9_MA5_n200_Linfmin50_9a <- ELEFAN_GA_boot(lfq=synLFQ7a, MA = 5,
                                                             seasonalised = seasonalised,
                                                             up_par = up_par, low_par = low_par,
                                                             popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
                                                             nresamp = 200, parallel = TRUE, no_cores = detectCores(),
                                                             seed = as.numeric(Sys.time()), resample = TRUE
)
t2 <- Sys.time()
t2 - t1



head(res_seed_Sys_time_expB9_MA5_n200_Linfmin50_9a$bootRaw)
# 

# define working directory

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)


### SAVE results ---------------
write.table(res_seed_Sys_time_expB9_MA5_n200_Linfmin50_9a$bootRaw  ,  
            row.names = FALSE, sep = ";",
            file = "res_seed_Sys_time_expB9_MA5_n200_Linfmin50_9a.csv"  )






# Used for actual analysis ----------

# USED for Analysis = SLOW, very wide Linf, search parameters: ----------------------------
MA <- 3 # MA = 3 is more realistic, MA = 5 gives less uncertainty...
low_par <- list(Linf = 30, K = 0.01, t_anchor = 0, C = 0, ts = 0)
up_par <- list(Linf = 200, K = 2, t_anchor = 1, C = 1, s = 1)
seasonalised <- TRUE
popSize <- 48
maxiter <- 100 #30 but start with 10
run <- 40
pmutation <- 0.2
nresamp <- 200 #Number of resamplings; 200 but start with 5



#  full bootstrap GA, runs III =====


low_par <- list(Linf = 30, K = 0.01, t_anchor = 0, C = 0, ts = 0)


set.seed(Sys.time())
t1 <- Sys.time()
res_seed_Sys_time_expC10_MA5_n200_Linfmin30_10 <- ELEFAN_GA_boot(lfq=synLFQ7a, MA = 5,
                                                               seasonalised = seasonalised,
                                                               up_par = up_par, low_par = low_par,
                                                               popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
                                                               nresamp = 200, parallel = TRUE, no_cores = detectCores(),
                                                               seed = as.numeric(Sys.time()), resample = TRUE
)
t2 <- Sys.time()
t2 - t1



head(res_seed_Sys_time_expC10_MA5_n200_Linfmin30_10$bootRaw)



# define working directory

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)


### SAVE results ---------------
write.table(res_seed_Sys_time_expC10_MA5_n200_Linfmin30_10$bootRaw  ,  
            row.names = FALSE, sep = ";",
            file = "res_seed_Sys_time_expC10_MA5_n200_Linfmin30_10.csv"  )


#  full bootstrap GA, runs IV =====


low_par <- list(Linf = 30, K = 0.01, t_anchor = 0, C = 0, ts = 0)


set.seed(Sys.time())
t1 <- Sys.time()
res_seed_Sys_time_expD11_MA3_n200_Linfmin30_11 <- ELEFAN_GA_boot(lfq=synLFQ7a, MA = 3,
                                                                 seasonalised = seasonalised,
                                                                 up_par = up_par, low_par = low_par,
                                                                 popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
                                                                 nresamp = 200, parallel = TRUE, no_cores = detectCores(),
                                                                 seed = as.numeric(Sys.time()), resample = TRUE
)
t2 <- Sys.time()
t2 - t1


head(res_seed_Sys_time_expD11_MA3_n200_Linfmin30_11$bootRaw)
# 

# define working directory

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)




### SAVE results ---------------
write.table(res_seed_Sys_time_expD11_MA3_n200_Linfmin30_11$bootRaw  ,  
            row.names = FALSE, sep = ";",
            file = "res_seed_Sys_time_expD11_MA3_n200_Linfmin30_11.csv"  )

# 
# 
# write.table(res_seed_Sys_time_exp5_MA3_n300_Linfmin30_5d$bootRaw  ,  
#             row.names = FALSE, sep = ";",
#             file = "res_seed_Sys_time_exp5_MA3_n300_Linfmin30_5d.csv"  )


# 
# write.table(res_seed_Sys_time_exp5_MA3_Linfmin30_5b$bootRaw  ,  
#             row.names = FALSE, sep = ";",
#             file = "res_seed_Sys_time_exp5_MA3_Linfmin30_5b.csv"  )
# 
# 
# write.table(res_seed_Sys_time_exp5_MA3_Linfmin30_5a$bootRaw  ,  
#             row.names = FALSE, sep = ";",
#             file = "res_seed_Sys_time_exp5_MA3_Linfmin30_5a.csv"  )

# res_seed1 # seed = 1
# 
# bootRaw
# Linf          K  t_anchor         C        ts     phiL
# 1 73.43908 0.33525641 0.5012003 0.2671175 0.3880342 3.257232
# 2 79.98777 0.13707359 0.4779149 0.3822207 0.1032816 2.943001
# 3 85.40818 0.08786159 0.3577160 0.9473927 0.1386260 2.806798
# 
# $seed
# [1] 2 3 4
# 
# 
# # seed = NULL 
# 
# res_seed_null_exp5
# 
# $bootRaw
# Linf          K  t_anchor         C        ts     phiL
# 1  80.36273 0.18170474 0.1497518 0.4121988 0.8830967 3.069476
# 2 132.82223 0.04219826 0.4361175 0.1123161 0.7172888 2.871836
# 3  85.76286 0.10497810 0.4543345 0.6175214 0.5183198 2.887697
# 
# $seed
# [1] 125272 125273 125274
# 
# attr(,"class")
# [1] "lfqBoot"


# 
# head(res_seed_Sys_time_exp5_MA3_Linfmin30_5a$bootRaw)
# Linf          K  t_anchor         C        ts     phiL
# 1  71.32504 0.24624722 0.8860452 0.6956757 0.4020117 3.097855
# 2 116.07617 0.06148727 0.5333962 0.3958513 0.1144714 2.918271
# 3  66.97340 0.40928635 0.6937918 0.8781590 0.6480501 3.263832
# 4  77.62680 0.19549713 0.5497919 0.3733314 0.2820126 3.071164
# 5  77.64766 0.19027992 0.5097190 0.5377912 0.3760502 3.059650
# 6 181.57447 0.02867642 0.2025031 0.6074693 0.7272168 2.975634
# > 
#   > head( res_seed_Sys_time_exp5_MA3_Linfmin30_5b$bootRaw)
# Linf          K  t_anchor         C        ts     phiL
# 1 85.68007 0.11693664 0.5064350 0.6958852 0.7310073 2.933710
# 2 75.17222 0.16641928 0.7599001 0.7584347 0.3280506 2.973318
# 3 94.23982 0.06465385 0.3886676 0.6194979 0.4452537 2.759063
# 4 87.76946 0.11827253 0.5434103 0.4820488 0.6648035 2.959571
# 5 78.51491 0.13166529 0.6949723 0.4434669 0.5596474 2.909376
# 6 76.35249 0.12061858 0.3590052 0.1817466 0.3950559 2.847061
# 





res_EXP5_v59_S_MA3_slowLinf30_200_n100_v59_5a <- res



#res_EXP3_S_MA5_slowLinf30_200_n200_3c <- res



# res_EXP_S_slow_n100_1b <- res


# define working directory

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)


# Save the results on the hard drive --------------------------
# Export (save) Results files ----------------
 
# 2.1 MA3 ----------
 
 # Experiment 1-------------------
# MA3, Southern stock, seasonal, slow, Males and Females, Linf 50 to 200 cm. nruns = 1000
#write.table(file = "res_EXP_S_slow_n100_1a.csv", res_EXP_S_slow_n100_1a$boot   )
#write.table(file = "res_EXP_S_slow_n100_1b.csv", res_EXP_S_slow_n100_1b$boot   )
#write.table(file = "res_EXP_S_slow_n200_1c.csv", res_EXP_S_slow_n200_1c$boot   )
#write.table(file = "res_EXP_S_slow_n200_1d.csv", res_EXP_S_slow_n200_1d$boot   )
#write.table(file = "res_EXP_eS_slow_n200_1e.csv", res_EXP_S_slow_n200_1e$boot   )
#write.table(file = "res_EXP_eS_slow_n200_1f.csv", res_EXP_S_slow_n200_1f$boot   )
 
 # 2.2 MA5 ------------- 
 
 # Experiment 2-------------------
 #MA5, Southern stock, seasonal, slow, Males and Females, Linf 50 to 200 cm. nruns = 800
 
 #write.table(file = "res_EXP_S_MA5_slow_n200_2a.csv", res_EXP_S_MA5_slow_n200_2a$boot   )

 #write.table(file = "res_EXP_S_MA5_slow_n300_2b.csv", res_EXP_S_MA5_slow_n300_2b$boot   )
 
# write.table(file = "res_EXP_S_MA5_slow_n300_2c.csv", res_EXP_S_MA5_slow_n300_2c$boot   )
 
 # Time difference of 2.09 hours for 300 runs (res_EXP_S_MA5_slow_n300_2c.csv)
 
 # write.table(file = "res_EXP_S_MA5_slow_n200_2d.csv", res_EXP_S_MA5_slow_n200_2d$boot   )
 
 
 
 
 # 2.3 MA5, 30 cm to 200 cm ----------
 

 # Experiment 3 -------------------
 
 
 #MA3, Southern stock, seasonal, slow, Males and Females, Linf 30 to 200 cm. nruns = 1000
 
 
 #write.table(file = "res_EXP3_S_MA5_slowLinf30_200_n300_3a.csv", 
 #      res_EXP3_S_MA5_slowLinf30_200_n500_3b$boot   )
 
# Time difference of 1.073997 hours
 
 
 #write.table(file = "res_EXP3_S_MA5_slowLinf30_200_n500_3b.csv", 
       #      res_EXP3_S_MA5_slowLinf30_200_n500_3b$boot   )
 
 #Time difference of 1.825781 hours
 
 
 # write.table(file = "res_EXP3_S_MA5_slowLinf30_200_n200_3c.csv", 
 #             res_EXP3_S_MA5_slowLinf30_200_n200_3c$boot   )
 # 
 
 
 
 # Time difference of ...
 
 
 
 # Experiment 4 -------------------
 
 
 #MA3, Southern stock, seasonal, slow, Males and Females, Linf 30 to 200 cm. nruns = 1000
 
 
 #write.table(file = "res_EXP3_S_MA5_slowLinf30_200_n300_3a.csv", 
 #      res_EXP3_S_MA5_slowLinf30_200_n500_3b$boot   )
 
 # Time difference of 1.073997 hours
 
 
 #write.table(file = "res_EXP3_S_MA5_slowLinf30_200_n500_3b.csv", 
 #      res_EXP3_S_MA5_slowLinf30_200_n500_3b$boot   )
 
 #Time difference of 1.825781 hours
 
 
 write.table(file = "res_EXP4_S_MA3_slowLinf30_200_n1000_4a.csv", 
             res_EXP4_S_MA3_slowLinf30_200_n1000_4a$boot   )
 
 # Time difference of ...
 
 
 
 
# retrieve the results from the hard drive --------------------------


test_object1b <-  read.table(file = "res_EXP_S_slow_n100_1a.csv" )

# use *.table, not *csv

# import the results into a new fishboot results object

res1b <- res


res1b$bootRaw <- test_object1b

# plot the imported results


# univariate density plot of bootstrapped pars
univariate_density(res1b)

# Linf / K scatterhist GA
LinfK_scatterhist(res1b, Linf.breaks=50, K.breaks=50, 
                  phi.contour = TRUE, phi.contour.lty = 1, phi.contour.col = 8,
                  pt.col = rgb(0,0,0,0.1), pt.pch = "."
)


# VBGF by time growth curve plot
CIinfo <- vbgfCI_time(
  res = res1b, 
  agemax = 12, CI = 95,
  perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
  ci.col = 1, ci.lty = 2, ci.lwd = 1, 
  maxd.col = 1, maxd.lty = 1, maxd.lwd = 2
)

# results:
#MA = 3 : higher Linf 74 to 120cm, slower growth  K 0.09 to 0.3, more realistic, but higher incertainty
# MA = 5 : lower Linf, 43 to 46 cm, faster growth, K approx. 0.8 y-1, narrower, less uncertainty 

# MA = 3
# Linf          K  t_anchor         C        ts     phiL
# 1  74.46511 0.31686723 0.3849602 0.1826848 0.5715933 3.244783
# 2  97.12960 0.09776887 0.4819707 0.5141707 0.7262876 2.964904
# 3 119.64557 0.09058452 0.5931490 0.6456655 0.2617472 3.112847
# 4  78.82207 0.28225487 0.3524015 0.3011222 0.5772155 3.243937


# # MA = 5
# Linf         K  t_anchor         C        ts     phiL
# 1 45.50524 0.8144299 0.7774926 0.3927460 0.5728192 3.226977
# 2 43.18455 0.7232855 0.7185665 0.6298900 0.7355327 3.129967
# 3 45.12517 0.7783987 0.4759445 0.4031592 0.3337805 3.200040
# 4 44.54400 0.7826411 0.4830835 0.4665224 0.4104989 3.191141
# 





# plot results GA
# univariate density plot of bootstrapped pars
univariate_density(res)

# Linf / K scatterhist GA
LinfK_scatterhist(res, Linf.breaks=50, K.breaks=50, 
                  phi.contour = TRUE, phi.contour.lty = 1, phi.contour.col = 8,
                  pt.col = rgb(0,0,0,0.1), pt.pch = "."
)


# VBGF by time growth curve plot
CIinfo <- vbgfCI_time(
  res = res, 
  agemax = 12, CI = 95,
  perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
  ci.col = 1, ci.lty = 2, ci.lwd = 1, 
  maxd.col = 1, maxd.lty = 1, maxd.lwd = 2
)



##############

# 3. Try other MA values -------------------------------
# BS = 2, fixed, variable MA

# adjust bin size (on the histogram) #Bin size 2! -----------------
synLFQ7a <- lfqModify(lfq, bin_size = 2) #mid-sizes are already at 0.5


# 3.1 MA vs time (fast settigs, wide Linf) -----------------------------

# FAST SETTINGS
# 3.2 USED for MA vs time Analysis = FAST ,  wide Linf, search parameters: ----------------------------

# TEST 1: MA3,  medium fast , medium wide, seasonal 
MA <- 3 #5 is best? MA = 3 is more realistic, MA = 5 gives less uncertainty...
low_par <- list(Linf = 50, K = 0.01, t_anchor = 0, C = 0, ts = 0)
up_par <- list(Linf = 180, K = 1, t_anchor = 1, C = 1, s = 1)
seasonalised <- TRUE
popSize <- 48
maxiter <- 50 #30 but start with 10
run <- 12
pmutation <- 0.2
nresamp <- 5 #Number of resamplings; 20 but start with 5



# 3.3. full bootstrap GA=====

set.seed(Sys.time())
t1 <- Sys.time()
res <- ELEFAN_GA_boot(lfq=synLFQ7a, MA = 3, seasonalised = seasonalised,
                      up_par = up_par, low_par = low_par,
                      popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
                      nresamp = nresamp, parallel = TRUE, no_cores = detectCores(),
                      seed = 1, resample = TRUE
)
t2 <- Sys.time()
t2 - t1
res$boot
head(res$seed)


res


# Linf          K  t_anchor         C         ts     phiL
# 1 86.94316 0.10948186 0.5046269 0.5512178 0.42988807 2.917813
# 2 91.63775 0.07315436 0.9391105 0.6350239 0.05898373 2.788389
# 3 85.80638 0.08778147 0.6186399 0.7472966 0.46776886 2.810442
# 4 73.26899 0.21786764 0.5718505 0.3689174 0.61883055 3.068033
# 5 76.19323 0.17617203 0.7264328 0.4276772 0.67591722 3.009770


# 1.161161 mins
# 1.16/5 runs= 0.23 min



# TEST 2: MA3,  medium fast , medium wide, non-seasonal (!)

MA <- 3 #5 is best? MA = 3 is more realistic, MA = 5 gives less uncertainty...
low_par <- list(Linf = 50, K = 0.01, t_anchor = 0, C = 0, ts = 0)
up_par <- list(Linf = 180, K = 1, t_anchor = 1, C = 1, s = 1)
seasonalised <- TRUE
popSize <- 48
maxiter <- 50 #30 but start with 10
run <- 12
pmutation <- 0.2
nresamp <- 5 #Number of resamplings; 20 but start with 5

?ELEFAN_GA_boot

# 3.3. full bootstrap GA=====

set.seed(Sys.time())
t1 <- Sys.time()
res <- ELEFAN_GA_boot(lfq=synLFQ7a, MA = MA, seasonalised = FALSE,
                      up_par = up_par, low_par = low_par,
                      popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
                      nresamp = nresamp, parallel = TRUE, no_cores = detectCores(),
                      seed = 1, resample = TRUE
)
t2 <- Sys.time()
t2 - t1
res$boot
head(res$seed)


res

# > res$boot
# Linf          K  t_anchor     phiL
# 1 89.72436 0.11160207 0.3978542 2.953493
# 2 90.56294 0.07072030 0.6623725 2.763445
# 3 99.48934 0.06309123 0.7233244 2.795522
# 4 85.58678 0.15998188 0.4987202 3.068884
# 5 93.48495 0.06887650 0.5872481 2.779554

# 
# 1.0449 mins
#  1.05/5 runs = 0.21 min/run



# TEST 3: MA5,  medium fast , medium wide, non-seasonal (!) OK, BEST?

MA <- 5 # MA = 3 is more realistic, MA = 5 gives less uncertainty...
low_par <- list(Linf = 50, K = 0.01, t_anchor = 0, C = 0, ts = 0)
up_par <- list(Linf = 180, K = 1, t_anchor = 1, C = 1, s = 1)
seasonalised <- TRUE
popSize <- 48
maxiter <- 50 #30 but start with 10
run <- 12
pmutation <- 0.2
nresamp <- 5 #Number of resamplings; 20 but start with 5

 
# $bootRaw
# Linf          K  t_anchor     phiL
# 1  99.06422 0.16211881 0.6196016 3.201667
# 2  83.43285 0.17004119 0.5874285 3.073228
# 3  94.64246 0.08019843 0.6453606 2.856338
# 4  83.13154 0.17607790 0.6762717 3.085237
# 5 101.09201 0.07560878 0.6744927 2.888006


# 33.59804 secs
# 33.6/5  runs = 6.7 sec/run



# TEST 4: MA7,  medium fast , medium wide, non-seasonal (!)  very variable

MA <- 7 # MA = 3 is more realistic, MA = 5 gives less uncertainty...
low_par <- list(Linf = 50, K = 0.01, t_anchor = 0, C = 0, ts = 0)
up_par <- list(Linf = 180, K = 1, t_anchor = 1, C = 1, s = 1)
seasonalised <- TRUE
popSize <- 48
maxiter <- 50 #30 but start with 10
run <- 12
pmutation <- 0.2
nresamp <- 5 #Number of resamplings; 20 but start with 5


# #bootRaw
# Linf          K   t_anchor     phiL
# 1  80.93718 0.15834517 0.04142287 3.015901
# 2 130.31661 0.05460174 0.69304545 2.967206
# 3  52.10047 0.56645317 0.31984525 3.186847
# 4 110.07744 0.18066886 0.43241455 3.340280
# 5 102.92340 0.14924205 0.65098206 3.198919


# 33.12 secs
# 33.12/5  runs = 6.6 sec/run


# TEST 5: MA 9,  medium fast , medium wide, non-seasonal (!) MA 9 is very variable! not good!

MA <- 9 # MA = 3 is more realistic, MA = 5 gives less uncertainty...
low_par <- list(Linf = 50, K = 0.01, t_anchor = 0, C = 0, ts = 0)
up_par <- list(Linf = 180, K = 1, t_anchor = 1, C = 1, s = 1)
seasonalised <- TRUE
popSize <- 48
maxiter <- 50 #30 but start with 10
run <- 12
pmutation <- 0.2
nresamp <- 5 #Number of resamplings; 20 but start with 5


# 
# > res
# $bootRaw
# Linf          K  t_anchor     phiL
# 1 110.40374 0.09220967 0.4861076 3.050744
# 2  62.80052 0.33931059 0.9767298 3.126524
# 3 122.74237 0.07018809 0.4925202 3.024252
# 4 143.14189 0.09256290 0.4284296 3.277970
# 5 111.38397 0.18767787 0.5349684 3.367058


# 31.14287 secs secs
# 31.14287 /5 runs = 6.2 sec/run