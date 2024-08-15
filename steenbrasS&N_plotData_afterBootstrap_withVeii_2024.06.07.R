#
### Analysis of west coast steenbras (Lithugnathus aureti) body growth
### using bootstrapped length-at-age (otolith) analyses
#
# Margit Wilhelm 
# August 2024
#
###
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
# Corresponding author: mwilhelm@unam.na
#
# This script is based on from Schwammborn et al. 2018. Analysis. 
# Read in Bootstrapped OTOLITH DATA
# THEN PLOTS VBGF
# THEN LFA WITH DENSITY PLOTS
# THEN T&R WITH DENSITY PLOTS
#
# 1. Initialization -----------------------------------------------
# 1.1. required packages =====
library(fishboot)
library(TropFishR) #for VBGF
library(tidyverse) #Uses ggplot for nice plots and saving
library(latticeExtra)

#To clean memory
gc()  
gc(reset=T)
rm(list = ls())

# 2. Load data =====
# For LFD data => ONLY FOR PLOTTING THE ADDITONAL PLOTS FROM FISHBOOT LATER

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

#age.length.rawS = read.csv("age-length-rawS.csv", header=TRUE) #This is Holzhausen and Kirchner (needed for some plots) 
#age.length.rawN = read.csv("age-length-rawN.csv", header=TRUE) #This is Holzhausen and Kirchner

age.length.rawVS = read.csv("age-length-rawVS.csv", header=TRUE)
age.length.rawVN = read.csv("age-length-rawVN.csv", header=TRUE)

age.length.rawVN_f <- subset(age.length.rawVN, Sex == "F")
age.length.rawVN_m <- subset(age.length.rawVN, Sex == "M")
age.length.rawVN_h <- subset(age.length.rawVN, Sex == "H")
age.length.rawVN_j <- subset(age.length.rawVN, Sex == "J")

age.length.rawVS_f <- subset(age.length.rawVS, Sex == "F")
age.length.rawVS_m <- subset(age.length.rawVS, Sex == "M")
age.length.rawVS_h <- subset(age.length.rawVS, Sex == "H")
age.length.rawVS_j <- subset(age.length.rawVS, Sex == "J")

#### Read results in from bootstrapped growth

resS <- read.csv("oto-S-all-1000results.csv", header = T)

#South            Linf           K         t0 Phi_vec
#CI_lower     55.47821 0.006688176 -11.414238 2.516683
#CI_upper    440.25647 0.118115194  -2.218605 2.737141
#max_density  68.59565 0.068655242  -4.718800 2.614907
#median       77.29374 0.066344910  -5.712392 2.617280

resN <- read.csv("oto-N-all-1000results.csv", header = T)

#North        Linf          K        t0   Phi_vec
#CI_lower    73.64609 0.05029839 -4.482484 2.704116
#CI_upper    98.90888 0.11711997 -1.631404 2.845078
#max_density 85.99679 0.07710472 -3.025265 2.765073
#median      85.63735 0.08062241 -3.027253 2.770901

resN_m <- read.csv("oto-N-mj-1000results.csv", header = T)

#             Linf          K         t0    Phi_vec
#CI_lower    38.96737 0.04638124 -2.6845132 2.746821
#CI_upper    93.43847 0.48536115  0.4127310 3.025234
#max_density 50.09587 0.24961268 -0.4503032 2.847746
#median      52.98438 0.26307843 -0.7814372 2.879455

resN_h <- read.csv("oto-N-herm-1000results.csv", header = T)
#             Linf            K          t0     Phi_vec
#CI_lower     31.54975 -0.005542514 -11.8398669 2.595118
#CI_upper    419.83599  0.171730413  -0.3357261 3.138356
#max_density  75.78489  0.033851470  -5.3905759 2.743273
#median       99.39827  0.050502819  -5.5812255 2.756027

resN_f <- read.csv("oto-N-female-1000results.csv", header = T)

#North female   Linf           K         t0     Phi_vec
#CI_lower     61.21894 -0.00463360 -26.072061 2.442285
#CI_upper    191.17809  0.17878492   3.156447 2.968343 
#max_density  88.14741  0.02992351  -3.658491 2.630163
#median       95.36352  0.05126335  -7.214562 2.678088

resS_m <- read.csv("oto-S-mj-1000results.csv", header = T)

#South male   Linf          K        t0         Phi_vec
#CI_lower     40.98460 -0.004164641 -9.5943805 2.501521
#CI_upper    139.13998  0.301085003  0.2707206 2.890752
#max_density  53.78748  0.139738763 -2.0701509 2.680395
#median       58.49556  0.138498330 -2.8824490 2.687267

resS_h <- read.csv("oto-S-herm-1000results.csv", header = T)

#South herm      Linf            K         t0    Phi_vec
#CI_lower       38.99467 -0.002958129 -15.077286 2.446454
#CI_upper    516.60980  0.121297456  -2.793909 3.072396
#max_density  85.05907  0.009607042 -10.861690 2.609255
#median      112.90811  0.030668667  -9.095498 2.643813

resS_f <- read.csv("oto-S-female-1000results.csv", header = T)

#South female Linf           K          t0      Phi_vec
#CI_lower     53.56263 0.006460748 -14.4253066 2.449335
#CI_upper    125.13152 0.177716936   0.9595037 2.848248
#max_density  65.71584 0.085637411  -2.8866989 2.637430
#median       68.63870 0.092537487  -4.6668980 2.650420

###################################
######## CALCULATE MAX DENS ETC
####################################
res = resN_h # repeat for all different results
names(res)

#Phi_prime = log10(K) + 2 log10(Linf)
res$Phi_vec <- log10(res$K) + 2 * log10(res$Linf)

############### 
# 'univariate_density' with output table for otoliths
##  Same  function as "univariate_density()", but also provides an output table, for your convenience -> see bootstrap
##############

tiff(filename = "univdens_N_herm.tiff",  height=1800, width=3200, res = 300) #Keep changing for diff. 

#THIS ONLY NEED TO RUN ONCE
CI = 95
use_hist = FALSE
nbreaks = 10
mar = c(1.5,2,2,0)
oma = c(1.5,0,0,0.5)
mgp = c(2,0.5,0)
tcl = -0.25
cex = 1

#PLOT DENSITY PLOT
op <- par(no.readonly = TRUE)
par(
  mfcol = c(1, ncol(res)),
  mar = mar, oma = oma,
  mgp = mgp, tcl = tcl, cex = cex
)

# Univariate plots -> Define expressions

VARS <- list(
  Linf = expression(italic(L)[infinity]),
  K = expression(italic(K)),
  t0 = expression(italic(t[0])),
  Phi_vec = expression(paste(phi,"'")))
  
Output_table <- as.data.frame(matrix(NA, nrow = 4, ncol =ncol(res))  ) 
names(Output_table) <-  (names(res)) 
rownames(Output_table) <- c("CI_lower","CI_upper", "max_density", "median")

# Univariate plots

for(i in seq(ncol(res))){
  x <- ks::kde(res[,i])
  
  h = hist(res[,i], plot=FALSE, breaks = nbreaks)
  
  xlim <- c(0, max(x$estimate))
  if(use_hist){
    xlim <- range(c(xlim, max(h$density)))
  }
  xlim <- xlim * c(0,1.1)
  
  plot(x$estimate, x$eval.points, t="n",
       xlim = xlim,
       xaxs = "i",
       ylab="", xlab="", col=1, lty=1
  )
  usr <- par()$usr
  
  CItxt <- paste0(round(100-CI), "%")
  inCI <- rle( x$estimate > x$cont[CItxt] )
  start.idx <- c(1, cumsum(inCI$lengths[-length(inCI$lengths)])+1)
  end.idx <- cumsum(inCI$lengths)
  limCI <- range(x$eval.points[start.idx[min(which(inCI$values))]:end.idx[max(which(inCI$values))]])
  
  in1 <- which(x$estimate > x$cont["99%"])
  mean1 <- mean(x$eval.points[in1])
  
  if(use_hist){
    rect(
      xleft = 0, ybottom = h$breaks[-length(h$breaks)],
      xright = h$density, ytop = h$breaks[-1],
      col = "grey90", border = "grey50"
    )
  }else{
    for(j in seq(inCI$lengths)){
      if(inCI$values[j]){
        polygon(
          y = c(x$eval.points[start.idx[j]:end.idx[j]], rev(x$eval.points[start.idx[j]:end.idx[j]])),
          x = c(x$estimate[start.idx[j]:end.idx[j]], rep(0, length(x$estimate[start.idx[j]:end.idx[j]]))),
          col = "grey90", #col = rgb(0,1,0,0.2),
          border = NA, lty = 3, lwd = 1
        )
      }
    }
  }
  
  lines(x$estimate, x$eval.points, lwd = 1, col = "grey50")
  
  # rug
  segments(x0 = 0, x1 = par()$cxy[1]*0.3, y0 = x$x, y1 = x$x, col=rgb(0,0,0,0.5), lwd=0.3)
  
  # range of CI
  abline(h = limCI, lty = 1, lwd=1, col = 1)
  text(y =c(limCI), x = mean(usr[1:2]),
       labels = paste(sprintf("%.2f", round(c(limCI),2))),
       pos = c(1,3), offset = 0.25, col = 1
  )
  abline(h = mean1, lty = 1, lwd=1, col = 1)
  text(y =  mean1, x = mean(usr[1:2]),
       labels = sprintf("%.2f", round(mean1,2)),
       pos = 3,
       offset = 0.25, col = 1
  )
  
  varlab <- VARS[[match(names(res)[i], names(VARS))]]
  mtext(varlab, line=0.25, side=3)
  
  box()
  
  
  Output_table[1,i] <- limCI[1]
  Output_table[2,i] <- limCI[2] 
  Output_table[3,i] <- mean1
  Output_table[4,i] <- median(res[,i])
  
}
mtext("Density", side = 1, line = 0, outer = TRUE)

dev.off()

Output_list <- list (Output_table) 

Output_table

###############################
#####
#Now calculate the 95 % confidence envelope AND plot 95% CI from matrix
#####  
##############################
t.seq.p = seq(0,50,0.1)

library(scales)

#plot.new()

#North ALL Swarm plot -> Need to initiate plot first

t.seq2 <- seq(0, 30, by = 0.1)
Linf.2 <-  round(resN$Linf[1])
K.2 <- round(resN$K[1])
t0.2 <- round(resN$t0[1])

L.seq.2 <- VBGF(t = t.seq2, list(Linf=Linf.2, K=K.2, C= 0, t0 = t0.2))  

plot(t.seq2, L.seq.2, t="l", col = "blue", xlim = c(0, 30), ylim = c(0, 120),
     xlab = "Age (y)", ylab = "Fork length (cm)")

#Then add swarm plots
for (w in  1:length(resN$K)) {
  
  Lt.p <- VBGF(list(Linf=resN$Linf[w],
                    K= resN$K[w],
                    t0= resN$t0[w]),
                    t=t.seq.p)
  
  lines(t.seq.p, Lt.p, col = alpha("grey", 0.1), lwd = 0.3 ) #grey transparent lines -> Only works on separate plot
  
}


# Now make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.
L.matN <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(resN$K)))

# Fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.matN))) {
  
  Lt <- VBGF(list(Linf=resN$Linf[w], 
                  K= resN$K[w], 
                  t0= resN$t0[w]),
             t=t.seq.p)
  
  
  L.matN[,w] <- round(Lt,2)
  
}

#Now for vector - North all
lo.ci.vecN <- as.numeric(apply(L.matN, 1, quantile, probs = 0.025))
up.ci.vecN <- as.numeric(apply(L.matN, 1, quantile, probs = 0.975))
med.vecN <-as.numeric(apply(L.matN, 1, quantile, probs = 0.5))

#Add line of median to the swarm plot
lines (t.seq.p, med.vecN, col = "red", lwd = 1)

###
###
#NOW FOR SOUTH 
#make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.
for (w in  1:length(resS$K)) {
  
  Lt.p <- VBGF(list(Linf=resS$Linf[w],
                    K= resS$K[w],
                    t0= resS$t0[w]),
               t=t.seq.p)
  
  lines(t.seq.p, Lt.p, col = alpha("grey", 0.1), lwd = 0.3 ) #grey transparent lines -> Only works on separate plot
  
}

L.matS <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(resS$K)))

# Fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.matS))) {
  
  Lt <- VBGF(list(Linf=resS$Linf[w], 
                  K= resS$K[w], 
                  t0= resS$t0[w]),
             t=t.seq.p)
  
  
  L.matS[,w] <- round(Lt,2)
  
}

#Now for vector - South all
lo.ci.vecS <- as.numeric(apply(L.matS, 1, quantile, probs = 0.025))
up.ci.vecS <- as.numeric(apply(L.matS, 1, quantile, probs = 0.975))
med.vecS <-as.numeric(apply(L.matS, 1, quantile, probs = 0.5))

#Add line of median to the swarm plot
lines (t.seq.p, med.vecS, col = "green", lwd = 1)

############################
#####
# Now calculate the 95 % confidence envelope AND plot 95% CI from matrix
#####  
###########################
t.seq.p = seq(0,50,0.1)

t.seq.p = as.array(t.seq.p)
med.vecS = as.array(med.vecS) 

lo.ci.vecS = as.array(lo.ci.vecS)
up.ci.vecS = as.array(up.ci.vecS)

med.vecN = as.array(med.vecN) 
lo.ci.vecN = as.array(lo.ci.vecN)
up.ci.vecN = as.array(up.ci.vecN)

library(scales)
view(up.ci.vecN)


#Add together
list_df = list(t.seq.p, med.vecS, lo.ci.vecS, up.ci.vecS, med.vecN, lo.ci.vecN, up.ci.vecN)
view(list_df)

areadf = list_df

#med.vecN, lo.ci.vecN, up.ci.vecN
#colnames(areadf)[1]<-"Age"
#colnames(areadf)[2]<-"med.S"
#colnames(areadf)[3]<-"lo.ci.S"
#colnames(areadf)[4]<-"up.ci.S"
#colnames(areadf)[5]<-"med.N"
#colnames(areadf)[6]<-"lo.ci.N"
#colnames(areadf)[7]<-"up.ci.N"

view(areadf)
write.csv(areadf, "otoNS_med_UL.csv")

#
#
# North female
# Swarm plots
for (w in  1:length(resN_f$K)) {
  
  Lt.p <- VBGF(list(Linf=resN_f$Linf[w],
                    K= resN_f$K[w],
                    t0= resN_f$t0[w]),
              t=t.seq.p)
  
  lines(t.seq.p, Lt.p, col = alpha("grey", 0.1), lwd = 0.3 ) #grey transparent lines -> Only works on separate plot
  
}

# Make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.

L.matNf <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(resN_f$K)))

# Fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.matNf))) {
  
  Lt <- VBGF(list(Linf=resN_f$Linf[w], 
                  K= resN_f$K[w], 
                  t0= resN_f$t0[w]),
             t=t.seq.p)
  
  
  L.matNf[,w] <- round(Lt,2)
  
}

# Now for vector - North female
lo.ci.vecNf <- as.numeric(apply(L.matNf, 1, quantile, probs = 0.025))
up.ci.vecNf <- as.numeric(apply(L.matNf, 1, quantile, probs = 0.975))
med.vecNf <-as.numeric(apply(L.matNf, 1, quantile, probs = 0.5))


#
#
# North male
#
# Make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.

L.matNm <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(resN_m$K)))

# Fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.matNm))) {
  
  Lt <- VBGF(list(Linf=resN_m$Linf[w], 
                  K= resN_m$K[w], 
                  t0= resN_m$t0[w]),
             t=t.seq.p)
  
  
  L.matNm[,w] <- round(Lt,2)
  
}

# Now for vector - North male
lo.ci.vecNm <- as.numeric(apply(L.matNm, 1, quantile, probs = 0.025))
up.ci.vecNm <- as.numeric(apply(L.matNm, 1, quantile, probs = 0.975))
med.vecNm <-as.numeric(apply(L.matNm, 1, quantile, probs = 0.5))

#
#
# North hermaphrodite
# Make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.

L.matNh <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(resN_h$K)))

# Fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.matNh))) {
  
  Lt <- VBGF(list(Linf=resN_h$Linf[w], 
                  K= resN_h$K[w], 
                  t0= resN_h$t0[w]),
             t=t.seq.p)
  
  
  L.matNh[,w] <- round(Lt,2)
  
}

# Now for vector - North hermaphrodite
lo.ci.vecNh <- as.numeric(apply(L.matNh, 1, quantile, probs = 0.025))
up.ci.vecNh <- as.numeric(apply(L.matNh, 1, quantile, probs = 0.975))
med.vecNh <-as.numeric(apply(L.matNh, 1, quantile, probs = 0.5))

# Add together
list_df = list(t.seq.p, med.vecNm, lo.ci.vecNm, up.ci.vecNm, med.vecNh, lo.ci.vecNh, up.ci.vecNh,med.vecNf, lo.ci.vecNf, up.ci.vecNf)
areadf = list_df
write.csv(areadf, "otoN_mfh_med_UL.csv")

#
#
#
# South female
# Make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.

L.matSf <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(resS_f$K)))

# Fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.matSf))) {
  
  Lt <- VBGF(list(Linf=resS_f$Linf[w], 
                  K= resS_f$K[w], 
                  t0= resS_f$t0[w]),
             t=t.seq.p)
  
  
  L.matSf[,w] <- round(Lt,2)
  
}

# Now for vector - South female
lo.ci.vecSf <- as.numeric(apply(L.matSf, 1, quantile, probs = 0.025))
up.ci.vecSf <- as.numeric(apply(L.matSf, 1, quantile, probs = 0.975))
med.vecSf <-as.numeric(apply(L.matSf, 1, quantile, probs = 0.5))


#
#
# South male
# Make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.

L.matSm <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(resS_m$K)))

# Fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.matSm))) {
  
  Lt <- VBGF(list(Linf=resS_m$Linf[w], 
                  K= resS_m$K[w], 
                  t0= resS_m$t0[w]),
             t=t.seq.p)
  
  L.matSm[,w] <- round(Lt,2)
  
}

# Now for vector - South male
lo.ci.vecSm <- as.numeric(apply(L.matSm, 1, quantile, probs = 0.025))
up.ci.vecSm <- as.numeric(apply(L.matSm, 1, quantile, probs = 0.975))
med.vecSm <-as.numeric(apply(L.matSm, 1, quantile, probs = 0.5))

#
#
# South hermaphrodite
# Make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.

L.matSh <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(resS_h$K)))

# Fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.matSh))) {
  
  Lt <- VBGF(list(Linf=resS_h$Linf[w], 
                  K= resS_h$K[w], 
                  t0= resS_h$t0[w]),
             t=t.seq.p)

  L.matSh[,w] <- round(Lt,2)
  
}

# Now for vector - South hermaphrodite
lo.ci.vecSh <- as.numeric(apply(L.matSh, 1, quantile, probs = 0.025))
up.ci.vecSh <- as.numeric(apply(L.matSh, 1, quantile, probs = 0.975))
med.vecSh <-as.numeric(apply(L.matSh, 1, quantile, probs = 0.5))

# Add together
list_df = list(t.seq.p, med.vecSm, lo.ci.vecSm, up.ci.vecSm, med.vecSh, lo.ci.vecSh, up.ci.vecSh,med.vecSf, lo.ci.vecSf, up.ci.vecSf)
areadf = list_df
write.csv(areadf, "otoS_mfh_med_UL.csv")


######################################
#######
#LFA RESULTS
########
######################################
## RESULTS OBTAINED FROM: 
## BOOTSTRAPPED USING SCRIPT: 
## "#SteenBras71_GAboot_runs_only_M3_MA5_Linfmin30_50ALL_Sizes.R"
## READ IN & ADDED TOGETHER & SCREENED USING SCRIPT: 
## "steenbrasS_v31_rs_reading_in_and_analyzing_comparing_wFINAL_plots.R"

path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

#path2 = path +'/' + "LFA_S_stock_New_analyses_Sys_time_results"

res <- read.table("res_MA3_Linfmin50_Sys.TimeOK.csv", header = T) 
#Need this for density plot, Kimura plot etc. But only use one for plotting with otoliths and T&R

L.mat <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(res$K)))

# fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.mat))) {
  
  Lt <- VBGF(list(Linf=res$Linf[w], 
                  K= res$K[w], 
                  t0= -1,
                  C=res$C[w],
                  ts=res$ts[w]),
                  t=t.seq.p)
  
  L.mat[,w] <- round(Lt,2)
  
}

# now for vector - South LFA
lo.ci.vec <- as.numeric(apply(L.mat, 1, quantile, probs = 0.025))
up.ci.vec <- as.numeric(apply(L.mat, 1, quantile, probs = 0.975))
med.vec <-as.numeric(apply(L.mat, 1, quantile, probs = 0.5))

view(lo.ci.vec)
 
#
#
#
## PLOT ALL TOGETHER
#
#
plot(t.seq.p, med.vecN, t="l", col = "red", lwd = 3, xlim = c(0,50), ylim = c(0,95), 
     xlab = "Age (y)", ylab = "Fork length (cm)") #North median

lines(t.seq.p, med.vecN, t="l", col = "red", lwd = 3)
     
points(age.length.rawVN$length_cm ~ age.length.rawVN$age_yr, col = "red") #raw north
lines(t.seq.p, lo.ci.vecN, lty = 3, col = "red", lwd = 1.5) #lower CI N
lines(t.seq.p, up.ci.vecN, lty = 3, col = "red", lwd = 1.5) #upper CI N

#col = alpha("grey", 0.1), lwd = 0.3 
plot(t.seq.p, med.vecS, t="l", col = "navy", lwd = 3, xlim = c(0,50), ylim = c(0,95), 
     xlab = "Age (y)", ylab = "Fork length (cm)")

points(age.length.rawVS$length_cm ~ age.length.rawVS$age_yr, col = "navy") #raw south
lines(t.seq.p, lo.ci.vecS, lty = 3, col = "navy", lwd = 1.5) #lower CI south
lines(t.seq.p, up.ci.vecS, lty = 3, col = "navy", lwd = 1.5) #upper CI south
lines(t.seq.p, med.vecS, t="l", col = "navy", lwd = 3) #median south

points(age.length.rawVN_f$length_cm ~ age.length.rawVN_f$age_yr, col = "orange", pch =2) #north female
lines(t.seq.p, lo.ci.vecNf, lty = 3, col = "orange", lwd = 1.5)
lines(t.seq.p, up.ci.vecNf, lty = 3, col = "orange", lwd = 1.5)
lines(t.seq.p, med.vecNf, t="l", col = "orange", lwd = 3)

points(age.length.rawVN_m$length_cm ~ age.length.rawVN_m$age_yr, col = "blue", pch = 2) #north male
lines(t.seq.p, lo.ci.vecNm, lty = 3, col = "blue", lwd = 1.5)
lines(t.seq.p, up.ci.vecNm, lty = 3, col = "blue", lwd = 1.5)
lines(t.seq.p, med.vecNm, t="l", col = "blue", lwd = 3)

points(age.length.rawVN_h$length_cm ~ age.length.rawVN_h$age_yr, col = "green", pch = 2) #north herm
lines(t.seq.p, lo.ci.vecNh, lty = 3, col = "green", lwd = 1.5)
lines(t.seq.p, up.ci.vecNh, lty = 3, col = "green", lwd = 1.5)
lines(t.seq.p, med.vecNh, t="l", col = "green", lwd = 3)

points(age.length.rawVN_j$length_cm ~ age.length.rawVN_j$age_yr, col = "black", pch = 2)

points(age.length.rawN$length_cm ~ age.length.rawN$age_yr, col = "green", pch = 2)
lines(t.seq.p, lo.ci.vecNHK, lty = 3, col = "green", lwd = 1.5)
lines(t.seq.p, up.ci.vecNHK, lty = 3, col = "green", lwd = 1.5)
lines(t.seq.p, med.vecNHK, t="l", col = "green", lwd = 3)

points(age.length.rawS$length_cm ~ age.length.rawS$age_yr, col = "blue", pch = 2)
lines(t.seq.p, lo.ci.vecSHK, lty = 3, col = "blue", lwd = 1.5)
lines(t.seq.p, up.ci.vecSHK, lty = 3, col = "blue", lwd = 1.5)
lines(t.seq.p, med.vecSHK, t="l", col = "blue", lwd = 3)


##
##
##
#EASIEST IS TO READ IT IN FROM EXCEL ==> LAA CI of tagging analysis
path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

laa95 = read.csv("laa95.csv", header = TRUE)
names(laa95)

##-----------------------------------------------
###
#Tagging median South
t <- seq(0,50,0.1)
L2 <- VBGF(list(Linf=49.7, K=0.185,  t0=-1, C=0, ts=0), t=t) 
#to = -1

p5la = xyplot(Lower.LAA.T.R ~ Age, data = laa95, type='l', pch = 5, lty = 2, lwd = 2, col = "orange") #T-R
p5ua = xyplot(Upper.LAA.T.R ~ Age, data = laa95, type='l', pch = 5, lty=2, lwd = 2, col = "orange")

p5ma = xyplot(med.vec  ~ t.seq.p, type='l', pch = 5, lty = 2, lwd = 2, col = "green") #LFA
p5ub = xyplot(up.ci.vec  ~ t.seq.p, type='l', pch = 5, lty = 2, lwd = 2, col = "green") #LFA
p5lb = xyplot(lo.ci.vec ~ t.seq.p, type='l', pch = 5, lty = 2, lwd = 2, col = "green") #Calculated above from LFA


plot(t.seq.p, med.vecN, t="l", col = "red", lwd = 3, xlim = c(0,50), ylim = c(0,95), 
     xlab = "Age (y)", ylab = "Fork length (cm)") #North median
points(age.length.rawVN$length_cm ~ age.length.rawVN$age_yr, col = "red") #raw north
lines(t.seq.p, lo.ci.vecN, lty = 3, col = "red", lwd = 1.5) #lower CI N
lines(t.seq.p, up.ci.vecN, lty = 3, col = "red", lwd = 1.5) #upper CI N
#col = alpha("grey", 0.1), lwd = 0.3 

#####
#create scatterplots of raw data -> library(latticeExtra)
#Red is north, navy is south, orange female, blue male, green herm. 
f_key <- list(x = .95, y = 0.01, corner = c(1, 0),
              text = list(c("North male fish", 
                            "North female fish",
                            "North hermaphrodite fish",
                            "North median", 
                            "CI = 95%",
                            "North max. dens.")), 
              lines = list(type = c("p", "p", "p", "l", "l", "l"), 
              col = c("blue", "orange", "black","red", "red", "black"),
              pch = c(1, 2, 3, 5, 5, 5),
              lty = c(1, 1, 2, 1, 2, 1),
              lwd = c(5, 5, 5, 5, 3, 4)))


# North maximum density
#North        Linf          K        t0
#max_density 85.99679 0.07710472 -3.025265
LNmd <- VBGF(list(Linf= 85.99679, K=0.07710472,  t0=-3.025265, C=0, ts=0), t=t.seq.p) 

# South            Linf           K         t0
#max_density  68.59565 0.068655242  -4.718800
LSmd <- VBGF(list(Linf=68.59565, K=0.068655242,  t0=-4.718800, C=0, ts=0), t=t.seq.p) 


p0 = xyplot(med.vecN ~ t.seq.p, xlab= "Age (years)", ylab = "Fork length (cm)",
          col.lab = "black", cex.lab = 20, panel.axis(text.fontface="bold"), 
          type=c('l'), lwd = 5,
          col = "red", pch=4, key=f_key,
          scales=list(x=list(limits=c(0,50), at = seq(0,50, by=2)), 
                      y= list(limits=c(0,95), at = seq(0, 95, by = 5) )))

pNVbl <- xyplot(lo.ci.vecN ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "red") #North 95% CI this study
pNVbu <- xyplot(up.ci.vecN ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "red") #North 95% CI this study
pNmd <- xyplot(LNmd ~ t.seq.p, type='l', lty = 1, lwd = 2, col = "grey") #North max dens.

NV_f = xyplot(age.length.rawVN_f$length_cm ~ age.length.rawVN_f$age_yr, col = "orange", pch =2) #north female
NV_m = xyplot(age.length.rawVN_m$length_cm ~ age.length.rawVN_m$age_yr, col = "blue", pch =1) #north female
NV_h = xyplot(age.length.rawVN_h$length_cm ~ age.length.rawVN_h$age_yr, col = "black", pch =3) #north female

NV = xyplot(age.length.rawVN$length_cm ~ age.length.rawVN$age_yr, col = "red", pch =1) #north female


tiff("VBGF_north_all.tiff", height=1800, width=3200, res = 300) #FOR POSTER - redone
opar <- par(mfrow = c(1,1), mar = c(3,3,2,1), oma = c(2,0,0,0))
p0 + as.layer(pNVbu) + as.layer(pNVbl) + as.layer(NV_f) + as.layer(NV_m) + as.layer(NV_h) + as.layer(pNmd)
dev.off()


##South
##
f1_key <- list(x = .95, y = 0.01, corner = c(1, 0),
               text = list(c("South male fish", 
                             "South female fish",
                             "South hermaphrodite fish",
                             "South median", 
                             "CI = 95%",
                             "South max. dens.")), 
               lines = list(type = c("p", "p", "p", "l", "l", "l"), 
                            col = c("blue", "orange", "black","navy", "navy", "black"),
                            pch = c(1, 2, 3, 5, 5, 5),
                            lty = c(1, 1, 2, 1, 2, 1),
                            lwd = c(5, 5, 5, 5, 3, 4)))

#Using f1_key
p1=xyplot(med.vecS ~ t.seq.p, xlab= "Age (years)", ylab = "Fork length (cm)",
          col.lab = "black", cex.lab = 20, panel.axis(text.fontface="bold"), 
          type=c('l'), lwd = 5,
          col = "navy", pch=4, key=f1_key,
          scales=list(x=list(limits=c(0,50), at = seq(0,50, by=2)), 
                      y= list(limits=c(0,95), at = seq(0, 95, by = 5) )))
p1

pSVbl <- xyplot(lo.ci.vecS ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "navy") #South 95% CI this study
pSVbu <- xyplot(up.ci.vecS ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "navy") #South 95% CI this study
pSmd <- xyplot(LSmd ~ t.seq.p, type='l', lty = 1, lwd = 4, col = "black") #South max dens.

SV_f = xyplot(age.length.rawVS_f$length_cm ~ age.length.rawVS_f$age_yr, col = "orange", pch =2) #north female
SV_m = xyplot(age.length.rawVS_m$length_cm ~ age.length.rawVS_m$age_yr, col = "blue", pch =1) #north female
SV_h = xyplot(age.length.rawVS_h$length_cm ~ age.length.rawVS_h$age_yr, col = "black", pch =3) #north female

SV = xyplot(age.length.rawVS$length_cm ~ age.length.rawVS$age_yr, col = "navy", pch =1) #north female

tiff("VBGF_south_all.tiff", height=1800, width=3200, res = 300) #FOR POSTER - redone
opar <- par(mfrow = c(1,1), mar = c(3,3,2,1), oma = c(2,0,0,0))
p1 + as.layer(pSVbu) + as.layer(pSVbl) + as.layer(SV_f) + as.layer(SV_m) + as.layer(SV_h) + as.layer(pSmd)
dev.off()

A = p0 + as.layer(pNVbu) + as.layer(pNVbl) + as.layer(NV_f) + as.layer(NV_m) + as.layer(NV_h) + as.layer(pNmd)
B = p1 + as.layer(pSVbu) + as.layer(pSVbl) + as.layer(SV_f) + as.layer(SV_m) + as.layer(SV_h) + as.layer(pSmd)


tiff("North&South_all.tiff", height=8, width=6, units = "in", res = 600, compression = "lzw")
par(mar=c(4.2, 3.8, 0.2, 0.2))#bottom, left, top, right
print(A, split=c(1,1,1,2), more=TRUE)
print(B, split=c(1,2,1,2), more=TRUE)
dev.off()

#
#
#
#North ALL
for (w in  1:length(resN$K)) {
  
  Lt.p <- VBGF(list(Linf=resN$Linf[w],
                    K= resN$K[w],
                    t0= resN$t0[w]),
               t=t.seq.p)
  
  lines(t.seq.p, Lt.p, col = alpha("grey", 0.1), lwd = 0.3 ) #grey transparent lines -> Only works on separate plot
  
}

#p0 + as.layer(pNV) + as.layer(pNVbu) + as.layer(pNVbl) + as.layer(psn) #doesn't work


########
########
#All overlay
########
########

fb_key <- list(x = .95, y = 0.01, corner = c(1, 0),
              text = list(c("Max dens LFA S (t0 = -1)", 
                            "Median LFA S (t0 = -1)",
                            "95% CI LFA S",
                            "Median T&R S (t0 = -1)", 
                            "95% CI T&R S",
                            "Bootstrapped otolith data N median", 
                            "Bootstrapped otolith data N 95% CI",
                            "Bootstrapped otolith data S median",
                            "Bootstrapped otolith data S 95% CI")), 
                            lines = list(type = c("l", "l", "l", "l", "l", "l", "l", "l", "l"), 
                           col = c("#009E73", "green", "green","orange", "orange","red", "red", "navy", "navy"),
                           lty = c(1, 1, 2, 1, 2, 1, 2, 1, 2),
                           pch = c(5, 5, 5, 5, 5, 5, 5, 5, 5), 
                           lwd = c(5, 5, 2, 5, 2, 5, 2, 5, 2)))

pSV <- xyplot(med.vecS ~ t.seq.p, type='l', lwd = 5, col = "navy") #Median South this study

#max density LFA MA3, Linf30
Lt1 <- VBGF(list(Linf=75.57374, K=0.11432142,  t0=-1, C=0.5638764, ts=0.51289124), t=t.seq.p) 

p0b=xyplot(Lt1 ~ t.seq.p, xlab= "Age (years)", ylab = "Fork length (cm)",
           col.lab = "black", cex.lab = 20, panel.axis(text.fontface="bold"), #, text.cex = 6), 
           type=c('l'), lwd = 5,
           col = "#009E73", pch=4, key=fb_key,
           scales=list(x=list(limits=c(0,50), at = seq(0,50, by=2)), 
                       y= list(limits=c(0,95), at = seq(0, 95, by = 5) )))

p0b + as.layer(pSV) + as.layer(pNV) + as.layer(pNVbl) + as.layer(pNVbu) + as.layer(pSVbl) + as.layer(pSVbu)

#Tagging median South
t <- seq(0,50,0.1)
L2 <- VBGF(list(Linf=49.7, K=0.185,  t0=-1, C=0, ts=0), t=t) 

p2 <- xyplot(L2 ~ t, type='l', lwd = 5, col = "orange") #T&R S median
p5la = xyplot(Lower.LAA.T.R ~ Age, data = laa95, type='l', pch = 5, lty = 2, lwd = 2, col = "orange") #T-R
p5ua = xyplot(Upper.LAA.T.R ~ Age, data = laa95, type='l', pch = 5, lty=2, lwd = 2, col = "orange")
#to = -1

p0b + as.layer(p2) + as.layer(p5la) + as.layer(p5ua) + as.layer(pSV) + as.layer(pNV) + as.layer(pNVbl) + as.layer(pNVbu) + as.layer(pSVbl) + as.layer(pSVbu)

p5ma = xyplot(med.vec  ~ t.seq.p, type='l', lwd = 5, col = "green") #LFA med
p5ub = xyplot(up.ci.vec  ~ t.seq.p, type='l', pch = 5, lty = 2, lwd = 2, col = "green") #LFA
p5lb = xyplot(lo.ci.vec ~ t.seq.p, type='l', pch = 5, lty = 2, lwd = 2, col = "green") #Calculated above from LFA

tiff("VBGF_overlay_all.tiff", height=1800, width=3200, res = 300) 
opar <- par(mfrow = c(1,1), mar = c(3,3,2,1), oma = c(2,0,0,0))
p0b + as.layer(p5ma) + as.layer(p5ub) + as.layer(p5lb) + as.layer(p2) + as.layer(p5la) + as.layer(p5ua) + as.layer(pSV) + as.layer(pNV) + as.layer(pNVbl) + as.layer(pNVbu) + as.layer(pSVbl) + as.layer(pSVbu)
dev.off()


##
#
#North & South Otoliths Overlay
f1b_key <- list(x = .95, y = 0.01, corner = c(1, 0),
               text = list(c("Raw otolith data N",
                             "Bootstrapped LAA Median N",
                             "95% CI N",
                             "Bootstrapped LAA Max. Dens. N",
                             "Raw otolith data S",
                             "Bootstrapped LAA Median S",
                             "95% CI S",
                             "Bootstrapped LAA Max. Dens. S")), 
                            lines = list(type = c("p", "l", "l", "l", "p", "l", "l", "l"), 
                            col = c("red", "red", "red", "grey", "navy", "navy", "navy", "black"),
                            pch = c(1, 5, 5, 5, 1, 5, 5, 5), 
                            lty = c(1, 1, 2, 1, 1, 1, 2, 1),
                            lwd = c(5, 5, 3, 3, 5, 5, 3, 4)))

#Using f1b_key
p1b=xyplot(med.vecN ~ t.seq.p, xlab= "Age (years)", ylab = "Fork length (cm)",
                  col.lab = "black", cex.lab = 20, panel.axis(text.fontface="bold"), 
                  type=c('l'), lwd = 5,
                  col = "red", pch=4, key=f1b_key,
                  scales=list(x=list(limits=c(0,50), at = seq(0,50, by=2)), 
                              y= list(limits=c(0,95), at = seq(0, 95, by = 5) )))


tiff("VBGF_North&south_overlay.tiff", height=1800, width=3200, res = 300) 
opar <- par(mfrow = c(1,1), mar = c(3,3,2,1), oma = c(2,0,0,0))
p1b + as.layer(pNVbu) + as.layer(pNVbl) + as.layer(pNmd) + as.layer(pSVbl) + as.layer(pSVbu)+ as.layer(pSV) + as.layer(pSmd) + as.layer(SV) + as.layer(NV)
dev.off()

#
#
###
##
fb2_key <- list(x = .95, y = 0.01, corner = c(1, 0),
               text = list(c("Max dens LFA S (t0 = -1)", 
                             "Median LFA S (t0 = -1)",
                             "95% CI LFA S",
                             "Median T&R S (t0 = -1)", 
                             "95% CI T&R S",
                             "Bootstrapped otolith data S median",
                             "Bootstrapped otolith data S 95% CI",
                             "Raw otolith data S female",
                             "Bootstrapped otolith data S female median", 
                             "Bootstrapped otolith data S female 95% CI", 
                             "Raw otolith data S male",
                             "Bootstrapped otolith data S male median", 
                             "Bootstrapped otolith data S male 95% CI", 
                             "Raw otolith data S hermaphrodite", 
                             "Bootstrapped otolith data S herm", 
                             "Bootstrapped otolith data S herm 95% CI")), 
                            lines = list(type = c("l", "l", "l", "l", "l", "l", "l", "p", "l", "l", "p", "l", "l", "p", "l", "l"), 
                            col = c("#009E73", "green", "green","orange", "orange", "black", "black", "red", "red", "red", "blue", "blue", "blue", "black", "grey", "grey"),
                            lty = c(1, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2),
                            pch = c(5, 5, 5, 5, 5, 5, 5, 2, 5, 5, 2, 5, 5, 3, 5, 5), 
                            lwd = c(5, 5, 2, 5, 2, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2)))

p1c=xyplot(Lt ~ t, xlab= "Age (years)", ylab = "Fork length (cm)",
           col.lab = "black", cex.lab = 20, panel.axis(text.fontface="bold"), #, text.cex = 6), 
           type=c('l'), lwd = 5,
           col = "#009E73", pch=4, key=fb2_key,
           scales=list(x=list(limits=c(0,50), at = seq(0,50, by=2)), 
                       y= list(limits=c(0,95), at = seq(0, 95, by = 5) )))
p1c
####

fb2d_key <- list(x = .95, y = 0.01, corner = c(1, 0),
                text = list(c("Max dens LFA S (t0 = -1)", 
                              "Median LFA S (t0 = -1)",
                              "95% CI LFA S",
                              "Bootstrapped otolith data N median",
                              "Bootstrapped otolith data N 95% CI",
                              "Raw otolith data N female",
                              "Bootstrapped otolith data N female median", 
                              "Bootstrapped otolith data N female 95% CI", 
                              "Raw otolith data N male",
                              "Bootstrapped otolith data N male median", 
                              "Bootstrapped otolith data N male 95% CI", 
                              "Raw otolith data N hermaphrodite", 
                              "Bootstrapped otolith data N herm", 
                              "Bootstrapped otolith data N herm 95% CI")), 
                              lines = list(type = c("l", "l", "l", "l", "l", "p", "l", "l", "p", "l", "l", "p", "l", "l"), 
                             col = c("#009E73", "green", "green","black", "black", "red", "red", "red", "blue", "blue", "blue", "black", "grey", "grey"),
                             lty = c(1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2),
                             pch = c(5, 5, 5, 5, 5, 2, 5, 5, 2, 5, 5, 3, 5, 5), 
                             lwd = c(5, 5, 2, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2)))

p1d = xyplot(Lt ~ t, xlab= "Age (years)", ylab = "Fork length (cm)",
           col.lab = "black", cex.lab = 20, panel.axis(text.fontface="bold"), #, text.cex = 6), 
           type=c('l'), lwd = 5,
           col = "#009E73", pch=4, key=fb2d_key,
           scales=list(x=list(limits=c(0,50), at = seq(0,50, by=2)), 
                       y= list(limits=c(0,95), at = seq(0, 95, by = 5) )))
p1d 

####
####
# ALL MEDIANS AND 95% CI
## create second scatterplot and add it to first plot
p2 <- xyplot(L2 ~ t, type='l', lwd = 5, col = "orange") #T&R S median
p2b <- xyplot(L2b ~ t, type='l', lwd = 5, col = "green") #LFA S median

pN <- xyplot(LN ~ t, type='l', lty = 2, lwd = 5, col = "red") #Median North (H&K)
pS <- xyplot(LS ~ t, type='l', lty = 2, lwd = 5, col = "navy") #Median South (H&K)

pNb <- xyplot(LN ~ t, type='l', lty = 1, lwd = 5, col = "blue") #Median North (H&K)
pSb <- xyplot(LS ~ t, type='l', lty = 1, lwd = 5, col = "black") #Median South (H&K)

pSc <- xyplot(LS ~ t, type='l', lty = 1, lwd = 5, col = "blue") #Median South (H&K)

pNd <- xyplot(LN ~ t, type='l', lty = 2, lwd = 5, col = "blue") #Median North (H&K)
pSd <- xyplot(LS ~ t, type='l', lty = 2, lwd = 5, col = "black") #Median South (H&K)

pScl <- xyplot(lo.ci.vecSHK ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "blue") #South (H&K) L 95%
pScu <- xyplot(up.ci.vecSHK ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "blue") #South (H&K) U 95%
pNcl <- xyplot(lo.ci.vecNHK ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "green") #South (H&K) L 95%
pNcu <- xyplot(up.ci.vecNHK ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "green") #South (H&K) U 95%

pNV <- xyplot(med.vecN ~ t.seq.p, type='l', lwd = 5, col = "red") #Median North this study
pNV2 <- xyplot(med.vecN ~ t.seq.p, type='l', lwd = 5, col = "black") #Median North this study

pSV <- xyplot(med.vecS ~ t.seq.p, type='l', lwd = 5, col = "navy") #Median South this study

pNVbl <- xyplot(lo.ci.vecN ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "red") #North 95% CI this study
pNVbl2 <- xyplot(lo.ci.vecN ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "black") #North 95% CI this study

pSVbl <- xyplot(lo.ci.vecS ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "black") #South 95% CI this study

pNVbu <- xyplot(up.ci.vecN ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "red") #North 95% CI this study
pNVbu2 <- xyplot(up.ci.vecN ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "black") #North 95% CI this study

pSVbu <- xyplot(up.ci.vecS ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "black") #South 95% CI this study

pSVf <- xyplot(med.vecSf ~ t.seq.p, type='l', lwd = 5, col = "red") #Median South female this study
pSVfl <- xyplot(lo.ci.vecSf ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "red") #South 95% CI female this study
pSVfu <- xyplot(up.ci.vecSf ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "red") #South 95% CI female this study

pSVm <- xyplot(med.vecSm ~ t.seq.p, type='l', lwd = 5, col = "blue") #Median South male this study
pSVml <- xyplot(lo.ci.vecSm ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "blue") #South 95% CI male this study
pSVmu <- xyplot(up.ci.vecSm ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "blue") #South 95% CI male this study

pSVh <- xyplot(med.vecSh ~ t.seq.p, type='l', lwd = 5, col = "grey") #Median South herm this study
pSVhl <- xyplot(lo.ci.vecSh ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "grey") #South 95% CI herm this study
pSVhu <- xyplot(up.ci.vecSh ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "grey") #South 95% CI herm this study

pNVf <- xyplot(med.vecNf ~ t.seq.p, type='l', lwd = 5, col = "red") #Median North female this study
pNVfl <- xyplot(lo.ci.vecNf ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "red") #South 95% CI female this study
pNVfu <- xyplot(up.ci.vecNf ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "red") #South 95% CI female this study

pNVm <- xyplot(med.vecNm ~ t.seq.p, type='l', lwd = 5, col = "blue") #Median North male this study
pNVml <- xyplot(lo.ci.vecNm ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "blue") #South 95% CI male this study
pNVmu <- xyplot(up.ci.vecNm ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "blue") #South 95% CI male this study

pNVh <- xyplot(med.vecNh ~ t.seq.p, type='l', lwd = 5, col = "grey") #Median North herm this study
pNVhl <- xyplot(lo.ci.vecNh ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "grey") #South 95% CI herm this study
pNVhu <- xyplot(up.ci.vecNh ~ t.seq.p, type='l', lty = 2, lwd = 2, col = "grey") #South 95% CI herm this study

#
#
#
### RAW DATA
#Now read in age-length data Holtzhausen & Kirchner
age.length.raw = read.csv("age-length-raw.csv", header=TRUE)
names(age.length.raw)

#p3 <- xyplot(H.K.South ~ xH.K, data = age.length.raw, type='p', pch = 1, col = "navy") #H&K South raw data
#p4 <- xyplot(yH.K.North ~ xH.K.North, data = age.length.raw, type='p', pch = 1, col = "red") #H&K North

#p32 <- xyplot(H.K.South ~ xH.K, data = age.length.raw, type='p', pch = 1, col = "black") #H&K South raw data
#p42 <- xyplot(yH.K.North ~ xH.K.North, data = age.length.raw, type='p', pch = 1, col = "blue") #H&K North

#p33 <- xyplot(H.K.South ~ xH.K, data = age.length.raw, type='p', pch = 2, col = "blue") #H&K South raw data
#p43 <- xyplot(yH.K.North ~ xH.K.North, data = age.length.raw, type='p', pch = 2, col = "green") #H&K North

p3b <- xyplot(length_cm ~ age_yr, data = age.length.rawVN, type='p', pch = 2, col = "red") #Veii raw North
p4b <- xyplot(length_cm ~ age_yr, data = age.length.rawVS, type='p', pch = 2, col = "navy") #Veii raw South

p3b1 <- xyplot(length_cm ~ age_yr, data = age.length.rawVN, type='p', pch = 1, col = "red") #Veii raw North
p4b1 <- xyplot(length_cm ~ age_yr, data = age.length.rawVS, type='p', pch = 1, col = "black") #Veii raw South

p4bf <-  xyplot(length_cm ~ age_yr, data = age.length.rawVS_f, type='p', pch = 2, col = "red") #Veii South female
p4bm <-  xyplot(length_cm ~ age_yr, data = age.length.rawVS_m, type='p', pch = 2, col = "blue") #Veii South male
p4bh <-  xyplot(length_cm ~ age_yr, data = age.length.rawVS_h, type='p', pch = 3, col = "black") #Veii South hermaphrodite

p5bf <-  xyplot(length_cm ~ age_yr, data = age.length.rawVN_f, type='p', pch = 2, col = "red") #Veii North female
p5bm <-  xyplot(length_cm ~ age_yr, data = age.length.rawVN_m, type='p', pch = 2, col = "blue") #Veii North male
p5bh <-  xyplot(length_cm ~ age_yr, data = age.length.rawVN_h, type='p', pch = 3, col = "black") #Veii North hermaphrodite

points(age.length.rawVN_h$length_cm ~ age.length.rawVN_h$age_yr, col = "green", pch = 2)
lines(t.seq.p, lo.ci.vecNh, lty = 3, col = "green", lwd = 1.5)
lines(t.seq.p, up.ci.vecNh, lty = 3, col = "green", lwd = 1.5)
lines(t.seq.p, med.vecNh, t="l", col = "green", lwd = 3)


##
##
##Create TIFF file 
#
#
tiff("VBGF_overlays_Veii1.tiff", height=1800, width=3200, res = 300)
opar <- par(mfrow = c(1,1), mar = c(3,3,2,1), oma = c(2,0,0,0))
p0b + #using fb_key - need black circled for otolith data S
  as.layer(p2) + as.layer(p3b) + as.layer(p4b1) + as.layer(p2b) + 
  as.layer(p5la) + as.layer(p5ua) + as.layer(p5lb) + as.layer(p5ub) + 
  as.layer(pNV) + as.layer(pSVb) + 
  as.layer(pNVbl) + as.layer(pSVbl) + as.layer(pNVbu) + as.layer(pSVbu)
dev.off()

#
#
tiff("VBGF_overlays_withVeii2.tiff", height=1800, width=3200, res = 300)
opar <- par(mfrow = c(1,1), mar = c(3,3,2,1), oma = c(2,0,0,0))
p1 + # use key: f1_key #max density in red
  as.layer(p2b) + as.layer(p5lb) + as.layer(p5ub) + #LFA median upper and lower
  as.layer(p32) + as.layer(p42) + as.layer(p3b) + as.layer(p4b) + #raw data
  as.layer(pNd) +as.layer(pSd) #north & south H$K fitted
dev.off()

#
tiff("VBGF_overlays_withVeii3.tiff", height=1800, width=3200, res = 300) # f2 key
p1bf2 + as.layer(p2b) + as.layer(p3) + as.layer(p3b) + as.layer(p4) + as.layer(p4b) + as.layer(p2b) + as.layer(pN) +as.layer(pS) + as.layer(pNV) +as.layer(pSV)
dev.off()
#
p1bf2 #DONE

tiff("VBGF_overlays_withVeii4.tiff", height=1800, width=3200, res = 300)
p1c + as.layer(p2) + #LFA max dens, T&R median
  as.layer(p2b) + as.layer(p5ub) + as.layer(p5lb) + #LFA median, upper and lower
  as.layer(p5la) + as.layer(p5ua) + #tag-rec lower & upper 
  as.layer(pSV) + as.layer(pSVbu) + as.layer(pSVbl) + #South all median, upper and lower
  as.layer(pSVf) + as.layer(pSVfl) + as.layer(pSVfu) + #median, upper and lower female
  as.layer(p4bm) + as.layer(p4bf) + as.layer(p4bh) + #raw data of male, female, herm South
  as.layer(pSVm) + as.layer(pSVml) + as.layer(pSVmu) + #median upper and lower male
  as.layer(pSVh) + as.layer(pSVhl) + as.layer(pSVhu) #median upper and lower herm
dev.off()

#
#
tiff("VBGF_overlays_withVeii5.tiff", height=1800, width=3200, res = 300) #f1b key DONE
p1b + as.layer(pNV) + as.layer(pSVb) + as.layer(pSc) + as.layer(p3b1) + 
  as.layer(p43) + as.layer(p4b1) + as.layer(p33) + 
  as.layer(pSVbl) + as.layer(pSVbu) + as.layer(pNVbl) + as.layer(pNVbu)+
  as.layer(pNcl) + as.layer(pNcu) + as.layer(pScl) + as.layer(pScu)  
dev.off()


tiff("VBGF_overlays_withVeii6.tiff", height=1800, width=3200, res = 300) #DONE
p1d  + #LFA max dens
  as.layer(p2b) + as.layer(p5ub) + as.layer(p5lb) + #LFA median, upper and lower
  as.layer(pNV2) + as.layer(pNVbu2) + as.layer(pNVbl2) + #North all median, upper and lower
  as.layer(pNVf) + as.layer(pNVfl) + as.layer(pNVfu) + #median, upper and lower female
  as.layer(p5bm) + as.layer(p5bf) + as.layer(p5bh) + #raw data of male, female, herm North
  as.layer(pNVm) + as.layer(pNVml) + as.layer(pNVmu) + #median upper and lower male
  as.layer(pNVh) + as.layer(pNVhl) + as.layer(pNVhu) #median upper and lower herm
dev.off()

#To plot curves on separate plots in one image
#opar <- par(mfrow = c(2,1), mar = c(3,3,1,1), oma = c(2,0,0,0))
#plot(t, Lt, t="l", xlab = "Relative age (years)", ylab = "Total length (cm)", mtext("LFD", side = 3, line =1))
#plot(t, L2, t="l", xlab = "Relative age (years)", ylab = "Total length (cm)", mtext("T&R", side = 3, line=2))
#par(opar)


###
###
# 3. PLOT ELEFAN_GA RESULTS ---------------------------------
###
###
path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

res <- read.table("res_MA3_Linfmin50_Sys.TimeOK.csv", header = T) 

#View(res)
names(res)

############### 
# univariate_density_w_output ------------------
##  Same  function as "univariate_density()", but also provides and output table

CI = 95
use_hist = FALSE
nbreaks = 10
mar = c(1.5,2,2,0)
oma = c(1.5,0,0,0.5)
mgp = c(2,0.5,0)
tcl = -0.25
cex = 1

tiff("LFA_linfmin50_ma3_univdens.tiff", height=1800, width=3200, res = 300) 

op <- par(no.readonly = TRUE)
par(
    mfcol = c(1, ncol(res)),
    mar = mar, oma = oma,
    mgp = mgp, tcl = tcl, cex = cex
  )
  
VARS <- list(
    Linf = expression(italic(L)[infinity]),
    K = expression(italic(K)),
    t_anchor = expression(italic(t)[anchor]),
    C = expression(italic(C)),
    ts = expression(italic(t)[s]),
    phiL = expression(paste(phi,"'"))
  )
  
Output_table <- as.data.frame(   matrix(NA, nrow = 4, ncol =ncol(res))  ) 
names(Output_table) <-  (names(res)) 
rownames(Output_table) <- c("CI_lower","CI_upper", "max_density", "median")
  
# univariate plots
for(i in seq(ncol(res))){
    x <- ks::kde(res[,i])
    
    h = hist(res[,i], plot=FALSE, breaks = nbreaks)
    
    xlim <- c(0, max(x$estimate))
    if(use_hist){
      xlim <- range(c(xlim, max(h$density)))
    }
    xlim <- xlim * c(0,1.1)
    
    plot(x$estimate, x$eval.points, t="n",
         xlim = xlim,
         xaxs = "i",
         ylab="", xlab="", col=1, lty=1
    )
    usr <- par()$usr
    
    CItxt <- paste0(round(100-CI), "%")
    inCI <- rle( x$estimate > x$cont[CItxt] )
    start.idx <- c(1, cumsum(inCI$lengths[-length(inCI$lengths)])+1)
    end.idx <- cumsum(inCI$lengths)
    limCI <- range(x$eval.points[start.idx[min(which(inCI$values))]:end.idx[max(which(inCI$values))]])
    
    in1 <- which(x$estimate > x$cont["99%"])
    mean1 <- mean(x$eval.points[in1])
    
    if(use_hist){
      rect(
        xleft = 0, ybottom = h$breaks[-length(h$breaks)],
        xright = h$density, ytop = h$breaks[-1],
        col = "grey90", border = "grey50"
      )
    }else{
      for(j in seq(inCI$lengths)){
        if(inCI$values[j]){
          polygon(
            y = c(x$eval.points[start.idx[j]:end.idx[j]], rev(x$eval.points[start.idx[j]:end.idx[j]])),
            x = c(x$estimate[start.idx[j]:end.idx[j]], rep(0, length(x$estimate[start.idx[j]:end.idx[j]]))),
            col = "grey90", #col = rgb(0,1,0,0.2),
            border = NA, lty = 3, lwd = 1
          )
        }
      }
    }
    
    # abline(v = x$cont[CItxt], lty=2, col="grey50")
    lines(x$estimate, x$eval.points, lwd = 1, col = "grey50")
    
    # rug
    segments(x0 = 0, x1 = par()$cxy[1]*0.3, y0 = x$x, y1 = x$x, col=rgb(0,0,0,0.5), lwd=0.3)
    
    # range of CI
    abline(h = limCI, lty = 1, lwd=1, col = 1)
    text(y =c(limCI), x = mean(usr[1:2]),
         labels = paste(sprintf("%.2f", round(c(limCI),2))),
         pos = c(1,3), offset = 0.25, col = 1
    )
    abline(h = mean1, lty = 1, lwd=1, col = 1)
    text(y =  mean1, x = mean(usr[1:2]),
         labels = sprintf("%.2f", round(mean1,2)),
         pos = 3,
         offset = 0.25, col = 1
    )
    
    varlab <- VARS[[match(names(res)[i], names(VARS))]]
    mtext(varlab, line=0.25, side=3)
    
    box()
    
    
    Output_table[1,i] <- limCI[1]
    Output_table[2,i] <- limCI[2] 
    Output_table[3,i] <- mean1
    Output_table[4,i] <- median(res[,i])
    
    
  }
mtext("Density", side = 1, line = 0, outer = TRUE)
par(op)

dev.off()

  
Output_table
  


########################################
####
# Now calculate the 95 % confidence envelope AND plot 95% CI from matrix
#####  
########################################

#library(scales)

#t.seq.p = seq(0,50,0.1)

# Make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.

L.mat <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(res$K)))

# fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.mat))) {
  
  Lt <- VBGF(list(Linf=res$Linf[w], 
                  K= res$K[w], 
                  t0= res$t_anchor[w], 
                  C = res$C[w],
                  ts = res$ts[w]),
  t=t.seq.p)
  
  
  L.mat[,w] <- round(Lt,2)
  
}

# Now for vector - LFA
lo.ci.vec <- as.numeric(apply(L.mat, 1, quantile, probs = 0.025))
up.ci.vec <- as.numeric(apply(L.mat, 1, quantile, probs = 0.975))
med.vec <-as.numeric(apply(L.mat, 1, quantile, probs = 0.5))

# Add together
list_df = list(t.seq.p, med.vec, lo.ci.vec, up.ci.vec)
areadf = list_df

view(areadf)
write.csv(areadf, "LFA_medUL_min30MA5.csv")

##NB Go to other file for plots and calcs

##########
#----------PLOT TAG-RECAPTURE CURVE - DENSITY PLOTS = DONE; NB: Linf first, then K
###########
path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

outp.mat.clean = read.csv("1000StoutputB1_9_from_Margit.csv", header =TRUE)

res <- outp.mat.clean

# Phi_prime = log10(K) + 2 log10(Linf)
res$Phi_vec <- log10(res$K) + 2 * log10(res$Linf)

#                K         Linf    Phi_vec
#CI_lower    -0.01158603  34.82568 2.404062
#CI_upper     0.41522135 191.21925 2.988057
#max_density  0.07118269  44.02530 2.723152
#median       0.18500000  49.70000 2.716971

#Get table of data from Tag-recapture data
# univariate plots
#THIS ONLY NEED TO RUN ONCE
CI = 95
use_hist = FALSE
nbreaks = 10
mar = c(1.5,2,2,0)
oma = c(1.5,0,0,0.5)
mgp = c(2,0.5,0)
tcl = -0.25
cex = 1

tiff(filename = "univdens_TR.tiff",  height=1800, width=3200, res = 300)

    op <- par(no.readonly = TRUE)
    par(
      # mfcol = c(floor(sqrt(ncol(res))), ceiling(sqrt(ncol(res)))),
      mfcol = c(1, ncol(res)),
      mar = mar, oma = oma,
      mgp = mgp, tcl = tcl, cex = cex
    )
    VARS <- list(
      Linf = expression(italic(L)[infinity]),
      K = expression(italic(K)),
      Phi_vec = expression(paste(phi,"'")))

    Output_table <- as.data.frame(matrix(NA, nrow = 4, ncol =ncol(res))  ) 
    names(Output_table) <-  (names(res)) 
    rownames(Output_table) <- c("CI_lower","CI_upper", "max_density", "median")
    
    # univariate plots
    for(i in seq(ncol(res))){
      x <- ks::kde(res[,i])
      
      h = hist(res[,i], plot=FALSE, breaks = nbreaks)
      
      xlim <- c(0, max(x$estimate))
      if(use_hist){
        xlim <- range(c(xlim, max(h$density)))
      }
      xlim <- xlim * c(0,1.1)
      
      plot(x$estimate, x$eval.points, t="n",
           xlim = xlim,
           xaxs = "i",
           ylab="", xlab="", col=1, lty=1
      )
      usr <- par()$usr
      
      CItxt <- paste0(round(100-CI), "%")
      inCI <- rle( x$estimate > x$cont[CItxt] )
      start.idx <- c(1, cumsum(inCI$lengths[-length(inCI$lengths)])+1)
      end.idx <- cumsum(inCI$lengths)
      limCI <- range(x$eval.points[start.idx[min(which(inCI$values))]:end.idx[max(which(inCI$values))]])
      
      in1 <- which(x$estimate > x$cont["99%"])
      mean1 <- mean(x$eval.points[in1])
      
      if(use_hist){
        rect(
          xleft = 0, ybottom = h$breaks[-length(h$breaks)],
          xright = h$density, ytop = h$breaks[-1],
          col = "grey90", border = "grey50"
        )
      }else{
        for(j in seq(inCI$lengths)){
          if(inCI$values[j]){
            polygon(
              y = c(x$eval.points[start.idx[j]:end.idx[j]], rev(x$eval.points[start.idx[j]:end.idx[j]])),
              x = c(x$estimate[start.idx[j]:end.idx[j]], rep(0, length(x$estimate[start.idx[j]:end.idx[j]]))),
              col = "grey90", #col = rgb(0,1,0,0.2),
              border = NA, lty = 3, lwd = 1
            )
          }
        }
      }
      
      # abline(v = x$cont[CItxt], lty=2, col="grey50")
      lines(x$estimate, x$eval.points, lwd = 1, col = "grey50")
      
      # rug
      segments(x0 = 0, x1 = par()$cxy[1]*0.3, y0 = x$x, y1 = x$x, col=rgb(0,0,0,0.5), lwd=0.3)
      
      # range of CI
      abline(h = limCI, lty = 1, lwd=1, col = 1)
      text(y =c(limCI), x = mean(usr[1:2]),
           labels = paste(sprintf("%.2f", round(c(limCI),2))),
           pos = c(1,3), offset = 0.25, col = 1
      )
      abline(h = mean1, lty = 1, lwd=1, col = 1)
      text(y =  mean1, x = mean(usr[1:2]),
           labels = sprintf("%.2f", round(mean1,2)),
           pos = 3,
           offset = 0.25, col = 1
      )
      
      varlab <- VARS[[match(names(res)[i], names(VARS))]]
      mtext(varlab, line=0.25, side=3)
      
      box()
      
      
      Output_table[1,i] <- limCI[1]
      Output_table[2,i] <- limCI[2] 
      Output_table[3,i] <- mean1
      Output_table[4,i] <- median(res[,i])
      
    }
    mtext("Density", side = 1, line = 0, outer = TRUE)
    par(op)
    
    dev.off()
    
    Output_list <- list (Output_table) 
    
  Output_table  
  
#---------------------------------
####
#Now calculate the 95 % confid. envelope AND plot 95% CI from matrix
#####
t.seq.p = seq(0,50,0.1)
  
library(scales)
  
# make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.

L.matTR <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(res$K)))
  
# fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.matTR))) {
    
    Lt <- VBGF(list(Linf=res$Linf[w], 
                    K= res$K[w], 
                    t0= 0),
               t=t.seq.p)
    
    
    L.matTR[,w] <- round(Lt,2)
    
}
  
#now for vector - T&R
lo.ci.vecTR <- as.numeric(apply(L.matTR, 1, quantile, probs = 0.025))
up.ci.vecTR <- as.numeric(apply(L.matTR, 1, quantile, probs = 0.975))
med.vecTR <-as.numeric(apply(L.matTR, 1, quantile, probs = 0.5))

#Add together
list_df = list(t.seq.p, med.vecTR, lo.ci.vecTR, up.ci.vecTR)
areadf = list_df

view(areadf)
write.csv(areadf, "T&R_med_UL.csv")
