#
### Analysis of west coast steenbras (Lithugnathus aureti) body growth
### using bootstrapped length-at-age (otolith) analyses

# Margit Wilhelm and Ralf Schwamborn 
# August 2024
#
###
#
# This script contains analyses conducted for the manuscript:
#
# "Extremely slow somatic growth and intermittent recruitment of west coast 
#steenbras, an over-exploited, longevous Sparid, analysed with novel bootstrapped 
#methods"
#
# Authors of the manuscript:
#
#Margit Wilhelm 
#Arariky Shikongo 
#Angelika Veii
#Ralf Schwamborn
#Corresponding author: mwilhelm@unam.na

## 0. Clean memory ------------------------------------------------------------
gc()  
gc(reset=T)
rm(list = ls()) 
#Ctrl+Shift+F10 clean memory

opar <- par() # save plot parameters

## 0. Load packages -----------------------------------------------------------

 library(TropFishR)
 library(ks)
 library(rfishbase)
 
 library(devtools)

 #install_github("rschwamborn/fishboot")
 library(fishboot)

# 0. Function definition ---------------------------------------------------------

#### Define the function "lenage_boot" ----------------------------------
# lenage_boot
# Bootstrapped length-at-age analysis 
# input: lengths and ages (e.g., from otolith readings) 
# output: bootstrap posteriors of VGBF parameters
# part of the "fishboot" R package 

lenage_boot <- function(input.data, nboot = 1000) {
  
  data.len.age <- input.data
  
  B = nboot ## number of bootstraps
  res = data.frame(Linf = numeric(B) , K = numeric(B), t0 = numeric(B)) ## vector to hold results
  n = length(data.len.age$length) # number of data pairs
  
  for(b in 1:B){
  
    tryCatch({   
      
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
      
    }, error=function(e){})
    
  }
  
  ret <- list()
  
  ret$bootRaw <- res
  
  ret$seed <- seed
  
  class(ret) <- "lfqBoot"
 
  return(ret)
}

# 5. Bootstrapped length_at_age (from otoliths)   ----------------------------
 
### Function Len.age_boot  V.06
## otolith length_at_age bootstrap V.05
# Copyright: Ralf Schwamborn, 2020
# includes Len.age_boot function

# Analysis of West coast steenbras, L. aureti, otolith data

#####  Bootstrap the VBGF fit to length-at-age data  #####

#### load  data ####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

Otol_len_age_VN <-read.csv("age-length-rawVN.csv", header = T)
Otol_len_age_VS <-read.csv("age-length-rawVS.csv", header = T)

#-------------------------------------------
##Subset data by males, females, etc.
Otol_len_age_VN_f = subset(Otol_len_age_VN, Sex=="F")
Otol_len_age_VN_m = subset(Otol_len_age_VN, Sex=="M")
Otol_len_age_VN_h = subset(Otol_len_age_VN, Sex=="H")

library(tidyverse)

Otol_len_age_VN<-Otol_len_age_VN %>%
  mutate(Sex2 = ifelse(Sex %in% c("F","H"), "FH","MJ")) #create combined sex variable

Otol_len_age_VS<-Otol_len_age_VS %>%
  mutate(Sex2 = ifelse(Sex %in% c("F","H"), "FH","MJ")) #create combined sex variable

Otol_len_age_VN_fh = subset(Otol_len_age_VN, Sex2=="FH")
Otol_len_age_VN_mj = subset(Otol_len_age_VN, Sex2=="MJ")
Otol_len_age_VS_fh = subset(Otol_len_age_VS, Sex2=="FH")
Otol_len_age_VS_mj = subset(Otol_len_age_VS, Sex2=="MJ")

Otol_len_age_VS_f = subset(Otol_len_age_VS, Sex=="F")
Otol_len_age_VS_m = subset(Otol_len_age_VS, Sex=="M")
Otol_len_age_VS_h = subset(Otol_len_age_VS, Sex=="H")

#View(Otol_len_age_VN_h)
length(Otol_len_age_VN_h$length_cm) #North: 30 females, 33 males, 13 H
length(Otol_len_age_VS_f$length_cm) #South: 48 females, 38 males, 67 H

attach(Otol_len_age_VN)
plot(length_cm ~ age_yr)

## adjust a VBGF function

# non linear least squares method , TropFishR::growth_length_age, method = "LSM"
# for this you need library(TropFishR)

#View(Otol_len_age_VN)
output1 <- growth_length_age(param = Otol_len_age_VS, method = "LSM", Linf_init = 50, CI = TRUE, age_plot=NULL) #Some of them need more iterations. Only get those that work when bootstrapping
output2 <- growth_length_age(param = Otol_len_age_VN, method = "LSM", Linf_init = 50, CI = TRUE, age_plot=NULL)
output3 <- growth_length_age(param = Otol_len_age_VN_f, method = "LSM", Linf_init = 50, CI = TRUE, age_plot=NULL)
output4 <- growth_length_age(param = Otol_len_age_VN_m, method = "LSM", Linf_init = 50, CI = TRUE, age_plot=NULL)
output5 <- growth_length_age(param = Otol_len_age_VN_h, method = "LSM", Linf_init = max(length_cm), CI = TRUE, age_plot=NULL)
output6 <- growth_length_age(param = Otol_len_age_VS_f, method = "LSM", Linf_init = max(length_cm), CI = TRUE, age_plot=NULL)
output7 <- growth_length_age(param = Otol_len_age_VS_m, method = "LSM", Linf_init = max(length_cm), CI = TRUE)
output8 <- growth_length_age(param = Otol_len_age_VS_h, method = "LSM", Linf_init = max(length_cm), CI = TRUE)

# 
summary(output2$mod)

?growth_length_age
# growth_length_age(steenbras_lenage, method = "GullandHolt")

# Bertalaffy plot
#growth_length_age(Lsynagris_lenage, method = "BertalanffyPlot", Linf_est = (max(length)))

# non linear least squares method
output <- growth_length_age(Otol_len_age_VN, method = "LSM",
                              CI = TRUE, age_plot=NULL,
                            Linf_init = (max(length_cm)))

summary(output$mod) 


## now bootstrap the VBGF with length-at-age  ###
#### DEFINE Function lenage_boot ####

lenage_boot <- function(input.data, nboot = 1000) {

  data.len.age <- input.data

  B = nboot ## number of bootstraps
  res = data.frame(Linf = numeric(B) , K = numeric(B), t0 = numeric(B)) ## vector to hold results
  n = length(data.len.age$length) # number of data pairs


  for(b in 1:B){

    tryCatch({   
  
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
  
  }, error=function(e){})
  
}

ret <- list()

ret$bootRaw <- res

ret$seed <- seed

class(ret) <- "lfqBoot"


return(ret)

}

####################
####################
###Start here every time and change data
####################
####################
resok <- lenage_boot(Otol_len_age_VN, nboot = 1000) #VN works best

res <- resok$bootRaw

hist(res$Linf)    #Check 0 values = undefined. Clean later (see below)
summary(res$Linf)
quantile(res$Linf, probs = c(0.025,0.5,0.975)) # 95% CIs

hist(res$K)
summary(res$K)
quantile(res$K, probs = c(0.025,0.5,0.975)) # 95% CIs

hist(res$t0)
summary(res$t0)
quantile(res$t0, probs = c(0.025,0.5,0.975)) # 95% CIs

output <- growth_length_age(param = Otol_len_age_VN, method = "LSM",Linf_init = 70, CI = TRUE, age_plot=NULL)

summary(output$mod)

# Works OK!
# Conclusion: Bootstrap 95% CI are much larger (CI width more than double) than nls (TropFishR) estimates

######
# plot otolith results, nboot = 1000 ###########

data.len.age <- Otol_len_age_VN

attach(Otol_len_age_VN)

###
# plot with grey curve swarms
###
# library(TropFishR)
# calculate Lt from t, with t0
EXP_Otol_len_age_400 <-resok$bootRaw #if not read in after already saved

# # clean (remove Linf = zero)
# 
EXP_Otol_len_age_400 <- EXP_Otol_len_age_400[!(EXP_Otol_len_age_400$Linf==0 ),]
# 
#View(EXP_Otol_len_age_400)
# 
length(EXP_Otol_len_age_400$Linf) #n=1000 if no non-fitted data (if not 1000, zero values removed. Need to bootstrap again and paste together)
res <- EXP_Otol_len_age_400

write.csv(res, "Veii_N_results1000.csv") #Save all results in new file. Merge files later if > 1 bootstrap

#View(res)

##### PLOT WITH GREAY CURVE SWARMS

t.seq2 <- seq(0, 30, by = 0.1)
Linf.2 <-  round((output2$Linf),2)
K.2 <- round((output2$K),3)
t0.2 <- output2$t0
L.seq.2 <- VBGF(t = t.seq2, list(Linf=Linf.2, K=K.2, C= 0, t0 = t0.2))  

plot(t.seq2, L.seq.2, t="l", col = "blue", xlim = c(0, 30), ylim = c(0, 120),
     xlab = "Age (y)", ylab = "Fork length (cm)")
points(data.len.age$length ~ data.len.age$age, col = "navyblue")


###
# plot with grey curve swarms
###

library(scales)

for (w in  1:length(res$K)) {

  t.seq.p = seq(0,max(data.len.age$age_yr),0.1)
  Lt.p <- VBGF(list(Linf=res$Linf[w],
                    K= res$K[w],
                    t0= res$t0[w]),
               t=t.seq.p)

  lines(t.seq.p, Lt.p, col = alpha("grey", 0.1), lwd = 0.3 ) #grey transparent lines

}

points(data.len.age$length_cm ~ data.len.age$age_yr, col = "navyblue")
lines(t.seq2, L.seq.2, t="l", col = "blue")

### nice plot
#t.seq.p

##### 
### calculate the 95 % confid. envelope 
#####

# make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.

L.mat <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(res$K)))

# fill the matrix with VBGF curves (length data for each age 't') from bootstrapping
for (w in 1:(ncol(L.mat))) {
  
  Lt <- VBGF(list(Linf=res$Linf[w], 
                  K= res$K[w], 
                  t0= res$t0[w]),
             t=t.seq.p)
  
  
  L.mat[,w] <- round(Lt,2)
  
}
View(L.mat)

# calculate 95% quantiles for each age 't'
#now for vector - all
lo.ci.vec <- as.numeric(apply(L.mat, 1, quantile, probs = 0.025))
up.ci.vec <- as.numeric(apply(L.mat, 1, quantile, probs = 0.975))
med.vec <-as.numeric(apply(L.mat, 1, quantile, probs = 0.5))

# insert CIs & median into the plot
lines(t.seq.p, med.vec, t="l", col = "red")
points(data.len.age$length ~ data.len.age$age, col = "navyblue")
lines(t.seq.p, lo.ci.vec, lty = 3, col = "navyblue", lwd = 1.5)
lines(t.seq.p, up.ci.vec, lty = 3, col = "navyblue", lwd = 1.5)

# Kimura plot ---------------------------

plot(res$K ~ res$Linf , ylim = c(0,0.15), xlim = c(70,100), col = alpha("blue", 0.2))

#fishboot style
#library(fishboot)

res2 <- list(bootRaw = res)
#res2$bootRaw$Linf

LinfK_scatterhist(res2)

univar_results <- univariate_density(res2)

univariate_density_w_output(res2) #Need to define below

# univariate_density_w_output ------------------
##  Same  function as "univariate_density()", but also provides and output table

univariate_density_w_output <- function(res, CI=95, use_hist = FALSE, nbreaks = 10,
                               mar = c(1.5,2,2,0), oma = c(1.5,0,0,0.5),
                               mgp = c(2,0.5,0), tcl = -0.25, cex = 1,
                               ...
){
  
  res <- res$bootRaw
 
  
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
  
 
 
   Output_list <- list ( Output_table) 
  
  
    return (Output_list)
  
  
}
