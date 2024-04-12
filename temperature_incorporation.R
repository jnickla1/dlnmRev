## Re-run computation with "fudged" inclusion of temperature copying the mean for days not observed
library("readxl")
my_data <- read.csv("~/Documents/dlnmRevData/PostCriteria_withTempPoll_DEIDENTIFIED3a.csv",sep=',',header=TRUE)
df <- my_data[ -c(1,144,147) ] #remove first ID col
#3-6 is blood prods
delcols <-c(3:6, 20, 60, 89, 138:139) #,18,40,42,47,89, 98:127,138:139)
fcols <-c(1, 9:12,14,17:20,22:24,28,30:31,33:69,71:73,78,80:81,83:86,89:92,132:143) #leave off DOW_month
df[fcols] <- lapply(df[fcols], factor)
df3 <- df[ -delcols ]


df3$tmean <- rowMeans(df3[91:120],na.rm=TRUE)
for (x in 91:120) {
  df3[,x] <- ifelse(is.na(df3[,x]), df3$tmean, df3[,x])
}

varper <- c(10,25,50,75,90)
lag <- 29
lagnk <- 4
argvar <- list(fun="bs",degree=2,knots=quantile(df3[91:120], varper/100,na.rm=T))
library("dlnm")
library("splines")
cb2 <- crossbasis(df3[91:120],lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)))
ns2 <- onebasis(df3$DOW_Month, fun='ns',df=6)

if(!exists("fit2")){
  source("~/Documents/dlnmRev/glm_wo_temperature_new.R")
}

model00 <- glm(paste("Composite_Readmit_Mort",fit2$final,"+ DOW_Discharge + ns2"), data = df3, family = "binomial",na.action=na.exclude)
summary(model00)
#season has no effect without considering temperature
modelcomb <- glm(paste("Composite_Readmit_Mort",fit2$final,"+ DOW_Discharge + cb2 + ns2"), data = df3, family = "binomial",na.action=na.exclude)
cp2 <- crosspred(cb2, modelcomb, cen=mean(df3$tmean,na.rm=T),from=270, to=300, lag=c(0,20), by=1)
dev.new()
cp2crop <- crosspred(cb2, modelcomb, cen=mean(df3$tmean,na.rm=T),from=270, to=300, lag=c(0,20), by=1)
plot(cp2crop, "contour", xlab="Temperature K", key.title=title("RR"))
#library('plot.matrix')
#plot(log(cp2crop$matRRfit))
dev.new()
plot(cp2,"overall",ci.level=0.95,log="y",ylab="RR", xlab= "Temperature K")
title("Overall with surgical predictors")
dev.new()
plot(cp2,"slices", lag=0,ci.level=0.95,ylim=c(0.33,3),log="y",ylab="RR", xlab= "Temperature K")
title("Lag=0 with surgical predictors")
dev.new()
plot(cp2,"slices", lag=1,ci.level=0.95,ylim=c(0.33,3),log="y",ylab="RR", xlab= "Temperature K")
title("Lag=1 with surgical predictors")
dev.new()
plot(cp2,"slices", lag=2,ci.level=0.95,ylim=c(0.33,3),log="y",ylab="RR", xlab= "Temperature K")
title("Lag=2 with surgical predictors")
dev.new()
np <- crosspred(ns2,modelcomb)
plot(np,ci.level=0.95,ylim=c(0.33,3),log="y",ylab="RR", xlab= "Month of Surgery")
title("Seasonal Effect")
#worse significance if we don't account for the seasonal effect

## Compute estimate of effect sizes with temperature on expanded model, dividing the logit values by day#+1

if(!exists("cblogit")){
  source("~/Documents/dlnmRev/fit_indep_days_since.R")
}


## Re-so the glm(start = ...) one each post-op day and print out each iteration to see where changes / improvements are made