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
library("dlnm")
library("splines")

varper <- c(2,10,25,50,75,90,98)
lag <- 29
lagnk <- 4
vb <- quantile(df3[91:120], varper/100,na.rm=T)
argvar0 <- list(fun="bs",degree=2,knots=vb)
arglag0 <- list(knots=logknots(lag,lagnk))

cb2 <- crossbasis(df3[91:120],lag=lag,argvar=argvar0, arglag=arglag0)
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

if(!exists("cblogit_preds")){
  source("~/Documents/dlnmRev/fit_indep_days_since.R")
}

varperf <- c(10,20,30,50,70,80,90)
vbf <- quantile(df3[91:120], varper/100,na.rm=T)

source("~/Documents/dlnmRev/custom_crossbasis_scatter.R")
environment(custom_crossbasis_scatter) <- asNamespace('dlnm')
totals = custom_crossbasis_scatter(final_df$curTemp,lag=lag,varbreaks=vbf,
                                   group=final_df$ID, pslags=(final_df$curFupDay+1),densReturn=TRUE)

xminmax0 = c(min(final_df$curTemp),max(final_df$curTemp),0.1)
#logknots(lag,lagnk)
arglag2 <- list(fun="bs",degree=2,knots=c( 0,1,2,5,12), intercept=FALSE)
source("~/Documents/dlnmRev/custom_crossbasis_gen.R")
environment(custom_crossbasis_gen) <- asNamespace('dlnm')
imprints = custom_crossbasis_gen(xminmax0,lag=lag, varbreaks = vbf, 
                                 argvar=argvar0, arglag=arglag2)
df_recomb_list <- data.frame(c(totals))
lagdiml <- length(arglag2$knots)+2
for (var in seq(dim(imprints)[1])){
  v_here = (var-1) %/% (lagdiml*lagdiml)
  l_here =((var-1) %% (lagdiml*lagdiml)) %/% lagdiml
  d_here = (var-1) %% (lagdiml)
  colname = paste0("v",v_here,".l",l_here,".d",d_here)
  df_recomb_list[colname] <- c(imprints[var,,,])
}

print("Crossbasis formed")
fitlm <- lm(c.totals. ~ ., data = df_recomb_list)
print("Interpolation 1 fitted")
#print(sum(is.na(fitlm$coefficients)))
coefs0 <- fitlm$coefficients[2:length(fitlm$coefficients)]
coefs0[is.na(coefs0)] <- 0

fittedarr = array(rep(fitlm$coefficients[1],length(totals)), dim = dim(totals))
for (var in seq(dim(imprints)[1])){
  fittedarr = fittedarr + imprints[var,,,]*coefs0[var]
}
dev.new()
library('plot.matrix')
dev.new()
plot(fittedarr[2,,],breaks=seq(90,120,2))
dev.new()
plot(totals[2,,],breaks=seq(90,120,2))


#events_df = final_df[final_df$Composite_Readmit_Mort == 1, ]

totals_counts = custom_crossbasis_scatter(final_df$curTemp,lag=lag,varbreaks=vbf,
                                   group=final_df$ID, pslags=(final_df$curFupDay+1))

totals_events_c = custom_crossbasis_scatter(final_df$curTemp,lag=lag,varbreaks=vbf,
                  group=final_df$ID, pslags=(final_df$curFupDay+1), dayweights = final_df$event)

totals_events = custom_crossbasis_scatter(final_df$curTemp,lag=lag,varbreaks=vbf,
      group=final_df$ID, pslags=(final_df$curFupDay+1), dayweights = final_df$event,densReturn=TRUE)
df_recomb_list$c.totals. = c(totals_events)
fitlm2 <- lm(c.totals. ~ ., data = df_recomb_list)
print("Interpolation 2 fitted")
coefsE <- fitlm2$coefficients[2:length(fitlm2$coefficients)]
coefsE[is.na(coefsE)] <- 0
fittedE = array(rep(fitlm2$coefficients[1],length(totals_events)), dim = dim(totals_events))
for (var in seq(dim(imprints)[1])){
  fittedE = fittedE + imprints[var,,,]*coefsE[var]
}
dev.new()
plot(fittedE[2,,],breaks=seq(0,8,1))
dev.new()
plot(totals_events[2,,],breaks=seq(0,8,1))

library("boot")
expect_pts0 = inv.logit(cblogit_preds$fit)
r30ratio = sum(expect_pts0,na.rm=TRUE)/sum(final_df$event)
expect_pts_high = inv.logit(cblogit_preds$fit+1.28*cblogit_preds$se.fit)/r30ratio
#using Kelly skewness to estimate fit distribution
expect_events_high =  custom_crossbasis_scatter(final_df$curTemp,lag=lag,varbreaks=vbf,
                              group=final_df$ID, pslags=(final_df$curFupDay+1), dayweights = expect_pts_high)
expect_pts_low = inv.logit(cblogit_preds$fit-1.28*cblogit_preds$se.fit)/r30ratio
expect_events_low =  custom_crossbasis_scatter(final_df$curTemp,lag=lag,varbreaks=vbf,
                              group=final_df$ID, pslags=(final_df$curFupDay+1),dayweights = expect_pts_low)
expect_pts = inv.logit(cblogit_preds$fit)/r30ratio
expect_events =  custom_crossbasis_scatter(final_df$curTemp,lag=lag,varbreaks=vbf,
                                               group=final_df$ID, pslags=(final_df$curFupDay+1),dayweights = expect_pts)

source("~/Documents/dlnmRev/prob_utils.R")
temp_shifts_cent = array(rep(0,length(totals)), dim = dim(totals_counts))
temp_shifts_high = array(rep(0,length(totals)), dim = dim(totals_counts))
temp_shifts_low = array(rep(0,length(totals)), dim = dim(totals_counts))
for(i in 1:dim(totals_counts)[1]) {
  for(j in 1:dim(totals_counts)[2]) {
    for(k in 1:dim(totals_counts)[3]) {
      if(totals_counts[i,j,k]>0){
      percdataHere = c(expect_events[i,j,k]/totals_counts[i,j,k], expect_events_low[i,j,k]/totals_counts[i,j,k], 
                       expect_events_high[i,j,k]/totals_counts[i,j,k])
      temp_shifts_cent[i,j,k] <- try(CmeasB_times_LNskewfit(percdataHere,1+totals_events_c[i,j,k],1+totals_counts[i,j,k]-totals_events_c[i,j,k]) )
      temp_shifts_high[i,j,k] <- try(CmeasB_times_LNskewfit(percdataHere,1+totals_events_c[i,j,k],1+totals_counts[i,j,k]-totals_events_c[i,j,k], P=0.975))
      temp_shifts_low[i,j,k] <- try(CmeasB_times_LNskewfit(percdataHere,1+totals_events_c[i,j,k],1+totals_counts[i,j,k]-totals_events_c[i,j,k], P=0.025))
    }}}}

for(j in 1:dim(totals_counts)[2]) {
  temp_shifts_cent[,j,] = temp_shifts_cent[,j,] / j #one patient will cover j number of lag days, so rate is averaged over all experienced
  temp_shifts_high[,j,] = temp_shifts_high[,j,] / j
  temp_shifts_low[,j,] = temp_shifts_low[,j,] / j
}

#dev.new()
#plot(temp_shifts_cent[2,,],breaks=seq(-0.2,0.2,0.1))
df_recomb_list$c.totals. = c(temp_shifts_cent)
shiftslm <- lm(c.totals. ~ ., data = df_recomb_list)
print("Interpolation to shifts fitted")
coefsS <- shiftslm$coefficients[2:length(shiftslm$coefficients)]
coefsS[is.na(coefsS)] <- 0
fittedS = array(rep(shiftslm$coefficients[1],length(totals_counts)), dim = dim(totals_counts))
for (var in seq(dim(imprints)[1])){
  fittedS = fittedS + imprints[var,,,]*coefsS[var]}

df_recomb_list$c.totals. = c(temp_shifts_high)
shiftslmhi <- lm(c.totals. ~ ., data = df_recomb_list)
print("Interpolation to shifts high fitted")
coefsShi <- shiftslmhi$coefficients[2:length(shiftslmhi$coefficients)]
coefsShi[is.na(coefsShi)] <- 0
fittedShi = array(rep(shiftslmhi$coefficients[1],length(totals_counts)), dim = dim(totals_counts))
for (var in seq(dim(imprints)[1])){
  fittedShi = fittedShi + imprints[var,,,]*coefsShi[var]}


df_recomb_list$c.totals. = c(temp_shifts_low)
shiftslmlo <- lm(c.totals. ~ ., data = df_recomb_list)
print("Interpolation to shifts low fitted")
coefsSlo <- shiftslmlo$coefficients[2:length(shiftslmlo$coefficients)]
coefsSlo[is.na(coefsSlo)] <- 0
fittedSlo = array(rep(shiftslmlo$coefficients[1],length(totals_counts)), dim = dim(totals_counts))
for (var in seq(dim(imprints)[1])){
  fittedSlo = fittedSlo + imprints[var,,,]*coefsSlo[var]}

#dev.new()
#plot(fittedS[2,,],breaks=seq(-0.2,0.2,0.1))

temp_day_cent = array(rep(0,dim(totals_counts)[1]*dim(totals_counts)[2]), dim =dim(totals_counts)[1:2])
temp_day_high = array(rep(0,dim(totals_counts)[1]*dim(totals_counts)[2]), dim =dim(totals_counts)[1:2])
temp_day_low = array(rep(0,dim(totals_counts)[1]*dim(totals_counts)[2]), dim =dim(totals_counts)[1:2])

for(j in 1:dim(totals_counts)[2]) {
  temp_shifts_cent[,j,] = temp_shifts_cent[,j,] / j #one patient will cover j number of lag days, so rate is averaged over all experienced
  temp_shifts_high[,j,] = temp_shifts_high[,j,] / j
  temp_shifts_low[,j,] = temp_shifts_low[,j,] / j
  for (i in 1:dim(totals_counts)[1]){
  temp_day_cent[i,j] = sum(temp_shifts_cent[i,j,])
  temp_day_high[i,j] = sum(temp_shifts_high[i,j,])
  temp_day_low[i,j] = sum(temp_shifts_low[i,j,])
  }
}

dev.new()
plot(temp_day_cent,breaks=seq(-0.6,0.6,0.1))

levels <- pretty(temp_day_cent, 20)
col1 <- colorRampPalette(c("blue", "white"))
col2 <- colorRampPalette(c("white", "red"))
col <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
varper2 <-c(5,(10+25)/2,(25+50)/2,(50+75)/2,(75+90)/2,95) #c(10,25,50,75,90)
vbf2= quantile(df3[91:120], varper2/100,na.rm=T)
(vbf-273.15) * 9/5 + 32
#filled.contour(x = vbf2, y =seq(1,30, 1),   z = temp_day_cent, col = col, levels = levels)

temp2_shifts_cent = array(rep(0,length(totals)/4), dim = c(dim(totals_counts)[1],dim(totals_counts)[2]/2,dim(totals_counts)[3]/2))
temp2_shifts_high = array(rep(0,length(totals)),dim = c(dim(totals_counts)[1],dim(totals_counts)[2]/2,dim(totals_counts)[3]/2))
temp2_shifts_low = array(rep(0,length(totals)), dim = c(dim(totals_counts)[1],dim(totals_counts)[2]/2,dim(totals_counts)[3]/2))
temp2_shifts_ps = array(rep(0,length(totals)), dim = c(dim(totals_counts)[1],dim(totals_counts)[2]/2,dim(totals_counts)[3]/2))
for(i in 1:dim(totals_counts)[1]) {
  for(j in 1:(dim(totals_counts)[2]/2)) {
    for(k in 1:(dim(totals_counts)[3]/2)) {
      if(totals_counts[i,j*2-1,k*2-1]>0){
      percdataHere = c(expect_events[i,j*2,k*2]+expect_events[i,j*2-1,k*2]+expect_events[i,j*2,k*2-1]+expect_events[i,j*2-1,k*2-1], 
                       expect_events_low[i,j*2,k*2]+expect_events_low[i,j*2-1,k*2]+expect_events_low[i,j*2,k*2-1]+expect_events_low[i,j*2-1,k*2-1], 
                       expect_events_high[i,j*2,k*2]+expect_events_high[i,j*2-1,k*2]+expect_events_high[i,j*2,k*2-1]+expect_events_high[i,j*2-1,k*2-1])
      tcount = (totals_counts[i,j*2,k*2]+totals_counts[i,j*2-1,k*2]+totals_counts[i,j*2,k*2-1]+totals_counts[i,j*2-1,k*2-1])
      percdataHere = percdataHere / tcount
      alp=1+totals_events_c[i,j*2,k*2]+totals_events_c[i,j*2-1,k*2]+totals_events_c[i,j*2,k*2-1]+totals_events_c[i,j*2-1,k*2-1]
      bet = 2 - alp + tcount
      temp2_shifts_cent[i,j,k] <- try(CmeasB_times_LNskewfit(percdataHere,alp,bet) )
      temp2_shifts_high[i,j,k] <- try(CmeasB_times_LNskewfit(percdataHere,alp,bet, P=0.975))
      temp2_shifts_low[i,j,k] <- try(CmeasB_times_LNskewfit(percdataHere,alp,bet, P=0.025))
      parms <- lognorm_params_kelly(percdataHere[1],percdataHere[2],percdataHere[3])
      temp2_shifts_ps[i,j,k]<- try(PmeasB_times_LNskew(1,parms[1],parms[2],parms[3],alp,bet,P=0))
    }}}}
temp2_day_cent = array(rep(0,dim(totals_counts)[1]*dim(totals_counts)[2]/2), dim =c(dim(totals_counts)[1],dim(totals_counts)[2]))
temp2_day_high = array(rep(0,dim(totals_counts)[1]*dim(totals_counts)[2]/2), dim =c(dim(totals_counts)[1],dim(totals_counts)[2]))
temp2_day_low = array(rep(0,dim(totals_counts)[1]*dim(totals_counts)[2]/2), dim =c(dim(totals_counts)[1],dim(totals_counts)[2]))

for(j in 1:dim(totals_counts)[2]/2) {
  temp2_shifts_cent[,j,] = temp2_shifts_cent[,j,] / (2*j-.5) #one patient will cover j number of lag days, so rate is averaged over all experienced
  temp2_shifts_high[,j,] = temp2_shifts_high[,j,] / (2*j-.5)
  temp2_shifts_low[,j,] = temp2_shifts_low[,j,] / (2*j-.5)
  for (i in 1:dim(totals_counts)[1]){
    temp2_day_cent[i,j] = sum(temp2_shifts_cent[i,j,]) *2
    temp2_day_high[i,j] = sum(temp2_shifts_high[i,j,]) *2
    temp2_day_low[i,j] = sum(temp2_shifts_low[i,j,]) *2
  }
}
temp2_shifts_ps[8,1:7,1:7]
temp2_shifts_cent[8,1:7,1:7]
mmx=c(min(df3[91:120],na.rm=TRUE),max(df3[91:120],na.rm=TRUE))
(mmx-273.15) * 9/5 + 32
sum(rowSums((df3[91:120]>mmx[1]) * (df3[91:120]<vbf[1]),  na.rm=TRUE) >0 , na.rm=TRUE)
## Re-run the glm(start = ...) one each post-op day and print out each iteration to see where changes / improvements are made