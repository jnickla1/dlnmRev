#source("/Documents/dlnmRev/combinedModel.R")


#calculate area
area <- (sum(cp2$matRRhigh < 1) + sum(cp2$matRRlow > 1))/length(cp2$matRRhigh)
areas= data.frame( nsamps = numeric(),areafrac= numeric())
areas[1, ] <- c(nrow(dftempall),area)


#resample
iseq <- seq(7, 12, by=0.05)
r <- 2
for (i in iseq) {
nsamps <- floor(exp(i))
dfresamp <- dftempall[sample(nrow(dftempall), nsamps,replace=T), ]

cb2r <- crossbasis(dfresamp[2:31],lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)))
ns2r <- onebasis(dfresamp$month, fun='ns',df=4)
modelcombr <- glm(paste("Composite_Readmit_Mort",fit2$final,"+ dow + cb2r + ns2r"), data = dfresamp, family = "binomial",na.action=na.exclude)
cp2r <- crosspred(cb2r, modelcombr, cen=mean(dfresamp$tmean,na.rm=T),from=20, to=80, lag=c(0,30), by=1)
arear <- (sum(cp2r$matRRhigh < 1) + sum(cp2r$matRRlow > 1))/length(cp2r$matRRhigh)
areas[r, ] <- c(nsamps,arear)
print(i)
r <- r+1
}