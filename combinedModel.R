#must run glm_wo_temperature and temps_trial first
df2 <- my_data[ -c(1,3:5) ]
df2[fcols] <- lapply(df2[fcols], factor)
df2 <- df2[ -delcols ]
df2$ID <- seq.int(nrow(df2))
dftempall <- merge(tempers,df2,by="ID")
cb2 <- crossbasis(dftempall[2:31],lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)))
modelcomb <- glm(paste("Composite_Readmit_Mort",fit2$final,"+ cb2"), data = dftempall, family = "binomial",na.action=na.exclude)
cp2 <- crosspred(cb2, modelcomb, cen=mean(dftemp$tmean,na.rm=T),from=20, to=80, lag=c(0,30), by=1)
#cp2crop <- crosspred(cb2, modelcomb, cen=mean(dftemp$tmean,na.rm=T),from=20, to=80, lag=c(0,20), by=1)
#plot(cp2crop, "contour", xlab="Temperature °F", key.title=title("RR"))
plot(cp2,"overall",ci.level=0.95,ylim=c(0.2,6),log="y",ylab="RR", xlab= "Temperature °F")
