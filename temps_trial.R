library("readxl")

tempers <- read.csv("~/Documents/dlnmRevData/Mar8_envs_tem.csv",sep=',',header=TRUE)

my_data <- read_excel("~/Documents/dlnmRevData/deidentified_Mar8.xlsx")
df <- my_data[ c(1:2) ]

dftemp <- merge(tempers,df,by="ID")
dftemp$tmean <- rowMeans(dftemp[2:31],na.rm=TRUE)
for (x in 2:31) {
  dftemp[,x] <- ifelse(is.na(dftemp[,x]), dftemp$tmean, dftemp[,x])
}

library("dlnm")
varper <- c(10,25,75,90)
lag <- 29
lagnk <- 4
argvar <- list(fun="bs",degree=2,knots=quantile(dftemp[2:31], varper/100,na.rm=T))
#argvar <- list(fun="bs",degree=2,knots=quantile(dftemp$tmean, varper/100,na.rm=T))


cb <- crossbasis(dftemp[2:31],lag=lag,argvar=argvar, arglag=list(knots=logknots(20,lagnk)))
summary(cb)
  
#plot(pred3.temp, xlab="Temperature", zlab="RR", theta=200, phi=40, lphi=30, main="3D graph of temperature effect")
#dev.new()
#fix(plot.crosspred); l155: levels <- exp(pretty(log(x$matfit), 10))

cblogit <- glm(dftemp$Composite_Readmit_Mort ~ cb, family = "binomial")
cp <- crosspred(cb, cblogit, cen=mean(dftemp$tmean,na.rm=T),from=20, to=80, lag=c(0,30), by=1)
cpcrop <- crosspred(cb, cblogit, cen=mean(dftemp$tmean,na.rm=T),from=20, to=80, lag=c(0,20), by=1)
dev.new()
plot(cpcrop, "contour", xlab="Temperature °F", key.title=title("RR"),
  plot.title=title("Contour plot",xlab="Temperature °F",ylab="Lag"))
dev.new()
plot(cp,"overall",ci.level=0.95,ylim=c(0.2,6),log="y",ylab="RR", xlab= "Temperature °F")
title("Overall")
dev.new()
plot(cp,"slices", lag=0,ci.level=0.95,ylim=c(0.2,6),log="y",ylab="RR", xlab= "Temperature °F")
title("Lag=0")
dev.new()
plot(cp,"slices", lag=1,ci.level=0.95,ylim=c(0.2,6),log="y",ylab="RR", xlab= "Temperature °F")
title("Lag=1")
dev.new()
plot(cp,"slices", var=74,ci.level=0.95,ylim=c(0.2,6),log="y",ylab="RR", xlab= "Lag")
title("Lagged Risk Profile at 74°F")
dev.new()
plot(cp,"slices", var=65,ci.level=0.95,ylim=c(0.2,6),log="y",ylab="RR", xlab= "Lag",add=TRUE, color="green")
title("Lagged Risk Profile at 65°F")
