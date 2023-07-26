library("readxl")

tempers <- read.csv("~/Documents/Mar8_envs_tem.csv",sep=',',header=TRUE)

my_data <- read_excel("~/Documents/deidentified_Mar8.xlsx")
df <- my_data[ c(1:2) ]

dftemp <- merge(tempers,df,by="ID")
dftemp$tmean <- rowMeans(dftemp[2:31],na.rm=TRUE)
for (x in 2:31) {
  dftemp[x,] <- ifelse(is.na(dftemp[x,]), dftemp$tmean, dftemp[x,])
}

library("dlnm")
varper <- c(10,25,50,60, 75,90)
lag <- 29
lagnk <- 7
argvar <- list(fun="bs",degree=2,knots=quantile(dftemp[2:31], varper/100,na.rm=T))



cb <- crossbasis(dftemp[2:31],lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)))
summary(cb)
  

cblogit <- glm(dftemp$Composite_Readmit_Mort ~ cb, family = "binomial")
cp <- crosspred(cb, cblogit, cen=mean(dftemp$tmean,na.rm=T),from=30, to=80, lag=c(0,10), by=1)
#plot(pred3.temp, xlab="Temperature", zlab="RR", theta=200, phi=40, lphi=30, main="3D graph of temperature effect")
#dev.new()
#fix(plot.crosspred); l155: levels <- exp(pretty(log(x$matfit), 10))
plot(cp, "contour", xlab="Temperature", key.title=title("RR"),
  plot.title=title("Contour plot",xlab="Temperature °F",ylab="Lag"))
dev.new()
plot(cp,"slices", lag=0,ci.level=0.95,ylim=c(0.2,6),log="y",ylab="RR", xlab= "Temperature °F")
title("Lag=0")
dev.new()
plot(cp,"slices", lag=1,ci.level=0.95,ylim=c(0.2,6),log="y",ylab="RR", xlab= "Temperature °F")
title("Lag=1")
dev.new()
plot(cp,"slices", var=77,ci.level=0.95,ylim=c(0.2,6),log="y",ylab="RR", xlab= "Lag")
title("Lagged Risk Profile at 77°F")
