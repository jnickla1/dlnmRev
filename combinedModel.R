library("readxl")
library("dlnm")
library("splines")
library("glmtoolbox")

my_data <- read_excel("~/Documents/deidentified_Jul28.xlsx")
delcols <-c(30,33,35)
fcols <-c(3:4,6,8:28,37:45,50,52:55,57:75)
df2 <- my_data[ -c(1,3:7) ]
df2[fcols] <- lapply(df2[fcols], factor)
df2 <- df2[ -delcols ]
df2$dow <-as.factor(my_data$Disch_Day_of_Week)
df2$month <- my_data$Disch_Month
df2$ID <- seq.int(nrow(df2))-1
alllogit <- glm(Composite_Readmit_Mort ~ ., data = df2, family = "binomial",na.action=na.exclude)
fit2 <- stepCriterion(alllogit, criterion="bic")
model2 <- glm(paste("Composite_Readmit_Mort",fit2$final), data = df2, family = "binomial",na.action=na.exclude)

tempers <- read.csv("~/Documents/Mar8_envs_tem.csv",sep=',',header=TRUE)
dftempall <- merge(tempers,df2,by="ID")
dftempall$tmean <- rowMeans(dftempall[2:31],na.rm=TRUE)
for (x in 2:31) {
  dftempall[,x] <- ifelse(is.na(dftempall[,x]), dftempall$tmean, dftempall[,x])
}

varper <- c(10,25,50,75,90)
lag <- 29
lagnk <- 4
argvar <- list(fun="bs",degree=2,knots=quantile(dftempall[2:31], varper/100,na.rm=T))

cb2 <- crossbasis(dftempall[2:31],lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)))
ns2 <- onebasis(dftempall$month, fun='ns',df=8)
modelcomb <- glm(paste("Composite_Readmit_Mort",fit2$final,"+ dow + cb2 + ns2"), data = dftempall, family = "binomial",na.action=na.exclude)
cp2 <- crosspred(cb2, modelcomb, cen=mean(dftempall$tmean,na.rm=T),from=20, to=80, lag=c(0,30), by=1)
dev.new()
cp2crop <- crosspred(cb2, modelcomb, cen=mean(dftempall$tmean,na.rm=T),from=20, to=80, lag=c(0,20), by=1)
#plot(cp2crop, xlab="Temperature °F", zlab="RR", ylab="days post-discharge", theta=200, phi=40, lphi=30, main="3D graph of temperature effect")
plot(cp2crop, "contour", xlab="Temperature °F", key.title=title("RR"))
dev.new()
plot(cp2,"overall",ci.level=0.95,ylim=c(0.17,6),log="y",ylab="RR", xlab= "Temperature °F")
title("Overall with surgical predictors")
#dev.new()
#plot(cp2,"slices", lag=2,ci.level=0.95,ylim=c(0.33,3),log="y",ylab="RR", xlab= "Temperature °F")
#title("Lag=0 with surgical predictors")
dev.new()
plot(cp2,"slices", lag=1,ci.level=0.95,ylim=c(0.33,3),log="y",ylab="RR", xlab= "Temperature °F")
title("Lag=1 with surgical predictors")
dev.new()
plot(cp2,"slices", lag=2,ci.level=0.95,ylim=c(0.33,3),log="y",ylab="RR", xlab= "Temperature °F")
title("Lag=2 with surgical predictors")
dev.new()
np <- crosspred(ns2,modelcomb)
plot(np,ci.level=0.95,ylim=c(0.33,3),log="y",ylab="RR", xlab= "Month of Surgery")
title("Seasonal Effect")

library(pROC)
dev.new()
test_prob1 = predict(model2, type = "response")
test_roc1 = roc(df2$Composite_Readmit_Mort ~ test_prob1, plot = TRUE, print.auc = FALSE,col='blue')
test_prob0 = predict(modelcomb, type = "response")
test_roc0 = roc(dftempall$Composite_Readmit_Mort ~ test_prob0, plot = TRUE, print.auc = TRUE,col='black',add=TRUE)

dev.new()
arv<-confint(modelcomb)
arv2<-merge(arv,coef(summary(modelcomb)),by='row.names')
rownames(arv2)<- arv2$Row.names
arv2 <- arv2[,-1]
rowsel <- coef(summary(model2))
rowsel <- rowsel[ ,c(1,4) ]
colnames(rowsel) <- c("Est0","Pinit")
arv2<-merge(arv2,rowsel,by='row.names')
rownames(arv2)<- arv2$Row.names
arv2 <- arv2[-1,]
sdHcT<-sd(df2$Last_Hematocrit,na.rm=T)
arv2['Last_Hematocrit',1] <- paste("Last_Hematocrit +1sd (", sprintf("%0.2f", sdHcT) ,")")
arv2['Last_Hematocrit',2:4] <-arv2['Last_Hematocrit',2:4]* sdHcT
arv2['Last_Hematocrit',8] <-arv2['Last_Hematocrit',8]* sdHcT
colnames(arv2)=c("Row.names","lower","upper","v","stderr","zval","pval","est0","pinit")
arv2 <- arv2[order(arv2$pval,decreasing=TRUE),]
arv2$ID <- seq.int(nrow(arv2))
library(ggplot2)
sp <- ggplot(arv2, aes(x=v, y=factor(ID))) +geom_point(aes(x = v, y = factor(ID)),col='black') + geom_point(aes(x = est0, y = factor(ID)),col='blue')+geom_errorbar(aes(xmin = lower, xmax = upper))
print(sp + scale_y_discrete(breaks=arv2$ID, labels=arv2$Row.names)+xlab("Logit Score") + ylab("") +geom_vline(xintercept=0,linetype=3))
