library("readxl")
my_data <- read_excel("~/Documents/dlnmRevData/deidentified_Mar8.xlsx")
df <- my_data[ -c(1,3:5) ]
delcols <-c(30,33,35)
fcols <-c(3:4,6,8:28,37:45,50,52:55,57:83)
df[fcols] <- lapply(df[fcols], factor)
df <- df[ -delcols ]
#colnames(df)
mylogit <- glm(Composite_Readmit_Mort ~ ., data = df, family = "binomial",na.action=na.exclude)


library("glmtoolbox")
fit1 <- stepCriterion(mylogit,direction="forward", criterion="p-value")
#10 covariates
fit2 <- stepCriterion(mylogit, criterion="bic")
#5 covariates, disagree on Valve, only RF-Last Hematocrit is non-categorical
model1 <- glm(paste("Composite_Readmit_Mort",fit1$final), data = df, family = "binomial",na.action=na.exclude)
summary(model1)
dev.new()
plot(model1, which=1)

library(ggplot2)
dev.new()
model2 <- glm(paste("Composite_Readmit_Mort",fit2$final), data = df, family = "binomial",na.action=na.exclude)
summary(model2)
plot(model2, which=1)
dev.new()
plot(df$Last_Hematocrit, df$Composite_Readmit_Mort)
#popvp<-mean(my_data$Post_Op_Pulm_Vent_Prolonged,na.rm=T)
#possi<-mean(my_data$Post_Op_Surgical_Site_Infection,na.rm=T)
#clda<-mean(my_data$Chronic_Lung_Disease_Adj,na.rm=T)
#va<-mean(my_data$Valve,na.rm=T)
tinymodel=glm(Composite_Readmit_Mort ~ Last_Hematocrit, data = df, family = "binomial")
fhist = hist(df[df$Composite_Readmit_Mort == 0,]$Last_Hematocrit, plot=FALSE)
thist = hist(df[df$Composite_Readmit_Mort == 1,]$Last_Hematocrit, plot=FALSE)
rect(xleft=fhist$mids-2,ybottom=0,xright=fhist$mids+2,ytop=fhist$counts/nrow(df),col=gray(0.5))
rect(xleft=thist$mids-2,ybottom=1-thist$counts/nrow(df),xright=thist$mids+2,ytop=1,col=gray(0.5))
curve(predict(tinymodel, data.frame(Last_Hematocrit=x), type="response"), col="red", lwd=3, add=TRUE)

dev.new()
library(pROC)
test_prob2 = predict(model2, type = "response")
#cresults=df[!is.na(df$Last_Hematocrit) & !is.na(df$Post_Op_Pulm_Vent_Prolonged) & !is.na(df$Post_Op_Surgical_Site_Infection)
#     & !is.na(df$Chronic_Lung_Disease_Adj) & !is.na(df$Valve), ]$Composite_Readmit_Mort 
test_roc2 = roc(df$Composite_Readmit_Mort ~ test_prob2, plot = TRUE, print.auc = FALSE,col='red')
print(test_roc2)
test_prob1 = predict(model1, type = "response")
test_roc1 = roc(df$Composite_Readmit_Mort ~ test_prob1, plot = TRUE, print.auc = FALSE,col='blue',add=TRUE)
print(test_roc1)
test_prob0 = predict(mylogit, type = "response")
test_roc0 = roc(df$Composite_Readmit_Mort ~ test_prob0, plot = TRUE, print.auc = TRUE,col='black',add=TRUE)
print(test_roc0)

dev.new()
arv<-confint(model2)
arv2<-merge(arv,coef(summary(model2)),by='row.names')
rownames(arv2) <- arv2$Row.names
#arv2 <- arv2[-1,]
sdHcT<-sd(df$Last_Hematocrit,na.rm=T)
arv2['Last_Hematocrit',1] <- paste("Last_Hematocrit +1sd (", sprintf("%0.2f", sdHcT) ,")")
arv2['Last_Hematocrit',2:4] <-arv2['Last_Hematocrit',2:4]* sdHcT
colnames(arv2)=c("Row.names","lower","upper","v","stderr","zval","pval")
arv2 <- arv2[order(arv2$pval,decreasing=TRUE),]
arv2$ID <- seq.int(nrow(arv2))
sp <- ggplot(arv2, aes(x=v, y=factor(ID))) +geom_point()+geom_errorbar(aes(xmin = lower, xmax = upper))
print(sp + scale_y_discrete(breaks=arv2$ID, labels=arv2$Row.names)+xlab("Logit Score") + ylab("") +geom_vline(xintercept=0,linetype=3))


