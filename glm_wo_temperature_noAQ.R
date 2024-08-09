library("readxl")

my_data <- read.csv("~/Documents/dlnmRevData/PostCriteria_withTempPoll_DEIDENTIFIED3ab.csv",sep=',',header=TRUE)
df <- my_data[ -c(1,144,147) ] #remove first ID col
#3-6 is blood prods
delcols <-c(3:6, 8, 20, 60, 89,93:97,  98:127, 138:139) #,18,40,42,47,89, 98:127,138:139) #remove temperature initially
fcols <-c(1, 9:12,14,17:20,22:24,28,30:31,33:69,71:73,78,80:81,83:86,89:92,132:144)
df[fcols] <- lapply(df[fcols], factor)
df <- df[ -delcols ]
df0 = df
#colnames(df)
#library(fastglm)
mylogit <- glm( Composite_Readmit_Mort ~ . ,  data=df0, family = binomial(),na.action=na.exclude)


library("glmtoolbox")
fit2 <- stepCriterion(mylogit,criterion="p-value",levels=c(0.05,0.05),test="lr")
#9 covariates
fit1 <- stepCriterion(mylogit, criterion="bic")
#4 covariates
fit3 <- stepCriterion(mylogit, criterion="aic")
#4 covariates, disagree on Valve, only RF-Last Hematocrit is non-categorical
if (0 == 0){ #sys.nframe() == 0){

model1 <- glm(paste("Composite_Readmit_Mort",fit1$final), data = df0, family = "binomial",na.action=na.exclude)
summary(model1)
dev.new()
plot(model1, which=1)

library(ggplot2)
dev.new()
model2 <- glm(paste("Composite_Readmit_Mort",fit2$final), data = df0, family = "binomial",na.action=na.exclude)
summary(model2)
plot(model2, which=1)
dev.new()
#plot(df$Post.Op.Surgical.Site.Infection1, df$Composite_Readmit_Mort)

model3 <- glm(paste("Composite_Readmit_Mort",fit3$final), data = df0, family = "binomial",na.action=na.exclude)
summary(model3)
test_prob3 = predict(model3, type = "response")


#popvp<-mean(my_data$Post_Op_Pulm_Vent_Prolonged,na.rm=T)
#possi<-mean(my_data$Post_Op_Surgical_Site_Infection,na.rm=T)
#clda<-mean(my_data$Chronic_Lung_Disease_Adj,na.rm=T)
#va<-mean(my_data$Valve,na.rm=T)
#tinymodel=glm(Composite_Readmit_Mort ~ Post.Op.Surgical.Site.Infection1, data = df, family = "binomial")
#fhist = hist(df[df$Composite_Readmit_Mort == 0,]$Post.Op.Surgical.Site.Infection1, plot=FALSE)
#thist = hist(df[df$Composite_Readmit_Mort == 1,]$Post.Op.Surgical.Site.Infection1, plot=FALSE)
#rect(xleft=fhist$mids-2,ybottom=0,xright=fhist$mids+2,ytop=fhist$counts/nrow(df),col=gray(0.5))
#rect(xleft=thist$mids-2,ybottom=1-thist$counts/nrow(df),xright=thist$mids+2,ytop=1,col=gray(0.5))
#curve(predict(tinymodel, data.frame(Last_Hematocrit=x), type="response"), col="red", lwd=3, add=TRUE)

dev.new()
library(pROC)
test_prob1 = predict(model1, type = "response")

#cresults=df[!is.na(df$Last_Hematocrit) & !is.na(df$Post_Op_Pulm_Vent_Prolonged) & !is.na(df$Post_Op_Surgical_Site_Infection)
#     & !is.na(df$Chronic_Lung_Disease_Adj) & !is.na(df$Valve), ]$Composite_Readmit_Mort 
test_roc1 = roc(df0$Composite_Readmit_Mort ~ test_prob1, plot = TRUE, print.auc = FALSE,col='red')
text(.55, .45, paste(model1$rank-1,"-predictor model",sep=""),col='red')
text(.55, .4, paste("AUC:",round(test_roc1$auc,4)),col='red')
text(.55, .35, paste("Events:",length(test_roc1$cases)),col='red')
text(.55, .3, paste("Patients:",length(test_roc1$controls)+length(test_roc1$cases)),col='red')
#print(test_roc2)
test_prob2 = predict(model2, type = "response")
test_roc2 = roc(df0$Composite_Readmit_Mort ~ test_prob2, plot = TRUE, print.auc = FALSE,col='blue',add=TRUE)
text(.15, .75, paste(model2$rank-1,"-predictor model (STEP 1)",sep=""),col='blue')
text(.15, .7, paste("AUC:",round(test_roc2$auc,4)),col='blue')
text(.15, .65, paste("Events:",length(test_roc2$cases)),col='blue')
text(.15, .6, paste("Patients:",length(test_roc2$controls)+length(test_roc2$cases)),col='blue')
#print(test_roc1)

test_prob3 = predict(model3, type = "response")
test_roc3 = roc(df0$Composite_Readmit_Mort ~ test_prob3, plot = TRUE, print.auc = FALSE,col='black',add=TRUE)
text(.85, .75, paste(model3$rank-1,"-predictor model",sep=""),col='black')
text(.85, .7, paste("AUC:",round(test_roc3$auc,4)),col='black')
text(.85, .65, paste("Events:",length(test_roc3$cases)),col='black')
text(.85, .6, paste("Patients:",length(test_roc3$controls)+length(test_roc1$cases)),col='black')

title(main = "Initial Patient-specific Predictors, Temperature not Considered",line=2.5)

#print(test_roc3)
#test_prob0 = predict(mylogit, type = "response")
#test_roc0 = roc(df$Composite_Readmit_Mort ~ test_prob0, plot = TRUE, print.auc = TRUE,col='black',add=TRUE)
#print(test_roc0)


arv2<-confint(model2)
arvcomb<-merge(arv2,coef(summary(model2)),by='row.names')
rownames(arvcomb) <- arvcomb$Row.names
arvcomb <- arvcomb[,c(-1)]
colnames(arvcomb)=c("lower2","upper2","v2","stderr2","zval2","pval2")
arv1<-confint(model1)
colnames(arv1)=c("lower1","upper1")
arvcomb12<-merge(arv1,arvcomb,by='row.names',all=T)
rownames(arvcomb12) <- arvcomb12$Row.names
arvcomb12 <- arvcomb12[,c(-1)]
cmod1=coef(summary(model1))
colnames(cmod1)=c("v1","stderr1","zval1","pval1")
arvcomb12<-merge(arvcomb12,cmod1,by='row.names',all=T)
rownames(arvcomb12) <- arvcomb12$Row.names
arvcomb12 <- arvcomb12[,c(-1)]
arv3<-confint(model3)
colnames(arv3)=c("lower3","upper3")
arvcomb13<-merge(arv3,arvcomb12,by='row.names',all=T)
rownames(arvcomb13) <- arvcomb13$Row.names
arvcomb13 <- arvcomb13[,c(-1)]
cmod3=coef(summary(model3))
colnames(cmod3)=c("v3","stderr3","zval3","pval3")
arvcomb13<-merge(arvcomb13,cmod3,by='row.names',all=T)
rownames(arvcomb13) <- arvcomb13$Row.names
#arvcomb13 <- arvcomb13[,c(-1)]
arvcomb13$pval1[is.na(arvcomb13$pval1)]=100
arvcomb13$pval3[is.na(arvcomb13$pval3)]=100
arvcomb13$pval2[is.na(arvcomb13$pval2)]=100
arvcomb13 <- arvcomb13[order(arvcomb13$pval2,arvcomb13$pval3,decreasing=TRUE),]
arvf <- arvcomb13[(nrow(arvcomb13)-model2$rank+1-6):nrow(arvcomb13),]

dev.new()
arvf$ID <- seq.int(nrow(arvf))
alterc=c(2:9,12:13,16:17)
arvf['RF.Last.Hematocrit',alterc] = arvf['RF.Last.Hematocrit',alterc] *15
arvf['RF.Last.Hematocrit','Row.names'] ="Last Hematocrit +15"
arvf['Lowest.Hematocrit',alterc] = arvf['Lowest.Hematocrit',alterc] *15
arvf['Lowest.Hematocrit','Row.names'] ="Lowest Hematocrit +15"

arvf['Total.ICU.Hours',alterc] = arvf['Total.ICU.Hours',alterc] *400
arvf['Total.ICU.Hours','Row.names'] ="Total.ICU.Hours +400"

#arvf['NO2_Post_Means',alterc] = arvf['NO2_Post_Means',alterc] *35
#arvf['NO2_Post_Means','Row.names'] ="NO2_Post_Means +35"
arvf['BMI',alterc] = arvf['BMI',alterc] *30
arvf['BMI','Row.names'] ="BMI +30"

arvf['Total.Postoperative.Ventilation.Hours',alterc] = arvf['Total.Postoperative.Ventilation.Hours',alterc] *200
arvf['Total.Postoperative.Ventilation.Hours','Row.names'] ="Total.Postoperative.Ventilation.Hours +200"

arvf['RF.Last.Creat.Level',alterc] = arvf['RF.Last.Creat.Level',alterc] *10
arvf['RF.Last.Creat.Level','Row.names'] = "Last Creat +10"
labelloc = match("(Intercept)",rownames(arvf))
sp1 <- ggplot() + geom_point(data=arvf, aes(x=v2, y=factor(ID)),color="blue")+geom_errorbar(data=arvf,aes(xmin = lower2, xmax = upper2, y=factor(ID)),color="blue",width=0.2) +
     geom_point(data=arvf, aes(x=v3, y=ID+0.2),color="black")+geom_errorbar(data=arvf,aes(xmin = lower3, xmax = upper3, y=ID+0.2),color="black",width=0.2) +
     geom_point(data=arvf, aes(x=v1, y=ID-0.2),color="red")+geom_errorbar(data=arvf,aes(xmin = lower1, xmax = upper1, y=ID-0.2),color="red",width=0.2)+
     annotate("text",x=4, y=labelloc,label=paste(model2$rank-1,"-predictor model (STEP 1)",sep=""),col='blue')+
    annotate("text",x=4, y=labelloc+.4,label=paste(model3$rank-1,"-predictor (",model2$rank+5," shown)",sep=""),col='black')+
    annotate("text",x=4, y=labelloc-.4,label=paste(model1$rank-1,"-predictor model",sep=""),col='red')
print(sp1+ scale_y_discrete(breaks=arvf$ID, labels=arvf$Row.names)+xlab("Logit Score Risk") + ylab("") +geom_vline(xintercept=0,linetype=3))
write.csv(arvcomb13, "~/Documents/dlnmRev/coefficients_raw_pvals.csv")
write.csv(arvf, "~/Documents/dlnmRev/coefficients_pvals_as_plotted.csv")
}