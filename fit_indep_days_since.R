library("readxl")
my_data <- read.csv("~/Documents/dlnmRevData/PostCriteria_withTempPoll_DEIDENTIFIED3a.csv",sep=',',header=TRUE)

df1 <- my_data[,c("Readmit_Day","Composite_Readmit_Mort","DOW_Discharge")]

num <- rep(NA, 30)
denom<- rep(NA, 30)
cur_denom=dim(df1)[1]
for (x in 0:29){
  num[x+1]=sum(df1[,1]==x)
  denom[x+1]=cur_denom
  cur_denom=cur_denom-num[x+1]
}
logfrac=log(num/denom)
#logfrac1=mean(logfrac[1],logfrac[2])
lmlogdf = data.frame(logfrac[3:30],c(2:29))
colnames(lmlogdf) = c("logfrac","day")
lmlogfrac = lm(logfrac ~ day,data=lmlogdf)

nd <- data.frame(day=seq(0,29))
p_conf2 <- predict(lmlogfrac,interval="confidence",newdata=nd)

low1=log(qbeta(0.025,num[1]+num[2]+1, denom[1]+denom[2]-num[1]-num[2]+1))
mid1=log(qbeta(0.5,num[1]+num[2]+1, denom[1]+denom[2]-num[1]-num[2]+1))
high1=log(qbeta(0.975,num[1]+num[2]+1, denom[1]+denom[2]-num[1]-num[2]+1))
pconf <-p_conf2
pconf[1,] <-c(mid1,low1,high1)
pconf[2,] <-c(mid1,low1,high1)

#plotting moved to later

##Calculate just the independent days since function - it is surviorship, so combines in a weird way with pt-level pred



#Compare readmission dow group vs any
passed<- rep(NA, 7)
for (x in 0:7){
  passed[x+1] = sum((df1$DOW_Discharge - x + df1$Readmit_Day) %/% 7 + 1*(x>df1$DOW_Discharge))
}

df1$Readmit_Day[df1$Readmit_Day==30]=NA
df1$DOW_Readmit = (df1$Readmit_Day+ df1$DOW_Discharge)%%7
readmitDOW <- table(df1$DOW_Readmit)
chidf = data.frame(readmitDOW,passed[1:7])
rownames(chidf)=chidf[,1]
chidf = chidf[,-1]
colnames(chidf) = c("event","nonevent")
chisq.test(chidf) #p-value = 0.004812

chidf$frac = chidf$event / (chidf$event +chidf$nonevent)


#CONVERT DATAFRAME TO SURVIVAL LONGFORM

df <- my_data[ -c(1,144) ] #remove first ID col

readmitseason <- table(df[df$Composite_Readmit_Mort==1,"DOW_Month"])
noreadmitseason <- table(df[df$Composite_Readmit_Mort==0,"DOW_Month"])
chidsea = data.frame(readmitseason,noreadmitseason)
chidsea <- chidsea[,-c(1,3)]
chisq.test(chidsea)
chidsea$frac = chidsea[,1] / (chidsea[,1]+chidsea[,2])
#3-6 is blood prods
delcols <-c(3:6, 20, 60, 89, 138:139) #,18,40,42,47,89, 98:127,138:139) #remove temperature initially
fcols <-c(1, 9:12,14,17:20,22:24,28,30:31,33:69,71:73,78,80:81,83:86,89:92,132:142)
df[fcols] <- lapply(df[fcols], factor)
df <- df[ -delcols ]
df <- df[df$Readmit_Day>0,]
# Create an empty list to store replicated data frames
replicated_dfs <- vector("list", nrow(df) )
k=1
for (i in 1:nrow(df)) {
  # Get the number of replications from "num_post_disc"
  num_replications <- df$Readmit_Day[i]
  # Create a data frame with the current row replicated by num_replications
  replicated_row <- df[i, -c(91:120)]
  replicated_row$curTemp=c(-100.01)
  replicated_row$curFupDay=c(-1)
  replicated_row$event=c(0)
  replicated_row$ID=i
  replicated_row$weight_numpd <- 1/30
  replicated_rows <- replicate(num_replications, replicated_row, simplify = FALSE)
  # Add a new column for the follow-up measurement
  
  # Add the follow-up column to each replicated data frame
  for (j in 1:num_replications) {
    replicated_rows[j][[1]]$curTemp <- df[i,(j+97-7)] #1+j if first col is ID
    replicated_rows[j][[1]]$curFupDay<- (j-1)
    if (j==num_replications) {
      replicated_rows[j][[1]]$event=df[i,]$Composite_Readmit_Mort
      replicated_rows[j][[1]]$weight_numpd= replicated_rows[j][[1]]$weight_numpd+ df[i,]$Composite_Readmit_Mort*29/30
    }
  }
  
  # Store the replicated data frames in the list
  replicated_dfs <- c(replicated_dfs, replicated_rows)
}
print("List compiled")
final_df <- as.data.frame(data.table::rbindlist(replicated_dfs))
#final_df <- final0_df[,!(names(final0_df) %in% c("Composite_Readmit_Mort"))]
#remember to remove this second spot that the answer is being stored
print("Dataframe converted")
rm(replicated_rows)
rm(replicated_dfs)

##Construct extra columns for converted dataframe, 3 basis functions for the lag days
final_df$curFupDb0 <- 1*(final_df$curFupDay<2)
final_df$curFupDb1 <- 1*(final_df$curFupDay>=2)
final_df$curFupDb2 <- final_df$curFupDay*(final_df$curFupDb1)
final_df$curDOW= (final_df$DOW_Discharge + final_df$curFupDay) %% 7
final_df[c("curDOW","DOW_Discharge")] <- lapply(final_df[c("curDOW","DOW_Discharge")], factor)

if(!exists("fit2")){
source("~/Documents/dlnmRev/glm_wo_temperature_new.R")
}
## ?? LOSS for each row using cross entropy loss?? 
modelspec = paste("final_df$event",fit2$final,"+ curDOW + curFupDb0 + curFupDb2")
cblogit <- glm(modelspec, data=final_df, family = binomial(),na.action=na.exclude, weights=final_df$weight_numpd)
summary(cblogit)
## ??? Improve these expected probs with GLM binomial, additional predictors DOW and season


dev.new() ##post-discharge day
plot(x=0:29,y=logfrac, xlab ="Post-Discharge Day",ylim=c(-8, -4.5))
polygon(c(0:29,rev(0:29)),c(pconf[,"lwr"],rev(pconf[,"upr"])),col=rgb(.2, 0, 1,0.3),lty=0)
lines(0:29,pconf[,c("fit")],col=1,lty=1)
title(main = "Fit of Risk on Post-Discharge Day")
text(4, -7, "Indep. 2.5-97.5% Conf. Interval",col=rgb(.2, 0, 1,0.6), pos=4)
text(4, -6.8, "Indep. Linear Fit for Days 2:29",col=1, pos=4)

avgplt = mean(logfrac)
ncbfit=data.frame(fit=double(30),lwr=double(30),upr=double(30))
ncbfit$fit = 0:29 * cblogit$coefficients["curFupDb2"]
ncbfit$lwr = ncbfit$fit - 2 * sqrt((0:29 -15.5)^2 * summary(cblogit)$coefficients["curFupDb2", 2]^2 
                                   + summary(lmlogfrac)$coefficients["(Intercept)", 2]^2 )
                                   #+ summary(cblogit)$coefficients["(Intercept)", 2]^2 )
ncbfit$upr = ncbfit$fit + 2 * sqrt((0:29 -15.5)^2 * summary(cblogit)$coefficients["curFupDb2", 2]^2 
                                   + summary(lmlogfrac)$coefficients["(Intercept)", 2]^2 )
                                   #+ summary(cblogit)$coefficients["(Intercept)", 2]^2 ) 
initconst = cblogit$coefficients["curFupDb0"]
initcombste = 2*sqrt( summary(cblogit)$coefficients["(Intercept)", 2]^2 + summary(cblogit)$coefficients["curFupDb2", 2]^2)
ncbfit[1,] <-c(initconst,initconst-initcombste,initconst+initcombste)
ncbfit[2,] <-c(initconst,initconst-initcombste,initconst+initcombste)
avgncb = mean(ncbfit[,1])
ncbfit <- ncbfit -avgncb + avgplt
#polygon(c(0:29,rev(0:29)),c(pconf$lwr,rev(pconf$upr)),col=rgb(.2, 0, 1,0.3),lty=0)
lines(0:29,ncbfit[,"fit"],col="red",lty=2,lwd = 2)
polygon(c(0:29,rev(0:29)),c(ncbfit[,"lwr"],rev(ncbfit[,"upr"])),col=rgb(1, 0, 0,0.2),lty=0)

text(4, -7.8, "Combined (no temp.) 2.5-97.5% CI",col=rgb(1, 0, 0,0.4), pos=4)
text(4, -7.4, "Combined (no temp.) Fit: ",col="red", pos=4)
text(4, -7.6, "     8-var (pt), DOW, Post-Disch. Day",col="red", pos=4)


dev.new() #current day

plot(0:6, log(chidf$frac), xlim=c(0,6), ylim=c(-7.5,-4.5), xlab="DOW", ylab="logfrac")
arrows(x0=0:6, y0=log(qbeta(0.025,chidf$event+1,chidf$nonevent+1)), x1=0:6, 
       y1=log(qbeta(0.975,chidf$event+1,chidf$nonevent+1)), code=3, angle=90, length=0.2, col="black", lwd=2)
#arrows(x0=x, y0=y-3, x1=x, y1=y+3, code=3, angle=90, length=0.5, col="blue", lwd=2)
avgdowplt=mean(log(chidf$frac))
dowfit = data.frame(fit=double(7),lwr=double(7),upr=double(7))
dowfit[1,"fit"] = 0
dowfit[2:7,"fit"] = cblogit$coefficients[10:15]
dowfit[2:7,"lwr"] = cblogit$coefficients[10:15] - 2 * summary(cblogit)$coefficients[10:15, 2]
dowfit[1,"lwr"] = -2*sqrt(mean(summary(cblogit)$coefficients[10:15, 2]^2))

dowfit[2:7,"upr"] = cblogit$coefficients[10:15] + 2 * summary(cblogit)$coefficients[10:15, 2]
dowfit[1,"upr"] = -dowfit[1,"lwr"]

avgdowncb = mean(dowfit[,1])
dowfit <- dowfit -avgdowncb + avgdowplt

points((0:6+0.2), dowfit[,1],col="red")
arrows(x0=(0:6+0.2), y0=dowfit[,2], x1=(0:6+0.2), y1=dowfit[,3], code=3, angle=90, length=0.2, col="red", lty=2,lwd = 2)
title(main = "Fit of Risk on Current Day of the Week (DOW)")
text(1, -6.8, "Indep. Bayesian CIs",col=1, pos=4)
text(1, -7.2, "Combined (no temp.) 2.5-97.5% CI",col="red", pos=4)
text(1, -7.4, "     8-var (pt), DOW, Post-Disch. Day",col="red", pos=4)


model2 <- glm(paste("Composite_Readmit_Mort",fit2$final), data = df0, family = "binomial",na.action=na.exclude)
arv2<-confint(model2)
arvcomb<-merge(arv2,coef(summary(model2)),by='row.names')
rownames(arvcomb) <- arvcomb$Row.names
arvcomb <- arvcomb[,c(-1)]
colnames(arvcomb)=c("lower2","upper2","v2","stderr2","zval2","pval2")
arv1<-confint(cblogit)
colnames(arv1)=c("lower1","upper1")
arvcomb12<-merge(arv1,arvcomb,by='row.names')
rownames(arvcomb12) <- arvcomb12$Row.names
arvcomb12 <- arvcomb12[,c(-1)]
cmod1=coef(summary(cblogit))
colnames(cmod1)=c("v1","stderr1","zval1","pval1")
arvcomb12<-merge(arvcomb12,cmod1,by='row.names')
rownames(arvcomb12) <- arvcomb12$Row.names

arvcomb12$pval2[is.na(arvcomb12$pval2)]=100
arvcomb12 <- arvcomb12[order(arvcomb12$pval2,decreasing=TRUE),]
arvf <- arvcomb12

dev.new() #adjusted other preds
arvf$ID <- seq.int(nrow(arvf))
alterc=c(2:7,10:11)
arvf['RF.Last.Hematocrit',alterc] = arvf['RF.Last.Hematocrit',alterc] *15
arvf['RF.Last.Hematocrit','Row.names'] ="Last Hematocrit +15"

arvf['NO2_Post_Means',alterc] = arvf['NO2_Post_Means',alterc] *35
arvf['NO2_Post_Means','Row.names'] ="NO2_Post_Means +35"
arvf['BMI',alterc] = arvf['BMI',alterc] *30
arvf['BMI','Row.names'] ="BMI +30"

arvf['Total.ICU.Hours',alterc] = arvf['Total.ICU.Hours',alterc] *400
arvf['Total.ICU.Hours','Row.names'] ="Total.ICU.Hours +400"
library(ggplot2)
labelloc = match("(Intercept)",rownames(arvf))
sp1 <- ggplot() + geom_point(data=arvf, aes(x=v2, y=factor(ID)),color="blue")+
  geom_errorbar(data=arvf,aes(xmin = lower2, xmax = upper2, y=factor(ID)),color="blue",width=0.2) +
  geom_point(data=arvf, aes(x=v1, y=ID-0.2),color="red")+geom_errorbar(data=arvf,aes(xmin = lower1, 
                                xmax = upper1, y=ID-0.2),width=0.5,color="red",lty=2,lwd = .6) +
  annotate("text",x=4, y=labelloc,label=paste(model2$rank-1,"-predictor model",sep=""),col='blue')+
  annotate("text",x=4, y=labelloc-.4,label="8-pred + DOW + post-Disch Day",col='red')
print(sp1+ scale_y_discrete(breaks=arvf$ID, labels=arvf$Row.names)+xlab("Logit Score Risk") + ylab("") +geom_vline(xintercept=0,linetype=3))

cblogit_preds <- predict(cblogit,se.fit=TRUE)
