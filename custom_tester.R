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


source("~/Documents/dlnmRev/custom_crossbasis.R")
environment(custom_crossbasis) <- asNamespace('dlnm')
cb1 <- crossbasis(dftempall[2:31],lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)))
cb2 <- custom_crossbasis(dftempall[2:31],lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)))
