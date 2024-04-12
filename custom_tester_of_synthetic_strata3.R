library("readxl")
library("dlnm")
library("splines")
library("glmtoolbox")
library("dplyr")
library("data.table")
library(R.utils)

tempers <- read.csv("~/Documents/dlnmRevData/sampled_corner_data.csv",sep=',',header=TRUE)
tempers$num_post_disc <- rowSums(!is.na(tempers)) - 2
tempers$tmean <- rowMeans(tempers[1:30],na.rm=TRUE)
dftempall <- tempers

# Create an empty list to store replicated data frames
replicated_dfs <- vector("list", nrow(dftempall))
for (i in 1:nrow(dftempall)) {
  # Get the number of replications from "num_post_disc"
  num_replications <- dftempall$num_post_disc[i]
  # Create a data frame with the current row replicated by num_replications
  replicated_row <- dftempall[i, -c(1:30)]
  replicated_row$curTemp=c(-100.01)
  replicated_row$curFup=c(0)
  replicated_row$event=c(0)
  replicated_row$weight_numpd <- 1/num_replications
  replicated_rows <- replicate(num_replications, replicated_row, simplify = FALSE)
  # Add a new column for the follow-up measurement
  
  # Add the follow-up column to each replicated data frame
  for (j in 1:num_replications) {
    replicated_rows[j][[1]]$curTemp <- dftempall[i,(j)] #1+j if first col is ID
    replicated_rows[j][[1]]$curFup<- (j)
    if (j==num_replications) {
      replicated_rows[j][[1]]$event=dftempall[i,]$Outcome
    }
  }
  
  # Store the replicated data frames in the list
  replicated_dfs <- c(replicated_dfs, replicated_rows)
}
print("List compiled")
final0_df <- as.data.frame(data.table::rbindlist(replicated_dfs))
final_df <- final0_df[,!(names(final0_df) %in% c("Outcome"))]
print("Dataframe converted")
rm(replicated_rows)
rm(replicated_dfs)


varper <- c(10,30,80)
lag <- 29
lagnk <- 2
#argvar <- list(fun="bs",degree=2,knots=quantile(dftempall[1:30], varper/100,na.rm=T))
argvar <- list(fun="strata",breaks=quantile(dftempall[1:30], varper/100,na.rm=T), intercept=TRUE)
arglag2 <- list(fun="bs",degree=2,knots=c( 1,5,10), intercept=FALSE)

vb = quantile(dftempall[1:30], varper/100,na.rm=T)
  #list(fun="strata", breaks = c(6), intercept=TRUE)
source("~/Documents/dlnmRev/custom_crossbasis_scatter.R")
environment(custom_crossbasis_scatter) <- asNamespace('dlnm')
totals = custom_crossbasis_scatter(final_df$curTemp,lag=lag,varbreaks=vb,group=final_df$ID, pslags=final_df$curFup)



argvar2 = list(fun="bs",degree=2,knots=quantile(dftempall[1:30], varper/100,na.rm=T))
xminmax0 = c(min(final_df$curTemp),maxfinal_df$curTemp),0.1)
#source("~/Documents/dlnmRev/custom_crossbasis.R")
#environment(custom_crossbasis) <- asNamespace('dlnm')
#cb2 <- custom_crossbasis(final_df$curTemp,lag=lag,argvar=argvar, arglag=arglag2,group=final_df$ID, pslags=final_df$curFup)

source("~/Documents/dlnmRev/custom_crossbasis_gen.R")
environment(custom_crossbasis_gen) <- asNamespace('dlnm')
imprints = custom_crossbasis_gen(xminmax0,lag=lag, varbreaks = vb, 
                                 argvar=argvar2, arglag=arglag2)
cb2a <-


print("Crossbasis formed")

A = matrix(
     c(   1 ,    0  ,   0  ,   0   ,  0    , 0   ,  0    , 0, 0,0,0,0,0,0,0,0,
         -1    , 1 ,    0   ,  0   ,  0   ,  0  ,   0   ,  0, 0,0,0,0,0,0,0,0,
         -1  ,   0  ,   1  ,   0   ,  0   ,  0   ,  0   ,  0, 0,0,0,0,0,0,0,0,
           1  ,  -1  ,  -1   ,  1   ,  0   ,  0   ,  0  ,   0, 0,0,0,0,0,0,0,0,
          -1  ,   0  ,   0 ,   0    , 1  ,   0  ,   0   ,  0, 0,0,0,0,0,0,0,0,
           1  ,  -1 ,    0  ,   0  ,  -1  ,   1   ,  0   ,  0, 0,0,0,0,0,0,0,0,
          1   ,  0   , -1   ,  0 ,   -1  ,   0   ,  1   ,  0, 0,0,0,0,0,0,0,0,
          -1  ,   1 ,    1  ,  -1  ,   1  ,  -1  ,  -1   ,  1, 0,0,0,0,0,0,0,0,
         -1  ,   0  ,   0 ,   0  , 0,0,0,0,  1  ,   0  ,   0   ,  0, 0,0,0,0,
         1  ,  -1 ,    0  ,   0 , 0,0,0,0,  -1  ,   1   ,  0   ,  0, 0,0,0,0,
         1   ,  0   , -1   ,  0 , 0,0,0,0,  -1  ,   0   ,  1   ,  0, 0,0,0,0,
         -1  ,   1 ,    1  ,  -1  , 0,0,0,0,   1  ,  -1  ,  -1   ,  1, 0,0,0,0,       
         -1  ,   0  ,   0 ,   0  , 0,0,0,0,0,0,0,0,  1  ,   0  ,   0   ,  0, 
         1  ,  -1 ,    0  ,   0 , 0,0,0,0,0,0,0,0,  -1  ,   1   ,  0   ,  0, 
         1   ,  0   , -1   ,  0 , 0,0,0,0,0,0,0,0,  -1  ,   0   ,  1   ,  0, 
         -1  ,   1 ,    1  ,  -1  , 0,0,0,0,0,0,0,0,   1  ,  -1  ,  -1   ,  1
         ), 
       nrow = 16, ncol = 16,byrow = TRUE )
hits <- cb2 %*% A
whits <- sweep(hits,1,final_df$weight_numpd,"*")
eventH <- hits[(final_df$event==1),]
uneventH <- hits[(final_df$event==0),]
print(sum(final_df$event))
print(colSums(eventH))
print(colSums(uneventH))

# cblogit <- glm(final_df$event ~ cb2, family = "binomial", weights=final_df$weight_numpd)
# print("Logit fitted, with NA values")
# print(sum(is.na(cblogit$coefficients)))
# #cblogit$coefficients[2:17] %*% A
# coefs0 <- cblogit$coefficients[2:17]
# coefs0[is.na(coefs0)] <- 0
# coefs0 %*% A

cbm <- cb2
rankifremoved <- sapply(1:ncol(cbm), function (x) qr(cbm[,-x])$rank)
remcols <- c()
while(max(rankifremoved) - min(rankifremoved) > 0){
  removing_now <- which.max(rankifremoved)
  cbm <- cbm[,-(removing_now)]
  print(dim(cbm))
  rankifremoved <- sapply(1:ncol(cbm), function (x) qr(cbm[,-x])$rank)
  remcols <- append(remcols, removing_now)}


Am = matrix(
  c(        1  ,   0  ,  0,   0  ,   0   ,  0  ,0 ,  0 , 0,0,0,0, 0,0,0,0,
            -1 ,    1  , 0,   0 ,    0  ,   0  ,0 ,  0 , 0,0,0,0, 0,0,0,0, 
            0  ,  -1   , 0,  1  ,   0   ,  0   ,0 ,  0 ,  0,0,0,0, 0,0,0,0,
            -1 ,    0  , 0,  0 ,    1  ,   0   ,0 ,   0,  0,0,0,0, 0,0,0,0,
            1  ,  -1   , 0,  0  ,  -1  ,   1   ,0,  0,  0,0,0,0, 0,0,0,0, 
            0  ,   1   , 0, -1  ,   0  ,  -1    ,0,  1 , 0,0,0,0, 0,0,0,0,
            -1 ,    0  , 0,  0 , 0,0,0,0,   1  ,   0   ,0 ,   0,  0,0,0,0, 
            1  ,  -1   , 0,  0  ,0,0,0,0,  -1  ,   1   ,0,  0,  0,0,0,0, 
            0  ,   1   , 0, -1  ,0,0,0,0,   0  ,  -1    ,0,  1 , 0,0,0,0, 
            -1 ,    0  , 0,  0 , 0,0,0,0, 0,0,0,0,   1  ,   0   ,0 ,   0,  
            1  ,  -1   , 0,  0  , 0,0,0,0, 0,0,0,0,  -1  ,   1   ,0,  0,  
            0  ,   1   , 0, -1  , 0,0,0,0, 0,0,0,0,  0  ,  -1    ,0,  1 
  ), 
  nrow = 12, ncol = 16,byrow = TRUE )

hitsm <- cbm %*% Am

cblogit <- glm(final_df$event ~ cbm, family = "binomial", weights=final_df$weight_numpd)
print("Logit fitted, with NA values")
print(sum(is.na(cblogit$coefficients)))
#cblogit$coefficients[2:17] %*% A
coefs0 <- cblogit$coefficients[2:length(cblogit$coefficients)]
coefs0[is.na(coefs0)] <- 0
coefs0 %*% Am

print("reduce degrees of freedom to cut back on NA values.")
#cb1 <- crossbasis(dftempall[2:31],lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)))
#cb2 <- custom_crossbasis(dftempall[2:31],lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)))
source("~/Documents/dlnmRev/custom_crosspred.R")
environment(custom_crosspred) <- asNamespace('dlnm')
cp <- custom_crosspred(cbm, cblogit,from=-2, to=90, lag=c(0,30), by=1)
crosspred.image.real(cp$lagRRfit,levels=seq(-5.5,5.5,0.5),title="lag")
crosspred.image.real(cp$dayRRfit,title="day")
points(x=tempers[,"tmean"],y=(tempers[,"num_post_disc"]+1))
