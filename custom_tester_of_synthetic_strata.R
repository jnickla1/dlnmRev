library("readxl")
library("dlnm")
library("splines")
library("glmtoolbox")
library("dplyr")
library("data.table")

tempers <- read.csv("~/Documents/dlnmRevData/sampled_flat_data.csv",sep=',',header=TRUE)
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
final_df <- as.data.frame(data.table::rbindlist(replicated_dfs))
print("Dataframe converted")
rm(replicated_rows)
rm(replicated_dfs)


varper <- c(10,30,80)
lag <- 29
lagnk <- 2
#argvar <- list(fun="bs",degree=2,knots=quantile(dftempall[1:30], varper/100,na.rm=T))
argvar <- list(fun="strata",breaks=quantile(dftempall[1:30], varper/100,na.rm=T))
arglag2 <-list(fun="strata", breaks = c(6), intercept=TRUE)

source("~/Documents/dlnmRev/custom_crossbasis.R")
environment(custom_crossbasis) <- asNamespace('dlnm')
cb2 <- custom_crossbasis(final_df$curTemp,lag=lag,argvar=argvar, arglag=arglag2,group=final_df$ID, pslags=final_df$curFup)
print("Crossbasis formed")
cblogit <- glm(final_df$Outcome ~ cb2, family = "binomial", weights=final_df$weight_numpd)
print("Logit fitted, with NA values")
print(sum(is.na(cblogit$coefficients)))
print("reduce degrees of freedom to cut back on NA values.")
#cb1 <- crossbasis(dftempall[2:31],lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)))
#cb2 <- custom_crossbasis(dftempall[2:31],lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)))
source("~/Documents/dlnmRev/custom_crosspred.R")
environment(custom_crosspred) <- asNamespace('dlnm')
cp <- custom_crosspred(cb2, cblogit,from=-2, to=90, lag=c(0,30), by=1)
crosspred.image.real(cp$lagRRfit,levels=seq(-5.5,5.5,0.5),title="lag")
crosspred.image.real(cp$dayRRfit,title="day")
points(x=tempers[,"tmean"],y=(tempers[,"num_post_disc"]+1))

collapsedD= vector(length = 30)
for (i in 1:30) { 
  a = hist(dftempall[,i],breaks=seq(-2.5,90.5,1),plot=FALSE)$density 
  collapsedD[i] = sum(a * log(cp$dayRRfit[,i]))
}
plot(collapsedD)

for (x in 2:30) { 
  tempers[,x] <- ifelse( (is.na(tempers[,x]) & !(tempers[,x-1] >= 1000)), 1000, tempers[,x])
  }
p_res=colSums(tempers[,1:30]==1000, na.rm=TRUE)/colSums(tempers[,1:30]>-10, na.rm=TRUE)
plot(p_res, type="l")
points(collapsedD/10+exp(cblogit$coefficients[1]))
arglag$intercept <- TRUE
basislag <- do.call("onebasis", modifyList(arglag, list(x = seq(0,29))))
ob <- glm(p_res ~ basislag+ 0)
op <- crosspred(basislag,ob,from=1, to=30, by=1)
plot(p_res, type="l")
points(op)
points((collapsedD-collapsedD[1])/10)
