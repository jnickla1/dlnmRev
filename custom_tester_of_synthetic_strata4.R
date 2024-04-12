library("readxl")
library("dlnm")
library("splines")
library("glmtoolbox")
library("dplyr")


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


varper <- c(10,20,30,50,70,80,90)
lag <- 29
lagnk <- 2
#argvar <- list(fun="bs",degree=2,knots=quantile(dftempall[1:30], varper/100,na.rm=T))
#argvar <- list(fun="strata",breaks=quantile(dftempall[1:30], varper/100,na.rm=T), intercept=TRUE)
arglag2 <- list(fun="bs",degree=2,knots=c( 0,1,5), intercept=FALSE)

vb = quantile(dftempall[1:30], varper/100,na.rm=T)
  #list(fun="strata", breaks = c(6), intercept=TRUE)
source("~/Documents/dlnmRev/custom_crossbasis_scatter.R")
environment(custom_crossbasis_scatter) <- asNamespace('dlnm')
totals = custom_crossbasis_scatter(final_df$curTemp,lag=lag,varbreaks=vb,
                                   group=final_df$ID, pslags=final_df$curFup)

varper2 <- c(10,30,50,80)
argvar2 = list(fun="bs",degree=2,knots=quantile(dftempall[1:30], varper2/100,na.rm=T))
xminmax0 = c(min(final_df$curTemp),max(final_df$curTemp),0.1)
#source("~/Documents/dlnmRev/custom_crossbasis.R")
#environment(custom_crossbasis) <- asNamespace('dlnm')
#cb2 <- custom_crossbasis(final_df$curTemp,lag=lag,argvar=argvar, arglag=arglag2,group=final_df$ID, pslags=final_df$curFup)

source("~/Documents/dlnmRev/custom_crossbasis_gen.R")
environment(custom_crossbasis_gen) <- asNamespace('dlnm')
imprints = custom_crossbasis_gen(xminmax0,lag=lag, varbreaks = vb, 
                                 argvar=argvar2, arglag=arglag2)

df_recomb_list <- data.frame(c(totals))
lagdiml <- length(arglag2$knots)+2
for (var in seq(dim(imprints)[1])){
  v_here = (var-1) %/% (lagdiml*lagdiml)
  l_here =((var-1) %% (lagdiml*lagdiml)) %/% lagdiml
  d_here = (var-1) %% (lagdiml)
  colname = paste0("v",v_here,".l",l_here,".d",d_here)
  df_recomb_list[colname] <- c(imprints[var,,,])
}

print("Crossbasis formed")


fitlm <- lm(c.totals. ~ ., data = df_recomb_list)
print("Interpolation 1 fitted")
#print(sum(is.na(fitlm$coefficients)))
coefs0 <- fitlm$coefficients[2:length(fitlm$coefficients)]
coefs0[is.na(coefs0)] <- 0
fittedarr = array(rep(fitlm$coefficients[1],length(totals)), dim = dim(totals))
for (var in seq(dim(imprints)[1])){
  fittedarr = fittedarr + imprints[var,,,]*coefs0[var]
}


#plot(totals[2,,],breaks=seq(55,68,1))
#plot(fittedarr[2,,],breaks=seq(55,68,1))
#plot(totals[,,10],breaks=seq(10,80,5))
#plot(fittedarr[,,10],breaks=seq(10,80,5))

events_df = final_df[final_df$event == 1, ]
totals_events = custom_crossbasis_scatter(events_df$curTemp,
                  lag=lag,varbreaks=vb,group=events_df$ID, pslags=events_df$curFup)
df_recomb_list$c.totals. = c(totals_events)
fitlm2 <- lm(c.totals. ~ ., data = df_recomb_list)
print("Interpolation 2 fitted")
coefsE <- fitlm2$coefficients[2:length(fitlm2$coefficients)]
coefsE[is.na(coefsE)] <- 0
fittedE = array(rep(fitlm2$coefficients[1],length(totals_events)), dim = dim(totals_events))
for (var in seq(dim(imprints)[1])){
  fittedE = fittedE + imprints[var,,,]*coefsE[var]
}

#plot(totals_events[4,,],breaks=seq(1,3,.25))
#plot(fittedE[4,,],breaks=seq(1,3,.25))
