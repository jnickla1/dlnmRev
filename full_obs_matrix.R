library("readxl")
library("dlnm")
library("splines")
library("glmtoolbox")
library("dplyr")
library("data.table")

my_data <- read_excel("~/Documents/deidentified_Jul28.xlsx")
delcols <-c(30,33,35)
fcols <-c(3:4,6,8:28,37:45,50,52:55,57:75)
df2 <- my_data[ -c(1,3:7) ]
df2[fcols] <- lapply(df2[fcols], factor)
df2 <- df2[ -delcols ]
df2$dow <-as.factor(my_data$Disch_Day_of_Week)
df2$month <- my_data$Disch_Month
df2$ID <- seq.int(nrow(df2))-1


tempers <- read.csv("~/Documents/Mar8_envs_tem.csv",sep=',',header=TRUE)
tempers$num_post_disc <- rowSums(!is.na(tempers)) - 1
tempers$tmean <- rowMeans(tempers[1:30],na.rm=TRUE)
dftempall <- merge(tempers,df2,by="ID")
dftempall <- dftempall[!is.na(dftempall$t0),]


# Create an empty list to store replicated data frames
replicated_dfs <- vector("list", nrow(dftempall))
for (i in 1:nrow(dftempall)) {
  # Get the number of replications from "num_post_disc"
  num_replications <- dftempall$num_post_disc[i]
  # Create a data frame with the current row replicated by num_replications
  replicated_row <- dftempall[i, -c(2:31)]
  replicated_row$curTemp=c(-100.01)
  replicated_row$curFup=c(0)
  replicated_row$event=c(0)
  replicated_row$weight_numpd <- 1/num_replications
  replicated_rows <- replicate(num_replications, replicated_row, simplify = FALSE)
  # Add a new column for the follow-up measurement
  
  # Add the follow-up column to each replicated data frame
  for (j in 1:num_replications) {
    replicated_rows[j][[1]]$curTemp <- dftempall[i,(1+j)]
    replicated_rows[j][[1]]$curFup<- (j)
    if (j==num_replications) {
      replicated_rows[j][[1]]$event=dftempall[i,]$Composite_Readmit_Mort
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

varper <- c(10,30,60,80)
lag <- 29
lagnk <- 2
argvar <- list(fun="bs",degree=2,knots=quantile(dftempall[2:31], varper/100,na.rm=T))

source("~/Documents/dlnmRev/custom_crossbasis.R")
environment(custom_crossbasis) <- asNamespace('dlnm')
cb2 <- custom_crossbasis(final_df$curTemp,lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk)),group=final_df$ID, pslags=final_df$curFup)
print("Crossbasis formed")
cblogit <- glm(final_df$event ~ cb2, family = "binomial", weights=final_df$weight_numpd)
print("Logit fitted, with NA values")
sum(is.na(cblogit$coefficients))
print("reduce degrees of freedom to cut back on NA values.")