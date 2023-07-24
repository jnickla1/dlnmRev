library("readxl")
my_data <- read_excel("~/Documents/deidentified_Mar8.xlsx")
df <- my_data[ -c(1,3:5) ]
delcols <-c(30,33,35)
fcols <-c(3:4,6,8:28,37:45,50,52:55,57:83)
df[fcols] <- lapply(df[fcols], factor)
df <- df[ -delcols ]
#colnames(df)
mylogit <- glm(Composite_Readmit_Mort ~ ., data = df, family = "binomial")

library("glmtoolbox")
fit1 <- stepCriterion(mylogit,direction="forward", criterion="p-value")
#10 covariates
fit2 <- stepCriterion(mylogit, criterion="bic")
#5 covariates, disagree on Valve, only RF-Last Hematocrit is non-categorical
model1 <- glm(paste("Composite_Readmit_Mort",fit1$final), data = df, family = "binomial")
summary(model1)
plot(model1, which=1)

dev.new()
model2 <- glm(paste("Composite_Readmit_Mort",fit2$final), data = df, family = "binomial")
summary(model2)
plot(model2, which=1)
