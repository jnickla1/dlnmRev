## ## Re-run the glm(start = ...) to see where changes / improvements are made
library("readxl")
library("dlnm")
library("splines")

if(!(exists("df3")&&exists("temp2_shifts_ps"))){
  source("~/Documents/dlnmRev/temperature_incorporation.R")
}


#logknots(lag,lagnk)
arglag2 <- list(fun="bs",degree=2,knots=c( 0,2,4,6,12), intercept=FALSE)
source("~/Documents/dlnmRev/custom_crossbasis_gen.R")
environment(custom_crossbasis_gen) <- asNamespace('dlnm')
imprints = custom_crossbasis_gen(xminmax0,lag=lag, varbreaks = vbf, 
                                 argvar=argvar0, arglag=arglag2)
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


#dev.new()
#plot(fittedarr[2,,],breaks=seq(90,120,2))
#dev.new()
#plot(totals[2,,],breaks=seq(90,120,2))


df_recomb_list$c.totals. = c(totals_events)
fitlm2 <- lm(c.totals. ~ ., data = df_recomb_list)
print("Interpolation 2 fitted")



#events_df = final_df[final_df$Composite_Readmit_Mort == 1, ]

coefsE <- fitlm2$coefficients[2:length(fitlm2$coefficients)]
coefsE[is.na(coefsE)] <- 0
fittedE = array(rep(fitlm2$coefficients[1],length(totals_events)), dim = dim(totals_events))
for (var in seq(dim(imprints)[1])){
  fittedE = fittedE + imprints[var,,,]*coefsE[var]
}
#dev.new()
#plot(fittedE[2,,],breaks=seq(0,8,1))
#dev.new()
#plot(totals_events[2,,],breaks=seq(0,8,1))

temp_shifts_cent22 = array(rep(0,length(totals)), dim = dim(totals_counts))
temp_shifts_high22 = array(rep(0,length(totals)), dim = dim(totals_counts))
temp_shifts_low22 = array(rep(0,length(totals)), dim = dim(totals_counts))

for(i in 1:dim(totals_counts)[1]) {
  temp_shifts_cent22[i,,] <- temp2_shifts_cent[i,,] %x% matrix(1, 2, 2)
  temp_shifts_high22[i,,] <- temp2_shifts_high[i,,] %x% matrix(1, 2, 2)
  temp_shifts_low22[i,,] <- temp2_shifts_low[i,,] %x% matrix(1, 2, 2)
}


#dev.new()
#plot(temp_shifts_cent[2,,],breaks=seq(-0.2,0.2,0.1))
df_recomb_list$c.totals. = c(temp_shifts_cent22)
shiftslm <- lm(c.totals. ~ ., data = df_recomb_list)
print("Interpolation to shifts fitted")
coefsS <- shiftslm$coefficients[2:length(shiftslm$coefficients)]
coefsS[is.na(coefsS)] <- 0
fittedS = array(rep(shiftslm$coefficients[1],length(totals_counts)), dim = dim(totals_counts))
for (var in seq(dim(imprints)[1])){
  fittedS = fittedS + imprints[var,,,]*coefsS[var]}




df_recomb_list$c.totals. = c(temp_shifts_high22)
shiftslmhi <- lm(c.totals. ~ ., data = df_recomb_list)
print("Interpolation to shifts high fitted")
coefsShi <- shiftslmhi$coefficients[2:length(shiftslmhi$coefficients)]
coefsShi[is.na(coefsShi)] <- 0
fittedShi = array(rep(shiftslmhi$coefficients[1],length(totals_counts)), dim = dim(totals_counts))
for (var in seq(dim(imprints)[1])){
  fittedShi = fittedShi + imprints[var,,,]*coefsShi[var]}


df_recomb_list$c.totals. = c(temp_shifts_low22)
shiftslmlo <- lm(c.totals. ~ ., data = df_recomb_list)
print("Interpolation to shifts low fitted")
coefsSlo <- shiftslmlo$coefficients[2:length(shiftslmlo$coefficients)]
coefsSlo[is.na(coefsSlo)] <- 0
fittedSlo = array(rep(shiftslmlo$coefficients[1],length(totals_counts)), dim = dim(totals_counts))
for (var in seq(dim(imprints)[1])){
  fittedSlo = fittedSlo + imprints[var,,,]*coefsSlo[var]}


source("~/Documents/dlnmRev/custom_crossbasis_point.R")
environment(custom_crossbasis_point) <- asNamespace('dlnm')

coeffs_basis <- shiftslm$coefficients[-1]
coeffs_basis[is.na(coeffs_basis)] <- 0
len_to_calc=7
all_weights = custom_crossbasis_point(rep(298,7),c(1,2,3,4,5,6,7), c(1,2,3,4,5,6,7), xminmax0,lag=lag, argvar=argvar0, arglag=arglag2)
point_values = array(rep(shiftslmlo$coefficients[1],len_to_calc)) + rowSums(matrix(rep(coeffs_basis, times=len_to_calc), nrow=len_to_calc, byrow=TRUE) * all_weights)

point_values
diag(fittedS[7,2:8,2:8])

#plot(temp2_shifts_cent[8,,])
#plot(fittedS[8,,])


##PLOTTING: slices
#col 1 - temperature on discharge day 0
#col 2 - temperature on discharge day +1
# y axis: lag since = current days since discharge
tempx=seq(xminmax0[1],xminmax0[2]+0.3,0.3)
crmat <- matrix(, nrow = length(tempx), ncol =lag+1 )
for (day in seq(0,lag)){
  all_weights_day = custom_crossbasis_point(tempx,rep(0,length(tempx)), rep(day,length(tempx)), xminmax0,lag=lag, argvar=argvar0, arglag=arglag2)
  len_to_calcD=dim(all_weights_day)[1]
  point_values = array(rep(shiftslmlo$coefficients[1],len_to_calcD)) + 
    rowSums(matrix(rep(coeffs_basis, times=len_to_calcD), nrow=len_to_calcD, byrow=TRUE) * all_weights_day)
  crmat[,day+1]=point_values
}
#par(mar=c(5, 4, 4, 4))
#plot(crmat)

noeff=0
levels <- pretty(crmat[,c(-1,-2)], 20)
col1 <- colorRampPalette(c("blue", "white"))
col2 <- colorRampPalette(c("white", "red"))
maxlevels <- max((sum(levels < noeff)),(sum(levels > noeff)))
minlevels <- min((sum(levels < noeff)),(sum(levels > noeff)))
col <- c(col1(maxlevels)[(maxlevels-minlevels):maxlevels], col2(maxlevels))
dev.new()
filled.contour(x = tempx -273.15, y = seq(2, lag), 
               z = crmat[,c(-1,-2)], col = col, levels = levels, 
               plot.title = title(main = "Smoothed Binned Effect (logit) \n of Today's Temperature (Lag=0d)",
                                  xlab = "Temperature °C", ylab = "Days post-Discharge"))

overlay_points <- function(x, y, top_value, bottom_value, xspacing, yspacing) {
  for (i in seq(1, length(x), by = xspacing)) {
    for (j in seq(1, length(y), by = yspacing)) {
      if (top_value[i, j] > 0 && bottom_value[i, j] < 0) {
        points(x[i], y[j], pch = 16, col = "grey")
      }}}}

hrmat <- matrix(, nrow = length(tempx), ncol =lag+1 )
lrmat <- matrix(, nrow = length(tempx), ncol =lag+1 )
coeffs_basis <- shiftslmhi$coefficients[-1]
coeffs_basis[is.na(coeffs_basis)] <- 0
for (day in seq(0,lag)){
  all_weights_day = custom_crossbasis_point(tempx,rep(0,length(tempx)), rep(day,length(tempx)), xminmax0,lag=lag, argvar=argvar0, arglag=arglag2)
  len_to_calcD=dim(all_weights_day)[1]
  point_values = array(rep(shiftslmlo$coefficients[1],len_to_calcD)) + 
    rowSums(matrix(rep(coeffs_basis, times=len_to_calcD), nrow=len_to_calcD, byrow=TRUE) * all_weights_day)
  hrmat[,day+1]=point_values
}
coeffs_basis <- shiftslmlo$coefficients[-1]
coeffs_basis[is.na(coeffs_basis)] <- 0
for (day in seq(0,lag)){
  all_weights_day = custom_crossbasis_point(tempx,rep(0,length(tempx)), rep(day,length(tempx)), xminmax0,lag=lag, argvar=argvar0, arglag=arglag2)
  len_to_calcD=dim(all_weights_day)[1]
  point_values = array(rep(shiftslmlo$coefficients[1],len_to_calcD)) + 
    rowSums(matrix(rep(coeffs_basis, times=len_to_calcD), nrow=len_to_calcD, byrow=TRUE) * all_weights_day)
  lrmat[,day+1]=point_values
}


#overlay_points(x = tempx -273.15, y = seq(2, lag), hrmat,lrmat, 2,2)
filled.contour(x = tempx -273.15, y = seq(2, 11), 
               z = crmat[,3:12], col = col, levels = levels, 
               plot.title = title(main = "Smoothed Binned Effect (logit) \n of Today's Temperature (Lag=0d)",
                                  xlab = "Temperature °C", ylab = "Days post-Discharge"),
               plot.axes = { axis(1); axis(2); box(); overlay_points(x = tempx -273.15, y = seq(2, 11), 
                                                                     hrmat[,3:12],lrmat[,c(-1,-2)], 8,1)})

#dev.new()
#filled.contour(x = tempx -273.15, y = seq(2, lag), 
#               z = hrmat[,c(-1,-2)], col = col, levels = levels, 
#               plot.title = title(main = "HIGH Effect (logit) \n of Today's Temperature (Lag=0d)",
#                                  xlab = "Temperature °C", ylab = "Days post-Discharge"))


#diag(x) - today's temperature
# y axis: current days since discharge



##PLOTTING: summaries

# string of constant temperatures - x axis
# y axis: current days since discharge
# height = added risk at that time

# string of constant temperatures - x axis
# height = total risk summing along all days


#if (ptype == "contour") {
#  if (x$lag[2] == 0) 
#    stop("contour plot not conceivable for unlagged associations")


#if (ptype == "3d") {
#  if (diff(x$lag) == 0) 
#    stop("3D plot not conceivable for unlagged associations")
#  plot.arg <- list(ticktype = "detailed", theta = 210, 
#                   phi = 30, xlab = "Var", ylab = "Lag", zlab = "Outcome", 
#                   col = "lightskyblue", zlim = c(min(x$matfit), max(x$matfit)), 
#                   ltheta = 290, shade = 0.75, r = sqrt(3), d = 5)
#  plot.arg <- modifyList(plot.arg, list(...))
#  plot.arg <- modifyList(plot.arg, list(x = x$predvar, 
#                                        y = seqlag(x$lag, x$bylag), z = x$matfit))
#  do.call("persp", plot.arg)
#}
