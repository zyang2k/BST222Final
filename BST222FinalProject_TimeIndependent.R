
################## Data Cleaning ################################

# load necessary libraries
library(KMsurv)
library(survival)
library(MASS)

# load burn data
data(burn)

# burn dataframe with new names
burn1 <- burn
burn1 <- data.frame(burn1,Treatment=factor(burn1$Z1,labels=c("Routine","Cleansing")))
burn1 <- data.frame(burn1,Gender=factor(burn1$Z2,labels=c("Male","Female")))
burn1 <- data.frame(burn1,Race=factor(burn1$Z3,labels=c("Nonwhite","White")))
burn1 <- data.frame(burn1,PercentBurned=burn1$Z4)
burn1 <- data.frame(burn1,SiteHead=factor(burn1$Z5,labels=c("NotBurned","Burned")))
burn1 <- data.frame(burn1,SiteButtock=factor(burn1$Z6,labels=c("NotBurned","Burned")))
burn1 <- data.frame(burn1,SiteTrunk=factor(burn1$Z7,labels=c("NotBurned","Burned")))
burn1 <- data.frame(burn1,SiteUpperLeg=factor(burn1$Z8,labels=c("NotBurned","Burned")))
burn1 <- data.frame(burn1,SiteLowerLeg=factor(burn1$Z9,labels=c("NotBurned","Burned")))
burn1 <- data.frame(burn1,SiteRespTract=factor(burn1$Z10,labels=c("NotBurned","Burned")))
burn1 <- data.frame(burn1,BurnType=factor(burn1$Z11,labels=c("Chemical","Scald","Electric","Flame")))

# remove unuseful columns
burn1 <- burn1[, -c(2:12)]

################## KM and NA surv ###############################

# Create the survival object
burn1.surv <- Surv(burn1$T3, burn1$D3)

# Fit KM and NA curves
KMcurves <- survfit(burn1.surv ~ Treatment, data = burn1)
NAcurves <- survfit(burn1.surv ~ Treatment, type = "fleming-harrington", data = burn1)

# Plot the KM curves
plot(KMcurves, col = c("black", "red"), lwd = 2, lty = 1, xlab = "Time", ylab = "Survival Probability", 
     main = "KM and NA Survival Curves for Two Types of Treatment")

# Add the NA curves
lines(NAcurves, col = c("black", "red"), lwd = 2, lty = 2)

# Add a legend
legend("bottomright", 
       legend = c("Routine Care (KM)", "Cleansing (KM)", "Routine Care (NA)", "Cleansing (NA)"),
       col = c("black", "red", "black", "red"), 
       lty = c(1, 1, 2, 2), 
       lwd = 2)

# Plot NA cumulative hazard vs time
plot(NAcurves, col = c("black", "red"), lwd = 2, fun = "cumhaz")
title("NA Cumulative Hazard vs. time for Two Treatment Groups")
legend("bottomright", c("Routine Care", "Cleansing"), col = c("black", "red"), lwd = 2)


# Plot complimentary log-log survival vs. log time
plot(NAcurves, fun = "cloglog", col = 1:2, lwd = 2, xlab = "time (days in log scale)", ylab = "log cumulative hazard")
legend("bottomright", legend = levels(burn1$Treatment), col = 1:2, lwd = 2)
title("Complementary Log-Log Survival Curves")

################## stepwise cox ###############################

# Construct Cox model on time-independent predictors
# null model
fit_0 <- coxph(burn1.surv ~ Treatment, data = burn1)

# Perform forward stepwise regression
step_model <- stepAIC(
  fit_0,
  scope = list(upper = ~Obs + Treatment + Gender + Race + PercentBurned + SiteHead + 
                 SiteButtock + SiteTrunk + SiteUpperLeg + SiteLowerLeg + SiteRespTract + 
                 BurnType, 
               lower = ~ Treatment),
  direction = "both",
  k = 2
)

# View the summary of the final model
summary(step_model)


# Based on the stepAIC result, we exclude burn site as a group




################## model checking ###############################

# step model checking

# test proportionality via cox.zph
# seems to be significant non-proportionality

burn1.zph <- cox.zph(step_model)
print(burn1.zph)

par(mfrow = c(2, 2))

plot(burn1.zph[1], 
     main = "Schoenfeld Residuals for Treatment Type")
plot(burn1.zph[2], 
     main = "Schoenfeld Residuals for Gender Type")
plot(burn1.zph[3], 
     main = "Schoenfeld Residuals for Race Type")
plot(burn1.zph[4], 
     main = "Schoenfeld Residuals for BurnType Type")



# cox-snell residual
par(mfrow = c(1, 1))

# calculate the various residuals 
burn1.mart <- residuals(step_model, type = "martingale")
burn1.dev <- residuals(step_model,type="deviance")
burn1.dfb <- residuals(step_model,type="dfbeta")
burn1.cs <- burn1$D3 - burn1.mart


# find linear predictor
burn1.preds <- predict(step_model)



# cumulative hazard of cox-snell residuals
surv.csr <- survfit(Surv(burn1.cs,burn1$D3)~1,
                    type="fleming-harrington")

plot(surv.csr, fun = "cumhaz")
abline(a = 0, 
       b = 1, 
       col = "red")
title("Cumulative Hazard of Cox-Snell Residuals")

# martingale residual vs. predictor
plot(burn1.preds, burn1.mart, 
     xlab = "linear predictor", 
     ylab = "Martingale Residual")
text(burn1.preds, burn1.mart, labels = rownames(burn1))
title("Martingale Residuals vs. Linear Predictor")


# deviance vs. predictor
plot(burn1.preds, burn1.dev, 
     xlab = "linear predictor", ylab = "Deviance Residual")
text(burn1.preds, burn1.dev, labels = rownames(burn1))
title("Deviance Residuals vs. Linear Predictor")

# dfbeta vs. observation order
par(mfrow = c(2, 2))

plot1A <- function() {
  plot(burn1.dfb[,1], xlab = "Observation Number",
       ylab = "dfbeta for Treatment Type")
  text(burn1.dfb[,1], labels = rownames(burn1))
  title("dfbeta Values for Treatment Type")
}

plot1A()

plot2A <- function() {
  plot(burn1.dfb[,2], xlab = "Observation Number",
       ylab = "dfbeta for Gender Type")
  text(burn1.dfb[,2], labels = rownames(burn1))
  title("dfbeta Values for Gender Type")
}

plot2A()

plot3A <- function() {
  plot(burn1.dfb[,3], xlab = "Observation Number",
       ylab = "dfbeta for Race Type")
  text(burn1.dfb[,3], labels = rownames(burn1))
  title("dfbeta Values for Race Type")
}

plot3A()

plot4A <- function() {
  plot(burn1.dfb[,4], xlab = "Observation Number",
       ylab = "dfbeta for Burn Type")
  text(burn1.dfb[,4], labels = rownames(burn1))
  title("dfbeta Values for Burn Type")
}

plot4A()



