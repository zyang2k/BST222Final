
################## Data Cleaning ################################

# load necessary libraries
library(KMsurv)
library(survival)
library(MASS)
library(car)

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
par(mfrow = c(1, 1))
burn1.surv <- Surv(burn1$T3, burn1$D3)
surv_diff <- survdiff(burn1.surv ~ Treatment, data = burn1)
surv_diff
# Fit KM and NA curves
KMcurves <- survfit(burn1.surv ~ Treatment, data = burn1)
NAcurves <- survfit(burn1.surv ~ Treatment, type = "fleming-harrington", data = burn1)

# Plot the KM curves
plot(KMcurves, col = c("black", "red"), lwd = 2, lty = 1, xlab = "Time (days)", ylab = "Percent Infection-free", 
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

################## Cox Model ###############################
# Construct Cox model on all time-independent predictors
fit_independent <- coxph(burn1.surv ~ Treatment + Gender + Race + PercentBurned +
                           SiteHead + SiteButtock + SiteTrunk +
                           SiteUpperLeg + SiteLowerLeg + SiteRespTract, 
                         data = burn1)

################## Cox Model ###############################
# Summarize the model
summary(fit_independent)


# Check proportional hazards assumption
cox.zph(fit_independent)

# assumption holds


# Plot Schoenfeld residuals to assess proportional hazards assumption
plot(cox.zph(fit_independent))
# schoenfeld residual not look too bad

# compare 2 models
AIC(fit_independent)

# stratified better

# use drop1 to test which to drop
drop1(fit_independent, test = "Chisq")

# drop 1 result indicate a only treatment and race matters, multicliearity? 

# use independent, drop unsignificant variables
fit_refined <- coxph(burn1.surv ~ Treatment + Gender + Race,
                     data = burn1)
summary(fit_refined)

# recheck proprtional assumption -- passed

################## Model Checking ###############################
# Schoenfeld residuals to test for proportionality via cox.zph
cox_zph <- cox.zph(fit_refined)

# Plot Schoenfeld residuals for each covariate
par(mfrow = c(2, 2))  # Adjust plot layout
plot(cox_zph[1], main = "Schoenfeld Residuals for Treatment")
plot(cox_zph[2], main = "Schoenfeld Residuals for Gender")
plot(cox_zph[3], main = "Schoenfeld Residuals for Race")
# Reset plotting layout
par(mfrow = c(1, 1))



# Predicted cumulative hazards for each individual
cum_haz <- predict(fit_refined, type = "expected")
# Cox-Snell residuals
cox_snell_residuals <- cum_haz
# Create a survival object for Cox-Snell residuals
surv_coxsnell <- Surv(cox_snell_residuals, burn1$D3)

# Fit the survival model
fit_coxsnell <- survfit(surv_coxsnell ~ 1)

# Plot cumulative hazard
plot(fit_coxsnell, fun = "cumhaz",
     xlab = "Cox-Snell Residuals", ylab = "Cumulative Hazard",
     main = "Cumulative Hazard of Cox-Snell Residuals")
abline(0, 1, col = "red", lty = 3)  # Reference line for goodness-of-fit




# Martingale residuals to assess nonlinearity
martingale_residuals <- residuals(fit_refined, type = "martingale")
linear_predictor <- predict(fit_refined, type = "lp")

plot(linear_predictor, martingale_residuals,
     xlab = "Linear Predictor", ylab = "Martingale Residuals",
     main = "Martingale Residuals vs. Linear Predictor")
abline(h = 0, col = "red", lty = 2)
text(linear_predictor, martingale_residuals,
     labels = 1:length(linear_predictor),  # Observation numbers as labels
     pos = 4, cex = 0.7, col = "black")    # Adjust position, size, and color



# Deviance residuals to detect outliers
deviance_residuals <- residuals(fit_refined, type = "deviance")

plot(linear_predictor, deviance_residuals,
     xlab = "Linear Predictor", ylab = "Deviance Residuals",
     main = "Deviance Residuals vs. Linear Predictor")
abline(h = 0, col = "red", lty = 2)
text(linear_predictor, deviance_residuals,
     labels = 1:length(linear_predictor),  # Observation numbers as labels
     pos = 4, cex = 0.7, col = "black")    # Adjust position, size, and color





# DFBETA values to check for influential observations
dfbeta_values <- residuals(fit_refined, type = "dfbeta")
# Extract covariate names from the model
covariate_names <- names(coef(fit_refined))
# Assign column names to dfbeta_values
colnames(dfbeta_values) <- covariate_names
# Plot DFBETA values by observation order for each covariate
par(mfrow = c(2, 2))  # Set up a 2x2 plotting grid
for (i in 1:ncol(dfbeta_values)) {
  # Plot DFBETA values
  plot(dfbeta_values[, i], 
       xlab = "Observation", ylab = "DFBETA",
       main = paste("DFBETA for", colnames(dfbeta_values)[i]))  # Customize points
  
  # Add reference line at 0
  abline(h = 0, col = "red", lty = 2)
  
  # Label each point with its observation number
  text(1:nrow(dfbeta_values), dfbeta_values[, i], 
       labels = 1:nrow(dfbeta_values), pos = 4, cex = 0.6, col = "black")
}

# identify interesting observations


