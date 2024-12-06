
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

################## Construct Time Dependent Covariate Dataset ################################

# Set up the initial dataset with the main survival information
burn2 <- tmerge(data1 = burn1, data2 = burn1, id = Obs, 
                tstop = T3, infection = event(T3, D3))

# Add surgical excision (T1) as a time-dependent covariate
burn2 <- tmerge(burn2, burn1, id = Obs, surgical = tdc(T1))

# Add antibiotic treatment (T2) as a time-dependent covariate
burn2 <- tmerge(burn2, burn1, id = Obs, antibiotics = tdc(T2))

# Define status to distinguish between the event and censoring
burn2$status <- as.integer(with(burn2, (tstop == T3 & D3 == 1)))

# Ensure the dataset is in the counting process format
burn2.surv <- with(burn2, Surv(time = tstart, time2 = tstop, event = status, type = "counting"))

# View the resulting dataset
head(burn2)

################## Fit Cox Models ################################

# Fit 4 cox models  with or w/o each time-dependent covariates
fit_time_indep <- coxph(burn2.surv ~ Treatment + Gender + Race, data = burn2)
fit_time_dep <- coxph(burn2.surv ~ Treatment + Gender + Race + surgical + antibiotics, data = burn2)
fit_time_dep_sur <- coxph(burn2.surv ~ Treatment + Gender + Race + surgical, data = burn2)
fit_time_dep_anti <- coxph(burn2.surv ~ Treatment + Gender + Race  + antibiotics, data = burn2)

AIC(fit_time_indep, fit_time_dep, fit_time_dep_sur, fit_time_dep_anti)
drop1(fit_time_indep, test = "Chisq")
drop1(fit_time_dep, test = "Chisq")
drop1(fit_time_dep_sur, test = "Chisq")
drop1(fit_time_dep_anti, test = "Chisq")

# only retain surgical as time dependent variable

# final model: fit_time_dep_sur
fit_time_dep_sur <- coxph(burn2.surv ~ Treatment + Race + surgical, data = burn2)


################## Model Checking ###############################


# Schoenfeld residuals
cox_zph <- cox.zph(fit_time_dep_sur)

par(mfrow = c(2, 2)) 
plot(cox_zph[1], main = "Schoenfeld Residuals for Treatment")
plot(cox_zph[2], main = "Schoenfeld Residuals for Gender")
plot(cox_zph[3], main = "Schoenfeld Residuals for Race")
plot(cox_zph[4], main = "Schoenfeld Residuals for surgical")



par(mfrow = c(1, 1))

# Martingale residuals
martingale_residuals <- residuals(fit_time_dep_sur, type = "martingale")
linear_predictor <- predict(fit_time_dep_sur, type = "lp")

plot(linear_predictor, martingale_residuals,
     xlab = "Linear Predictor", ylab = "Martingale Residuals",
     main = "Martingale Residuals vs. Linear Predictor")
abline(h = 0, col = "red", lty = 2)
text(linear_predictor, martingale_residuals,
     labels = 1:length(linear_predictor), pos = 4, cex = 0.6)

# Deviance residuals
deviance_residuals <- residuals(fit_time_dep_sur, type = "deviance")

plot(linear_predictor, deviance_residuals,
     xlab = "Linear Predictor", ylab = "Deviance Residuals",
     main = "Deviance Residuals vs. Linear Predictor")
abline(h = 0, col = "red", lty = 2)
text(linear_predictor, deviance_residuals,
     labels = 1:length(linear_predictor), pos = 4, cex = 0.6)




# Cox-Snell Residuals: Assess Goodness of Fit
cox_snell_residuals <- -log(survfit(fit_time_dep_sur)$surv)  # Compute Cox-Snell residuals
cox_snell_residuals <- cox_snell_residuals[1:length(burn2$status)]  # Match length to data

# Create a survival object for Cox-Snell residuals
burn2$status <- as.integer(burn2$status)
surv_coxsnell <- Surv(time = cox_snell_residuals, event = burn2$status)

# Fit a survival model for Cox-Snell residuals
fit_coxsnell <- survfit(surv_coxsnell ~ 1)

# Plot cumulative hazard for Cox-Snell residuals
plot(fit_coxsnell, fun = "cumhaz",
     xlab = "Cox-Snell Residuals", ylab = "Cumulative Hazard",
     main = "Cumulative Hazard of Cox-Snell Residuals")
abline(0, 1, col = "red", lty = 3)  # Reference line for goodness-of-fit


# dfbeta
par(mfrow = c(2, 2))
dfbeta <- residuals(fit_time_dep_sur, type = "dfbeta")

covariate_names <- c("Treatment", "Gender", "Race", "surgical", "antibiotics")


for (i in 1:5) {
  # Plot DFBETA values for the i-th covariate
  plot(dfbeta[, i], 
       xlab = "Observation", ylab = "DFBETA", 
       main = paste("DFBETA for", covariate_names[i]))
  
  # Add observation number labels to points
  text(x = 1:nrow(dfbeta), y = dfbeta[, i], 
       labels = 1:nrow(dfbeta), pos = 4, cex = 0.6)
}



########### Questions ###########

# shall i use time dependent covariates? no improvement in AIC

### examine TD variables individually

# > AIC(fit_independent, fit_refined, fit_time_dep)
# df      AIC
# fit_independent 10 438.5335
# fit_refined      3 367.8405
# fit_time_dep     5 427.9481

# how to adjust for colinearity? Does drop1 handle colinearity?
# construct model with each TD variable and do drop1.
# 4 models 
# 1. no td var
# 2. td var 1
# 3. td var 2
# 4. tdvar 1, 2

# other ways to assess goodness of fit other than cox-snell?
# check assumptions is the most important
# check general assumption if violated


# language in final report? is current ok?

# outliers -- patient covariant -- report to burn surgon -- 

# seperate 



# Martingale residuals
martingale_residuals <- residuals(fit_time_dep_sur, type = "martingale")

# Deviance residuals
deviance_residuals <- residuals(fit_time_dep_sur, type = "deviance")

# Linear predictor (log hazard)
linear_predictor <- predict(fit_time_dep_sur, type = "lp")

# Identify extreme deviance residuals (e.g., threshold > 2 standard deviations)
threshold_dev <- 2 * sd(deviance_residuals)
extreme_deviance <- which(abs(deviance_residuals) > threshold_dev)

# Identify extreme martingale residuals (e.g., threshold > 2 standard deviations)
threshold_mart <- 2 * sd(martingale_residuals)
extreme_martingale <- which(abs(martingale_residuals) > threshold_mart)

# DFBETA values for each predictor
dfbeta_values <- residuals(fit_time_dep_sur, type = "dfbeta")

# Identify extreme DFBETA values (threshold > 2 standard deviations for any predictor)
threshold_dfbeta <- apply(dfbeta_values, 2, function(x) 2 * sd(x))
extreme_dfbeta <- which(apply(abs(dfbeta_values), 1, function(x) any(x > threshold_dfbeta)))


# Define a function to extract rows for any residual type
extract_extreme <- function(dataset, indices) {
  indices <- unique(indices)  # Ensure no duplicates
  return(dataset[indices, ])
}

# Apply the function for martingale, deviance, and DFBETA residuals
extreme_martingale_data <- extract_extreme(burn2, extreme_martingale)
extreme_deviance_data <- extract_extreme(burn2, extreme_deviance)
extreme_dfbeta_data <- extract_extreme(burn2, extreme_dfbeta)




# Martingale Residuals
extreme_martingale_values <- martingale_residuals[extreme_martingale]  # Extract values
extreme_martingale_data <- burn2[extreme_martingale, ]                 # Extract rows
extreme_martingale_data$Martingale <- extreme_martingale_values        # Add as new column

# Deviance Residuals
extreme_deviance_values <- deviance_residuals[extreme_deviance]        # Extract values
extreme_deviance_data <- burn2[extreme_deviance, ]                    # Extract rows
extreme_deviance_data$Deviance <- extreme_deviance_values             # Add as new column

# DFBETA Values
extreme_dfbeta_values <- dfbeta_values[extreme_dfbeta, ]              # Extract values (matrix)
extreme_dfbeta_data <- burn2[extreme_dfbeta, ]                        # Extract rows
extreme_dfbeta_data <- cbind(extreme_dfbeta_data, extreme_dfbeta_values)  # Add all DFBETA columns

