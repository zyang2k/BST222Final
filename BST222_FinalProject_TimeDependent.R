
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

# Fit the Cox model with time-dependent covariates
fit_time_dep <- coxph(burn2.surv ~ Treatment + Gender + Race + surgical + antibiotics, data = burn2)

# Summarize the results
summary(fit_time_dep)

drop1(fit_time_dep, test = "Chisq")
vif(fit_time_dep)


################## Model Checking ###############################


# Schoenfeld residuals
cox_zph <- cox.zph(fit_time_dep)

par(mfrow = c(2, 3)) 
plot(cox_zph[1], main = "Schoenfeld Residuals for Treatment")
plot(cox_zph[2], main = "Schoenfeld Residuals for Gender")
plot(cox_zph[3], main = "Schoenfeld Residuals for Race")
plot(cox_zph[4], main = "Schoenfeld Residuals for surgical")
plot(cox_zph[5], main = "Schoenfeld Residuals for antibiotics")



par(mfrow = c(1, 1))

# Martingale residuals
martingale_residuals <- residuals(fit_time_dep, type = "martingale")
linear_predictor <- predict(fit_time_dep, type = "lp")

plot(linear_predictor, martingale_residuals,
     xlab = "Linear Predictor", ylab = "Martingale Residuals",
     main = "Martingale Residuals vs. Linear Predictor")
abline(h = 0, col = "red", lty = 2)
text(linear_predictor, martingale_residuals,
     labels = 1:length(linear_predictor), pos = 4, cex = 0.6)

# Deviance residuals
deviance_residuals <- residuals(fit_time_dep, type = "deviance")

plot(linear_predictor, deviance_residuals,
     xlab = "Linear Predictor", ylab = "Deviance Residuals",
     main = "Deviance Residuals vs. Linear Predictor")
abline(h = 0, col = "red", lty = 2)
text(linear_predictor, deviance_residuals,
     labels = 1:length(linear_predictor), pos = 4, cex = 0.6)




# Cox-Snell Residuals: Assess Goodness of Fit
cox_snell_residuals <- -log(survfit(fit_time_dep)$surv)  # Compute Cox-Snell residuals
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
par(mfrow = c(2, 3))
dfbeta <- residuals(fit_time_dep, type = "dfbeta")

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
# > AIC(fit_independent, fit_stratified, fit_refined, fit_time_dep)
# df      AIC
# fit_independent 10 438.5335
# fit_stratified   9 376.4142
# fit_refined      3 367.8405
# fit_time_dep     5 427.9481

# how to adjust for colinearity? Does drop1 handle colinearity?

# other ways to assess goodness of fit other than cox-snell?

# language in final report? is current ok?

