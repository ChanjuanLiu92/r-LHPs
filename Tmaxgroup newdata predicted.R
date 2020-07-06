
rm(list=ls())

library(lattice)
library("R2WinBUGS")	# Load the R2WinBUGS library
library("lme4")

### 1. Data generation
# Generate two samples of body mass measurements of male peregrines
dat=read.csv("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\data\\162data\\rdata162.csv")
r=dat$r
Tmax=dat$Tmax
group=dat$group
n=162

library(xlsx)
data_test=read.xlsx("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\data\\162data\\new 10 stocks.xlsx",1,title=T)
m=nrow(data_test)
r.test.ob=data_test$r
Tmax.test=data_test$Tmax
group.test=data_test$group
group.test
### 2.Analysis using WinBUGS
#R work directory
setwd("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxgroup") # May have to adapt that

# Save BUGS description of the model to working directory
sink("model.txt")
cat("
model {
# hyper-prior
mu.bTmax ~ dnorm(0, 0.01)
tau.bTmax~ dgamma (0.5,0.01)

#prior
for(gp in 1:3){
bTmax[gp]~dnorm(mu.bTmax, tau.bTmax)
}


tau.Tmax ~ dgamma(0.01,0.01)
tau.r~dgamma(0.01,0.01)


# Likelihood	
for (i in 1:n) 	{
LgTmax[i]<- log(Tmax[i])
estTmax[i] ~ dlnorm(LgTmax[i], tau.Tmax) 

r[i] ~ dnorm(est.r[i], tau.r)
est.r[i] <-bTmax[group[i]]/estTmax[i]

predicted[i]<-est.r[i]
residual[i] <- r[i] - est.r[i]
sq[i] <- pow(residual[i], 2)                         # Squared residuals for observed data

# Generate replicate data and compute fit stats for them
r.new[i] ~ dnorm(est.r[i], tau.r)                      # one new squared residuals for actual data set
sq.new[i] <- pow(r.new[i] - predicted[i], 2)       # Squared residuals for new data

}
fit <- sum(sq[ ])                                              # Sum of squared residuals for actual data
fit.new <- sum(sq.new[ ])                                # Sum of squared residuals for new data set
test <- step(fit.new - fit)                                 # Test whether new data set more exteme
bpvalue <- mean(test) 
# Bayesian p - value

#predicted new data
for (j in 1:m) 	{
LgTmax.test[j]<- log(Tmax.test[j])
estTmax.test[j] ~ dlnorm(LgTmax.test[j], tau.Tmax)
rnew.test[j]~ dnorm(est.rnew.test[j], tau.r)
est.rnew.test[j]<-bTmax[group.test[j]]/estTmax.test[j]
}


}
",fill=TRUE)
sink()


# Package all the stuff to be handed over to WinBUGS
# Bundle data
win.data <- list("r","Tmax","group","n","Tmax.test","m","group.test")

# Function to generate starting values
inits <- function()
  list (
    mu.bTmax=1,tau.bTmax=1,
    bTmax=c(1,1,1),
    tau.Tmax=1,tau.r=1,
    estTmax=dat$Tmax)


# Parameters to be monitored (= to estimate)
# params <- c("bTmax","mu.bTmax","tau.bTmax","tau.Tmax","tau.r","bpvalue","predicted",
#             "estTmax","residual")
params <- c("bTmax","rnew.test","estTmax.test","bpvalue","mu.bTmax","tau.bTmax","tau.Tmax","tau.r",
            "estTmax","residual")

# MCMC settings
nc <- 2				# Number of chains
ni <- 50000				# Number of draws from posterior (for each chain)
nb <- 40000					# Number of draws to discard as burn-in
nt <- 1					# Thinning rate

# Start Gibbs sampler: Run model in WinBUGS and save results in object called out
out <- bugs(data = win.data, inits = inits, parameters.to.save = params, model.file = "model.txt", 
            n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, debug = F, DIC = TRUE, 
            bugs.directory ="F:\\winbugs143_unrestricted\\winbugs14_full_patched\\WinBUGS14\\")


#！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
#results
#！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
ls()
out					# Produces a summary of the object
names(out)

stat=data.frame(out$pD,out$DIC)
coe=out$summary

# errors on test data
Bias=data.frame(r.test.ob-out$mean$rnew.test)
MAE<- mean(abs(r.test.ob-out$mean$rnew.test))
MSE<-mean(((r.test.ob-out$mean$rnew.test)^2))
RMSE<-sqrt(mean(((r.test.ob-out$mean$rnew.test)^2)))
RE<-data.frame((r.test.ob-out$mean$rnew.test)/(r.test.ob))
MARE=mean(abs(RE[,1]))
BHEIV_errors=data.frame(stat,MAE,MSE,RMSE,MARE)
bias_re=cbind(Bias,RE)


# Fishlife estimated results
FL_r=data_test$r.FishLife
Bias=data.frame(r.test.ob-FL_r)
MAE<- mean(abs(r.test.ob-FL_r))
MSE<-mean(((r.test.ob-FL_r)^2))
RMSE<-sqrt(mean(((r.test.ob-FL_r)^2)))
RE<-data.frame((r.test.ob-FL_r)/(r.test.ob))
MARE=mean(abs(RE[,1]))
FL_errors=data.frame(MAE,MSE,RMSE,MARE)

# 
write.xlsx(stat,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\data\\162data\\M6 new 10 stocks.xlsx",sheetName="stat",append=T)
write.xlsx(coe,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\data\\162data\\M6 new 10 stocks.xlsx",sheetName="coe",append=T)
write.xlsx(BHEIV_errors,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\data\\162data\\M6 new 10 stocks.xlsx",sheetName="BHEIV_errors",append=T)
write.xlsx(FL_errors,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\data\\162data\\M6 new 10 stocks.xlsx",sheetName="FL_errors",append=T)


