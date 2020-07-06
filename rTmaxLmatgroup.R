rm(list=ls())

library(lattice)
library("R2WinBUGS")	# Load the R2WinBUGS library
library("lme4")

### 1. Data generation
# Generate two samples of body mass measurements of male peregrines
dat=read.csv("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\data\\162data\\rdata162.csv")
r=dat$r
Tmax=dat$Tmax
Lmat=dat$Lmat
group=dat$group
n=162


### 2.Analysis using WinBUGS
#R work directory
setwd("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLmatgroup") # May have to adapt that

# Save BUGS description of the model to working directory
sink("model.txt")
cat("
model {
# hyper-prior
mu.bTmax ~ dnorm(0, 0.0001)
tau.bTmax~ dgamma (0.5,0.01)

mu.bLmat ~ dnorm(0, 0.0001)
tau.bLmat~ dgamma (0.5,0.1)

#prior
for(gp in 1:3){
bTmax[gp]~dnorm(mu.bTmax, tau.bTmax)
bLmat[gp]~dnorm(mu.bLmat, tau.bLmat)
}

tau.Tmax ~ dgamma(0.01,0.1)
tau.Lmat ~ dgamma(50,0.1)
tau.r~dgamma(0.01,0.01)

# Likelihood	
for (i in 1:n) 	{

LgTmax[i]<- log(Tmax[i])
estTmax[i] ~ dlnorm(LgTmax[i], tau.Tmax) 

LgLmat[i]<- log(Lmat[i])
estLmat[i] ~ dlnorm(LgLmat[i], tau.Lmat) 


r[i] ~ dnorm(est.r[i], tau.r)

est.r[i] <- bTmax[group[i]]/estTmax[i]+pow(estLmat[i],bLmat[group[i]])



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
bpvalue <- mean(test)                                    # Bayesian p - value
}
",fill=TRUE)
sink()


# Package all the stuff to be handed over to WinBUGS
# Bundle data
win.data <- list("r","Tmax","Lmat","group","n")

# Function to generate starting values
inits <- function()
  list (
    mu.bTmax=1,tau.bTmax=1,
    mu.bLmat=0,tau.bLmat=1,
    bTmax=c(1,1,1),bLmat=c(1,1,1),
    tau.Tmax=1,estTmax=dat$Tmax,
    tau.Lmat=1,estLmat=dat$Lmat,
    tau.r=1)


# Parameters to be monitored (= to estimate)
params <- c("bTmax","bLmat",
            "mu.bTmax","tau.bTmax","tau.Tmax",
            "mu.bLmat","tau.bLmat","tau.Lmat",
            "tau.r","bpvalue","predicted", 
            "estTmax","estLmat","residual")
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
#winbugs should be closed after running the bus function
ls()
out					# Produces a summary of the object
names(out)

stat=data.frame(out$pD,out$DIC,out$mean$bpvalue)
coe=out$summary
coe_stat=as.matrix(summary(out$sims.list$bTmax,
                           out$sims.list$bLmat))

dic = out$DIC
pd = out$pD

require(devEMF)
emf("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLmatgroup\\residual plot mix without beta0.emf",
    width = 6.7, height = 3.35)
par(mfrow=c(1,2))
plot(out$mean$predict~ dat$r,xlab="Observed r",ylab="Predicted r",cex=0.8); abline(0,1,col="red")
plot(out$mean$residual~ dat$r, xlab="Observed r",ylab="Residuals",cex=0.8,pch='x',ylim = c(-0.5,0.5));abline(0,0,col="red")
dev.off()


Bias<-mean(dat$r-out$mean$predicted)
MAE<- mean(abs(out$mean$predicted-dat$r))
MSE<-mean(((out$mean$predicted-dat$r)^2))
RMSE<-sqrt(mean(((out$mean$predicted-dat$r)^2)))
PE<-(dat$r-out$mean$predicted)/(dat$r)
MPE=mean(PE)
MAPE=mean(abs(PE))
R=cor(out$mean$predicted,dat$r)
result=data.frame(Bias,MAE,MSE,RMSE,PE,MPE,MAPE,R)



library(loo)
LLmat=out$sims.matrix
dim(LLmat)
LLmat=out$sims.matrix[1:20000,1:170]
dim(LLmat)
rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1:2, each = 10000))
loo_results=loo(LLmat, r_eff = rel_n_eff)
loo_mat=t(loo_results$estimates)
loo_mat

WAIC_mat= waic(LLmat)
waic_mat=t(WAIC_mat$estimates)
waic_mat

library(xlsx)
write.xlsx(stat,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLmatgroup\\mix without beta0.xlsx",sheetName = "stat",row.names=T)
write.xlsx(coe,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLmatgroup\\mix without beta0.xlsx",sheetName = "coe",row.names=T,append=T)
write.xlsx(result,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLmatgroup\\mix without beta0.xlsx",sheetName = "result",row.names=T,append=T)
write.xlsx(loo_mat,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLmatgroup\\mix without beta0.xlsx",sheetName = "loo_estimates",row.names=T,append=T)



#！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
# require(devEMF)
emf("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLmatgroup\\mix without beta0.emf",
    width = 3.35, height = 3.35)
observed=r
BHEIV_predicted=out$mean$predicted#BHEIV model predicted r
plot(observed,col="red",pch=1, main = "rTmaxLinfgroup mix without beta0",cex.main=0.8,
     xlab = "index",ylab="r")
points(BHEIV_predicted,col="blue",pch=2)
legend("topright",legend = c("observed r","BHEIV_predicted r"),col=c("red","blue"),pch=c(1,2),cex=0.7)
dev.off()


sigma_r=sqrt(1/out$mean$tau.r)
r_com=data.frame(observed,BHEIV_predicted,sigma_r)
write.xlsx(r_com,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLinfgroup\\mix without beta0.xlsx",sheetName = "r_com",row.names=T,append=T)



est_Tmax=out$mean$estTmax
sigma_Tmax=sqrt(1/out$mean$tau.Tmax)

est_Lmat=out$mean$estLmat
sigma_Lmat=sqrt(1/out$mean$tau.Lmat)

p_com=data.frame(Tmax,est_Tmax,sigma_Tmax,
                 Lmat,est_Lmat,sigma_Lmat)
write.xlsx(p_com,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLinfgroup\\mix without beta0.xlsx",sheetName = "p_com",row.names=T,append=T)







# # # BGR convergence diagnostic called Rhat, and an
# #BGR喘栖登僅頁倦辺漸��自峙頁1.1匆辛恬葎辺漸
# hist(out$summary[,8])			# Rhat values in the eighth column of the summary
# which(out$summary[,8] > 1.1)		# None in this case
# # # # For trace plots for the entire chains, do
# matplot(out$sims.array[1000:5000,1:3,1], type = "l")
# # # # We can also produce graphical summaries, e.g., histograms of the posterior
# # # # distributions for each parameter:
# data=data.frame(out$sims.list$bM);names(data)=c("b1M","b2M","b3M")
# library(ggplot2)
# cols=c('Invertebrate'="#f04546", 'Elasmobranch'="#33CC33", 'Teleost'="#3591d1")
# shapes=c('Invertebrate'=1,'Elasmobranch'=3,'Teleost'=5)
# density_plot=
#   ggplot(data)+
#   geom_density(aes(b1M,color="Invertebrate",linetype="Invertebrate"))+
#   geom_density(aes(b2M,color="Elasmobranch",linetype="Elasmobranch"))+
#   geom_density(aes(b3M,color="Teleost",linetype="Teleost"))+
#   labs(y="Density",x=expression(italic(beta[M])))+
#   theme_bw()+
#   scale_colour_manual(name="group",values=cols)+
#   scale_linetype_manual(name="group",values = shapes)+
#   theme(legend.position=c(.99, .99),legend.justification = c("right", "top"))+
#   theme(legend.key = element_blank())+theme(legend.background = element_blank())
# plot(density_plot)
# 
# library(devEMF)
# emf("G:\\02-r_LHPs paper\\results20200222\\rMgroup\\density_plot.emf",width = 3.35, height = 2.54)
# plot(density_plot)
# dev.off()
# 
# 
# # pdf("G:\\02-r_LHPs paper\\results20200222\\rMgroup\\density_plot.pdf",width = 3.35, height = 2.54)
# # plot(density_plot)
# # dev.off()
# 
# # tiff("G:\\02-r_LHPs paper\\results20200222\\rMgroup\\density_plot.tif",width = 335, height = 254,compression = c("none"))
# # png("G:\\02-r_LHPs paper\\results20200222\\rMgroup\\density_plot.png",width = 335, height = 254)
# 
# 
