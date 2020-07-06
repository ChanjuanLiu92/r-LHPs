rm(list=ls())

library(lattice)
library("R2WinBUGS")	# Load the R2WinBUGS library
library("lme4")

### 1. Data generation
# Generate two samples of body mass measurements of male peregrines
dat=read.csv("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\data\\162data\\rdata162.csv")
r=dat$r
Tmax=dat$Tmax
Linf=dat$Linf
group=dat$group
n=162


### 2.Analysis using WinBUGS
#R work directory
setwd("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLinfgroup") # May have to adapt that

# Save BUGS description of the model to working directory
sink("model.txt")
cat("
model {
# hyper-prior
mu.bTmax ~ dnorm(0, 0.0001)
tau.bTmax~ dgamma (0.5,0.01)

mu.bLinf ~ dnorm(0, 0.0001)
tau.bLinf~ dgamma (0.5,0.1)

#prior
for(gp in 1:3){
bTmax[gp]~dnorm(mu.bTmax, tau.bTmax)
bLinf[gp]~dnorm(mu.bLinf, tau.bLinf)
}

tau.Tmax ~ dgamma(0.01,0.1)
tau.Linf ~ dgamma(50,0.1)
tau.r~dgamma(0.01,0.01)

sigma_r<-sqrt(1/tau.r)
sigma_Tmax<-sqrt(tau.Tmax)
sigma_Linf<-sqrt(tau.Linf)


# Likelihood	
for (i in 1:n) 	{

LgTmax[i]<- log(Tmax[i])
estTmax[i] ~ dlnorm(LgTmax[i], tau.Tmax) 

LgLinf[i]<- log(Linf[i])
estLinf[i] ~ dlnorm(LgLinf[i], tau.Linf) 


r[i] ~ dnorm(est.r[i], tau.r)
est.r[i] <- bTmax[group[i]]/estTmax[i]+pow(estLinf[i],bLinf[group[i]])

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
win.data <- list("r","Tmax","Linf","group","n")

# Function to generate starting values
inits <- function()
  list (
    mu.bTmax=1,tau.bTmax=1,
    mu.bLinf=0,tau.bLinf=1,
    bTmax=c(1,1,1),
    bLinf=c(1,1,1),
    tau.Tmax=1,estTmax=dat$Tmax,
    tau.Linf=1,estLinf=dat$Linf,
    tau.r=1)


# Parameters to be monitored (= to estimate)
params <- c("bTmax","bLinf",
            "mu.bTmax","tau.bTmax","tau.Tmax",
            "mu.bLinf","tau.bLinf","tau.Linf",
            "tau.r","sigma_r","sigma_Tmax","sigma_Linf",
            "bpvalue","predicted", 
            "estTmax","estLinf","residual")

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
                           out$sims.list$bLinf))



dic = out$DIC
pd = out$pD

# require(devEMF)
# emf("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLinfgroup\\residual plot mix without beta0.emf",
#     width = 6.7, height = 3.35)

pdf('E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\figure\\figure 10.pdf',width=7, height=3.8,onefile = FALSE)

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
LLmat=out$sims.matrix[1:20000,1:176]
dim(LLmat)
rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1:2, each = 10000))
loo_results=loo(LLmat, r_eff = rel_n_eff)
loo_mat=t(loo_results$estimates)
loo_mat

WAIC_mat= waic(LLmat)
waic_mat=t(WAIC_mat$estimates)
waic_mat




library(xlsx)
write.xlsx(stat,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLinfgroup\\mix without beta0.xlsx",sheetName = "stat",row.names=T)
write.xlsx(coe,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLinfgroup\\mix without beta0.xlsx",sheetName = "coe",row.names=T,append=T)
write.xlsx(result,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLinfgroup\\mix without beta0.xlsx",sheetName = "result",row.names=T,append=T)
write.xlsx(loo_mat,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLinfgroup\\mix without beta0.xlsx",sheetName = "loo_estimates",row.names=T,append=T)

# require(devEMF)
emf("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLinfgroup\\mix without beta0.emf",
    width = 3.35, height = 3.35)
observed=r
BHEIV_predicted=out$mean$predicted#BHEIV model predicted r
plot(observed,col="red",pch=1, main = "rTmaxLinfgroup mix without beta0",cex.main=0.8,
     xlab = "index",ylab="r")
points(BHEIV_predicted,col="blue",pch=2)
legend("topright",legend = c("observed r","BHEIV_predicted r"),col=c("red","blue"),pch=c(1,2),cex=0.7)
dev.off()


r_com=data.frame(observed,BHEIV_predicted)
write.xlsx(r_com,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLinfgroup\\mix without beta0.xlsx",sheetName = "r_com",row.names=T,append=T)



est_Tmax=out$mean$estTmax
est_Linf=out$mean$estLinf
p_com=data.frame(Tmax,est_Tmax,
                 Linf,est_Linf)
write.xlsx(p_com,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmaxLinfgroup\\mix without beta0.xlsx",sheetName = "p_com",row.names=T,append=T)




#！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
data=data.frame(out$sims.list$bTmax);names(data)=c("b1Tmax","b2Tmax","b3Tmax")
library(ggplot2)
cols=c('Invertebrate'="#f04546", 'Elasmobranch'="#33CC33", 'Teleost'="#3591d1")
shapes=c('Invertebrate'=1,'Elasmobranch'=3,'Teleost'=5)
density_plot1=
  ggplot(data)+
  geom_density(aes(b1Tmax,color="Invertebrate",linetype="Invertebrate"),cex=0.8)+
  geom_density(aes(b2Tmax,color="Elasmobranch",linetype="Elasmobranch"),cex=0.8)+
  geom_density(aes(b3Tmax,color="Teleost",linetype="Teleost"),cex=0.8)+
  labs(y="Density",x=expression(italic(beta[T[max]])))+
  theme_bw()+
  scale_colour_manual(name="group",values=cols)+
  scale_linetype_manual(name="group",values = shapes)+
  # theme(legend.position=c(.4, .99),legend.justification = c("right", "top"))+
  # theme(legend.key = element_blank())+theme(legend.background = element_blank())+
  xlim(-2,5)
plot(density_plot1)



data=data.frame(out$sims.list$bLinf);names(data)=c("b1Linf","b2Linf","b3Linf")
library(ggplot2)
cols=c('Invertebrate'="#f04546", 'Elasmobranch'="#33CC33", 'Teleost'="#3591d1")
shapes=c('Invertebrate'=1,'Elasmobranch'=3,'Teleost'=5)
density_plot2=
  ggplot(data)+
  geom_density(aes(b1Linf,color="Invertebrate",linetype="Invertebrate"),cex=0.8)+
  geom_density(aes(b2Linf,color="Elasmobranch",linetype="Elasmobranch"),cex=0.8)+
  geom_density(aes(b3Linf,color="Teleost",linetype="Teleost"),cex=0.8)+
  labs(y="Density",x=expression(italic(beta[L[inf]])))+
  theme_bw()+
  scale_colour_manual(name="group",values=cols)+
  scale_linetype_manual(name="group",values = shapes)+
  # theme(legend.position=c(0.45, .99),legend.justification = c("right", "top"))+
  # theme(legend.key = element_blank())+theme(legend.background = element_blank())+
  xlim(-1.5,1)
plot(density_plot2)


library(magrittr)
library(ggpubr)
pdf('E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\figure\\figure 11.pdf',width=8, height=4,onefile = FALSE)
# png( paste0('E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\figure\\figure 1.png'), width= 20, height=12, unit='cm', res=300)
ggarrange(density_plot1,density_plot2,labels = c("A", "B"), 
          font.label=list(size=10,face="plain"),ncol = 2, nrow =1,common.legend = TRUE, legend = "top")

dev.off()