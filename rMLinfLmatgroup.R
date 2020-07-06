rm(list=ls())
library(lattice)
library("R2WinBUGS")	# Load the R2WinBUGS library
library("lme4")

### 1. Data generation
# Generate two samples of body mass measurements of male peregrines
dat=read.csv("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\data\\162data\\rdata162.csv")
r=dat$r
M=dat$M
k=dat$k
Linf=dat$Linf
Lmat=dat$Lmat
group=dat$group
n=162

### 2.Analysis using WinBUGS
setwd("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLinfLmatgroup") # May have to adapt that

# Save BUGS description of the model to working directory
sink("model.txt")
cat("
model {
# hyper-prior
mu.bM ~ dnorm(0, 0.01)
tau.bM~ dgamma (0.5, 0.01)

mu.bLinf ~ dnorm(0, 0.1)
tau.bLinf~ dgamma (0.5,0.01)

mu.bLmat ~ dnorm(0, 0.1)
tau.bLmat~ dgamma (0.5,0.01)


#prior
for(gp in 1:3){
bM[gp]~dnorm(mu.bM, tau.bM)
bLinf[gp]~dnorm(mu.bLinf, tau.bLinf)
bLmat[gp]~dnorm(mu.bLmat, tau.bLmat)
}


tau.M ~ dgamma(0.01,0.01)
tau.Linf ~ dgamma(50,.1)
tau.Lmat ~ dgamma(50,.1)
tau.r~dgamma(0.01,0.01)


# Likelihood	
for (i in 1:n) 	{
LgM[i]<- log(M[i])
estM[i] ~ dlnorm(LgM[i], tau.M)

LgLinf[i]<- log(Linf[i])
estLinf[i] ~ dlnorm(LgLinf[i], tau.Linf)

LgLmat[i]<- log(Lmat[i])
estLmat[i] ~ dlnorm(LgLmat[i], tau.Lmat)

r[i] ~ dnorm(est.r[i], tau.r)


est.r[i] <- bM[group[i]]*estM[i]+
pow(estLinf[i],bLinf[group[i]])+pow(estLmat[i],bLmat[group[i]])


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
win.data <- list("r","Linf","M","Lmat","group","n")

# Function to generate starting values
inits <- function()
  list (
    mu.bM=0,tau.bM=0.01,
    mu.bLinf=0,tau.bLinf=0.01,
    mu.bLmat=0,tau.bLmat=0.01,
    bM=c(1,1,1),bLinf=c(1,1,1),bLmat=c(1,1,1),
    tau.M=0.01,estM=dat$M,
    tau.Linf=0.01,estLinf=dat$Linf,
    tau.Lmat=0.01,estLmat=dat$Lmat,
    tau.r=0.01)


# Parameters to be monitored (= to estimate)
params <- c("bM","bLinf","bLmat","predicted","bpvalue")


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

stat =data.frame(out$pD,out$DIC,out$mean$bpvalue)
coe=out$summary
coe_stat=as.matrix(summary(out$sims.list$bM),
                   summary(out$sims.list$bLinf),
                   ummary(out$sims.list$bLmat))


dic = out$DIC
pd = out$pD

plot(out$median$residual~ dat$r)
points(out$mean$residual~ dat$r, pch='x', col=2)
plot(out$median$predict~ dat$r) ; abline(0,1)


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
rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1:2, each = 10000))
loo_results=loo(LLmat, r_eff = rel_n_eff)
loo_mat=loo_results$estimates
loo_mat

WAIC_mat= waic(LLmat)
waic_mat=WAIC_mat$estimates
waic_mat


library(xlsx)
write.xlsx(stat,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLinfLmatgroup\\original data without beta0-1.xlsx",sheetName = "stat",row.names=T)
write.xlsx(coe,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLinfLmatgroup\\original data without beta0-1.xlsx",sheetName = "coe",row.names=T,append=T)
write.xlsx(result,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLinfLmatgroup\\original data without beta0-1.xlsx",sheetName = "result",row.names=T,append=T)
write.xlsx(loo_mat,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLinfLmatgroup\\original data without beta0-1.xlsx",sheetName = "loo_estimates",row.names=T,append=T)


#！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
# require(devEMF)
emf("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLinfLmatgroup\\original data without beta0-1.emf",
    width = 3.35, height = 3.35)
observed=r
BHEIV_predicted=out$mean$predicted#BHEIV model predicted r
plot(observed,col="red",pch=1, main = "rMLinfLmatgroup original data without beta0-1",cex.main=0.8,
     xlab = "index",ylab="r")
points(BHEIV_predicted,col="blue",pch=2)
legend("topright",legend = c("observed r","BHEIV_predicted r"),col=c("red","blue"),pch=c(1,2),cex=0.7)
dev.off()


sigma_r=sqrt(1/out$mean$tau.r)
r_com=data.frame(observed,BHEIV_predicted,sigma_r)
write.xlsx(r_com,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLinfLmatgroup\\original data without beta0-1.xlsx",sheetName = "r_com",row.names=T,append=T)


est_M=out$mean$estM
sigma_M=sqrt(1/out$mean$tau.M)

est_Linf=out$mean$estLinf
sigma_Linf=sqrt(1/out$mean$tau.Linf)

est_Lmat=out$mean$estLmat
sigma_Lmat=sqrt(1/out$mean$tau.Lmat)

p_com=data.frame(M,est_M,sigma_M,
                 Linf,est_Linf,sigma_Linf,
                 Lmat,est_Lmat,sigma_Lmat)
write.xlsx(p_com,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLinfLmatgroup\\original data without beta0-1.xlsx",sheetName = "p_com",row.names=T,append=T)


# 
# # # BGR convergence diagnostic called Rhat, and an
# #BGR喘栖登僅頁倦辺漸��自峙頁1.1匆辛恬葎辺漸
# hist(out$summary[,8])			# Rhat values in the eighth column of the summary
# which(out$summary[,8] > 1.1)		# None in this case
# # # # For trace plots for the entire chains, do
# matplot(out$sims.array[1000:5000,1:3,1], type = "l")
# # # # We can also produce graphical summaries, e.g., histograms of the posterior
# # # # distributions for each parameter:
# data1=data.frame(out$sims.list$bM);names(data1)=c("b1M","b2M","b3M")
# data2=data.frame(out$sims.list$bk);names(data2)=c("b1k","b2k","b3k")
# data3=data.frame(out$sims.list$bTmat);names(data3)=c("b1Tmat","b2Tmat","b3Tmat")
# 
# library(ggplot2)
# cols=c('Invertebrate'="#f04546", 'Elasmobranch'="#33CC33", 'Teleost'="#3591d1")
# shapes=c('Invertebrate'=1,'Elasmobranch'=3,'Teleost'=5)
# bM_plot=
#   ggplot(data1)+
#   geom_density(aes(b1M,color="Invertebrate",linetype="Invertebrate"))+
#   geom_density(aes(b2M,color="Elasmobranch",linetype="Elasmobranch"))+
#   geom_density(aes(b3M,color="Teleost",linetype="Teleost"))+
#   labs(y="Density",x=expression(italic(beta[M])))+
#   theme_bw()+
#   scale_colour_manual(name="group",values=cols)+
#   scale_linetype_manual(name="group",values = shapes)+
#   theme(legend.position=c(.99, .99),legend.justification = c("right", "top"))+
#   theme(legend.key = element_blank())+theme(legend.background = element_blank())
# plot(bM_plot)
# 
# library(devEMF)
# emf("G:\\02-r_LHPs paper\\results20200222\\rMkTmatgroup\\bM_plot.emf",width = 3.35, height = 2.54)
# plot(bM_plot)
# dev.off()
# 
# 
# bk_plot=
#   ggplot(data2)+
#   geom_density(aes(b1k,color="Invertebrate",linetype="Invertebrate"))+
#   geom_density(aes(b2k,color="Elasmobranch",linetype="Elasmobranch"))+
#   geom_density(aes(b3k,color="Teleost",linetype="Teleost"))+
#   labs(y="Density",x=expression(italic(beta[k])))+
#   theme_bw()+
#   scale_colour_manual(name="group",values=cols)+
#   scale_linetype_manual(name="group",values = shapes)+
#   theme(legend.position=c(.99, .99),legend.justification = c("right", "top"))+
#   theme(legend.key = element_blank())+theme(legend.background = element_blank())
# plot(bk_plot)
# 
# library(devEMF)
# emf("G:\\02-r_LHPs paper\\results20200222\\rMkTmatgroup\\bk_plot.emf",width = 3.35, height = 2.54)
# plot(bk_plot)
# dev.off()
# 
# bTmat_plot=
#   ggplot(data3)+
#   geom_density(aes(b1Tmat,color="Invertebrate",linetype="Invertebrate"))+
#   geom_density(aes(b2Tmat,color="Elasmobranch",linetype="Elasmobranch"))+
#   geom_density(aes(b3Tmat,color="Teleost",linetype="Teleost"))+
#   labs(y="Density",x=expression(italic(beta[T[mat]])))+
#   theme_bw()+
#   scale_colour_manual(name="group",values=cols)+
#   scale_linetype_manual(name="group",values = shapes)+
#   theme(legend.position=c(.99, .99),legend.justification = c("right", "top"))+
#   theme(legend.key = element_blank())+theme(legend.background = element_blank())
# plot(bTmat_plot)
# 
# library(devEMF)
# emf("G:\\02-r_LHPs paper\\results20200222\\rMkTmatgroup\\bTmat_plot.emf",width = 3.35, height = 2.54)
# plot(bTmat_plot)
# dev.off()
# 
# 
# # pdf("G:\\02-r_LHPs paper\\results20200222\\rMgroup\\bTmat_plot.pdf",width = 3.35, height = 2.54)
# # plot(density_plot)
# # dev.off()
# 
# # tiff("G:\\02-r_LHPs paper\\results20200222\\rMgroup\\density_plot.tif",width = 335, height = 254,compression = c("none"))
# # png("G:\\02-r_LHPs paper\\results20200222\\rMgroup\\density_plot.png",width = 335, height = 254)
# 
# 
