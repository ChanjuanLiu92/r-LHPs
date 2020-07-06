
rm(list=ls())

library(lattice)
library("R2WinBUGS")	# Load the R2WinBUGS library
library("lme4")

### 1. Data generation
# Generate two samples of body mass measurements of male peregrines
dat=read.csv("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\data\\162data\\rdata162.csv")
r=dat$r
Tmat=dat$Tmat
group=dat$group
n=162

### 2.Analysis using WinBUGS
#R work directory
setwd("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmatgroup") # May have to adapt that

# Save BUGS description of the model to working directory
sink("model.txt")
cat("
model {
# hyper-prior
mu.bTmat ~ dnorm(0, 0.01)
tau.bTmat~ dgamma (0.5,0.01)

#prior
for(gp in 1:3){
bTmat[gp]~dnorm(mu.bTmat, tau.bTmat)
}

tau.Tmat ~ dgamma(0.01,0.01)
tau.r~dgamma(0.01,0.01)


# Likelihood	
for (i in 1:n) 	{
LgTmat[i]<- log(Tmat[i])
estTmat[i] ~ dlnorm(LgTmat[i], tau.Tmat) 

r[i] ~ dnorm(est.r[i], tau.r)
est.r[i] <- bTmat[group[i]]/estTmat[i]


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
win.data <- list("r","Tmat","group","n")

# Function to generate starting values
inits <- function()
  list (
    mu.bTmat=1,tau.bTmat=1,
    bTmat=c(1,1,1),
    tau.Tmat=1,tau.r=1,
    estTmat=dat$Tmat)


# Parameters to be monitored (= to estimate)
params <- c("bTmat","mu.bTmat","tau.bTmat","tau.Tmat","tau.r","bpvalue","predicted",
            "estTmat","residual")

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

stat=data.frame(out$pD,out$DIC,out$mean$bpvalue)
coe=out$summary
coe_stat=as.matrix(summary(out$sims.list$bTmat))

dic = out$DIC
pd = out$pD

require(devEMF)
emf("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmatgroup\\residual plot inverse function.emf",
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
# LLmat=out$sims.matrix
# dim(LLmat)
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
write.xlsx(stat,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmatgroup\\inverse function.xlsx",sheetName = "stat",row.names=T)
write.xlsx(coe,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmatgroup\\inverse function.xlsx",sheetName = "coe",row.names=T,append=T)
write.xlsx(result,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmatgroup\\inverse function.xlsx",sheetName = "result",row.names=T,append=T)
write.xlsx(waic_mat,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmatgroup\\inverse function.xlsx",sheetName = "waic_mat",row.names=T,append=T)




# ！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
# require(devEMF)
emf("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmatgroup\\inverse function.emf",width = 3.35, height = 3.35)
observed=r
BHEIV_predicted=out$mean$predicted#BHEIV model predicted r
plot(observed,col="red",pch=1, main = "rTmatgroup inverse function",cex.main=0.8,
     xlab = "index",ylab="r")
points(BHEIV_predicted,col="blue",pch=2)
legend("topright",legend = c("observed r","BHEIV_predicted r"),col=c("red","blue"),pch=c(1,2),cex=0.7)
dev.off()


sigma_r=sqrt(1/out$mean$tau.r)
r_com=data.frame(observed,BHEIV_predicted,sigma_r)
write.xlsx(r_com,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmatgroup\\inverse function.xlsx",
           sheetName = "r_com",row.names=T,append=T)



est_Tmat=out$mean$estTmat
sigma_Tmat=sqrt(1/out$mean$tau.Tmat)
Tmat_com=data.frame(Tmat,est_Tmat,sigma_Tmat)
write.xlsx(Tmat_com,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rTmatgroup\\inverse function.xlsx",
           sheetName = "Tmat_com",row.names=T,append=T)





# # # BGR convergence diagnostic called Rhat, and an
# #BGR喘栖登僅頁倦辺漸��自峙頁1.1匆辛恬葎辺漸
# hist(out$summary[,8])			# Rhat values in the eighth column of the summary

data=data.frame(out$sims.list$bTmat);names(data)=c("b1Tmat","b2Tmat","b3Tmat")
library(ggplot2)
cols=c('Invertebrate'="#f04546", 'Elasmobranch'="#33CC33", 'Teleost'="#3591d1")
shapes=c('Invertebrate'=1,'Elasmobranch'=3,'Teleost'=5)
density_plot=
  ggplot(data)+
  geom_density(aes(b1Tmat,color="Invertebrate",linetype="Invertebrate"),cex=0.8)+
  geom_density(aes(b2Tmat,color="Elasmobranch",linetype="Elasmobranch"),cex=0.8)+
  geom_density(aes(b3Tmat,color="Teleost",linetype="Teleost"),cex=0.8)+
  labs(y="Density",x=expression(italic(beta[T[mat]])))+
  theme_bw()+
  scale_colour_manual(name="group",values=cols)+
  scale_linetype_manual(name="group",values = shapes)+
  theme(legend.position=c(0.4, .99),legend.justification = c("right", "top"))+
  theme(legend.key = element_blank())+theme(legend.background = element_blank())+
  xlim(-4,2)
plot(density_plot)



# emf("G:\\02-r_LHPs paper\\results20200222\\rTmatgroup\\density_plot nlm without beta0.emf",width = 3.35, height = 2.54)
# plot(density_plot)
# dev.off()

# 
# 
# 
# r_predicted_mean=exp(out$mean$predicted)
# ob_r=dat$r
# rdata=data.frame(r_predicted_mean,ob_r)
# library(ggplot2)
# p1=ggplot(rdata,aes(x=ob_r,y=r_predicted_mean))+
#   geom_point(size=2,shape=1)+
#   geom_abline(slope = 1,intercept = 0,color="red")+
#   theme_bw()+
#   labs(y=expression(paste("Predicted mean of ",r)),x=expression(paste("Observed ",r)),
#        title=expression(paste("BHEIV model M(11,1)")))+
#   theme(plot.title = element_text(hjust = 0.5,size=10))
# plot(p1)
# 
# residuals=ob_r-r_predicted_mean
# residual_data=data.frame(residuals,r_predicted_mean)
# p2=ggplot(residual_data,aes(x=r_predicted_mean,y=residuals))+
#   ylim(-1,1)+
#   geom_point(size=2,shape=1)+
#   geom_hline(yintercept=0,color="red")+
#   theme_bw()+
#   labs(y="Residuals",x=expression(paste("Predicted mean of ",r)),
#        title=expression(paste("BHEIV model M(11,1)")))+
#   theme(plot.title = element_text(hjust = 0.5,size=10))
# plot(p2)
# 
# library(magrittr)
# library(ggpubr)
# ggarrange(p1,p2,ncol = 2, nrow = 1,labels = c("A", "B"),font.label=list(size=10,face="plain"))
# 
