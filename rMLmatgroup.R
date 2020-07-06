

rm(list = ls())

library(lattice)
library("R2WinBUGS")	# Load the R2WinBUGS library
library("lme4")

### 1. Data generation
# Generate two samples of body mass measurements of male peregrines
dat=read.csv("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\data\\162data\\rdata162.csv")
group=dat$group
r=dat$r
M=dat$M
Lmat=dat$Lmat
n=162

### 2.Analysis using WinBUGS
#R work directory
setwd("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLmatgroup") # May have to adapt that

# Save BUGS description of the model to working directory
sink("model.txt")
cat("
model {
# hyper-prior
mu.bM ~ dnorm(0, 0.01)
tau.bM~ dgamma (0.5, 0.01)	

mu.bLmat ~ dnorm(0, 0.1)
tau.bLmat~ dgamma (0.5,0.01)

#prior
for(gp in 1:3){
bM[gp]~dnorm(mu.bM, tau.bM)
bLmat[gp]~dnorm(mu.bLmat, tau.bLmat)
}

tau.M ~ dgamma(0.01,0.01)
tau.Lmat ~ dgamma(5,.01) 
tau.r~dgamma(0.01,0.01)

# Likelihood	
for (i in 1:n) 	{

LgM[i]<- log(M[i])
estM[i] ~ dlnorm(LgM[i], tau.M)

LgLmat[i]<- log(Lmat[i])
estLmat[i] ~ dlnorm(LgLmat[i], tau.Lmat)

r[i] ~ dnorm(est.r[i], tau.r)
est.r[i] <- bM[group[i]]*estM[i]+pow(estLmat[i],bLmat[group[i]])



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
win.data <- list("r","M","Lmat","n","group")

# Function to generate starting values
inits <- function()
  list (
    mu.bM=1,tau.bM=1,
    mu.bLmat=-2,tau.bLmat=1,
    bM=c(1,1,1),  bLmat=c(-2,-2,-2),
    tau.M=1, tau.Lmat=1, tau.r=1,
    estM=dat$M, estLmat=dat$Lmat, 
    r.new=dat$r)


# Parameters to be monitored (= to estimate)

params <- c("bM", "bLmat",
            "mu.bM","tau.bM","tau.M",
            'mu.bLmat','tau.bLmat', 'tau.Lmat',
            "tau.r","bpvalue","predicted", 
            "estM","estLmat","est.r","residual")

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
                   summary(out$sims.list$bLmat))


dic = out$DIC
pd = out$pD

require(devEMF)
emf("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLmatgroup\\residual plot mix without beta0.emf",
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
LLmat=out$sims.matrix[1:20000,1:177]
dim(LLmat)
rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1:2, each = 10000))
loo_results=loo(LLmat, r_eff = rel_n_eff)
loo_mat=loo_results$estimates
loo_mat

WAIC_mat= waic(LLmat)
waic_mat=WAIC_mat$estimates
waic_mat


library(xlsx)
write.xlsx(stat,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLmatgroup\\mix without beta0.xlsx",sheetName = "stat",row.names=T)
write.xlsx(coe,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLmatgroup\\mix without beta0.xlsx",sheetName = "coe",row.names=T,append=T)
write.xlsx(result,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLmatgroup\\mix without beta0.xlsx",sheetName = "result",row.names=T,append=T)
write.xlsx(waic_mat,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLmatgroup\\mix without beta0.xlsx",sheetName = "waic_mat",row.names=T,append=T)



#！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
# require(devEMF)
emf("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLmatgroup\\mix without beta0.emf",
    width = 3.35, height = 3.35)
observed=r
BHEIV_predicted=out$mean$predicted#BHEIV model predicted r
plot(observed,col="red",pch=1, main = "rMLmatgroup mix without beta0",cex.main=0.8,
     xlab = "index",ylab="r")
points(BHEIV_predicted,col="blue",pch=2)
legend("topright",legend = c("observed r","BHEIV_predicted r"),col=c("red","blue"),pch=c(1,2),cex=0.7)
dev.off()


sigma_r=sqrt(1/out$mean$tau.r)
r_com=data.frame(observed,BHEIV_predicted,sigma_r)
write.xlsx(r_com,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLmatgroup\\mix without beta0.xlsx",sheetName = "r_com",row.names=T,append=T)


est_M=out$mean$estM
sigma_M=sqrt(1/out$mean$tau.M)

est_Lmat=out$mean$estLmat
sigma_Lmat=sqrt(1/out$mean$tau.Lmat)

p_com=data.frame(M,est_M,sigma_M,
                 Lmat,est_Lmat,sigma_Lmat)
write.xlsx(p_com,"E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLmatgroup\\mix without beta0.xlsx",sheetName = "p_com",row.names=T,append=T)




# 
# # # # BGR convergence diagnostic called Rhat, and an
# # #BGR喘栖登僅頁倦辺漸��自峙頁1.1匆辛恬葎辺漸
# # hist(out$summary[,8])			# Rhat values in the eighth column of the summary
# # which(out$summary[,8] > 1.1)		# None in this case
# # # # # For trace plots for the entire chains, do
# # matplot(out$sims.array[1000:5000,1:3,1], type = "l")
# # # # We can also produce graphical summaries, e.g., histograms of the posterior
# # # # distributions for each parameter:
# 
# data1=data.frame(out$sims.list$bM);names(data1)=c("b1M","b2M","b3M")
# data2=data.frame(out$sims.list$bLmat);names(data2)=c("b1Lmat","b2Lmat","b3Lmat")
# 
# 
# 
# library(ggplot2)
# cols=c('Invertebrate'="#f04546", 'Elasmobranch'="#33CC33", 'Teleost'="#3591d1")
# shapes=c('Invertebrate'=1,'Elasmobranch'=3,'Teleost'=5)
# bM_plot=
#   ggplot(data1)+
#   geom_density(aes(b1M,color="Invertebrate",linetype="Invertebrate"),cex=0.8)+
#   geom_density(aes(b2M,color="Elasmobranch",linetype="Elasmobranch"),cex=0.8)+
#   geom_density(aes(b3M,color="Teleost",linetype="Teleost"),cex=0.8)+
#   labs(y="Density",x=expression(italic(beta[M])))+
#   theme_bw()+
#   scale_colour_manual(name="group",values=cols)+
#   scale_linetype_manual(name="group",values = shapes)
# # theme(legend.position=c(.99, .99),legend.justification = c("right", "top"))+
# # theme(legend.key = element_blank())+theme(legend.background = element_blank())
# plot(bM_plot)
# 
# 
# bLmat_plot=
#   ggplot(data2)+
#   geom_density(aes(b1Lmat,color="Invertebrate",linetype="Invertebrate"),cex=0.8)+
#   geom_density(aes(b2Lmat,color="Elasmobranch",linetype="Elasmobranch"),cex=0.8)+
#   geom_density(aes(b3Lmat,color="Teleost",linetype="Teleost"),cex=0.8)+
#   labs(y="Density",x=expression(italic(beta[L[mat]])))+
#   theme_bw()+
#   scale_colour_manual(name="group",values=cols)+
#   scale_linetype_manual(name="group",values = shapes)+
#   xlim(-0.001,0.004)
# # theme(legend.position=c(.99, .99),legend.justification = c("right", "top"))+
# # theme(legend.key = element_blank())+theme(legend.background = element_blank())+
# plot(bLmat_plot)
# 
# 
# library(magrittr)
# library(ggpubr)
# ggarrange(bM_plot,bLmat_plot,labels = c("A", "B"), font.label=list(size=10,face="plain"),
#           ncol = 2, nrow = 1,common.legend = TRUE,legend="top")
# 
# 
# pdf("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLmatgroup\\density plot1.pdf",
#     width = 7, height = 3.5,onefile = FALSE)
# ggarrange(bM_plot,bLmat_plot,labels = c("A", "B"), font.label=list(size=10,face="plain"),
#           ncol = 2, nrow = 1,common.legend = TRUE,legend="top")
# dev.off()
# 
# 
# #！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
# r_predicted_mean=out$mean$predicted
# ob_r=dat$r
# rdata=data.frame(r_predicted_mean,ob_r)
# library(ggplot2)
# p1=ggplot(rdata,aes(x=ob_r,y=r_predicted_mean))+
#   geom_point(size=2,shape=1)+
#   geom_abline(slope = 1,intercept = 0,color="red")+
#   theme_bw()+
#   labs(y=expression(paste("Predicted mean of ",r)),x=expression(paste("Observed ",r)),
#        title=expression(paste("BHEIV model M19")))+
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
#        title=expression(paste("BHEIV model M19")))+
#   theme(plot.title = element_text(hjust = 0.5,size=10))
# plot(p2)
# 
# library(magrittr)
# library(ggpubr)
# ggarrange(p1,p2,ncol = 2, nrow = 1,labels = c("A", "B"),font.label=list(size=10,face="plain"))
# 
# pdf("E:\\01-doctor_LCJ\\09-fishery\\02-r_LHPs paper\\20200312BHEIV\\rMLmatgroup\\residual plot.pdf",
#     width = 7, height = 3.5,onefile = FALSE)
# ggarrange(p1,p2,ncol = 2, nrow = 1,labels = c("A", "B"),font.label=list(size=10,face="plain"))
# 
# dev.off()
# 
