rm(list=ls())
library(survival)
library(MASS)
library(matrixcalc)

args <- commandArgs(trailingOnly = TRUE)
#setwd("C:/Users/ycho/Box/Yeying/Mediation analysis/survival mediaion analysis")
#setwd("C:/Users/ycho/Box/Current Projects (Methods)/Mediation analysis/Mediation analysis/survival mediaion analysis")

n <- 500
nrun <- 1
mm <- 50
cen.rate <- rep(0,nrun)
nde_lasso <- rep(0,nrun)
nie.prod_lasso <- rep(0,nrun)
nie.diff_lasso <- rep(0,nrun)

nde.alasso <- rep(0,nrun)
nde.alasso.se <- rep(0,nrun)

nie.prod_alasso <- rep(0,nrun)
nie.diff_alasso <- rep(0,nrun)
nie.prod.alasso.se <- rep(0,nrun)
cover.nde.alasso <- rep(0,nrun)
cover.nie.prod.alasso <- rep(0,nrun)
beta_lasso_m <- matrix(0,nrow=nrun,ncol=mm+3)
sigma2_lasso_m <- rep(0,nrun)

beta_alasso_m <- matrix(0,nrow=nrun,ncol=mm+3)
sigma2_alasso_m <- rep(0,nrun)
se_beta_alasso <- matrix(0,nrow=nrun,ncol=mm+3)
se_sigma2_alasso <- rep(0,nrun) 

beta_scad_m <- matrix(0,nrow=nrun,ncol=mm+3)
sigma2_scad_m <- rep(0,nrun)
se_beta_scad <- matrix(0,nrow=nrun,ncol=mm+3)
se_sigma2_scad <- rep(0,nrun) 

sel_scad_m <- matrix(0, ncol=mm, nrow=nrun)
sel_lasso_m <- matrix(0, ncol=mm, nrow=nrun)
sel_alasso_m <- matrix(0, ncol=mm, nrow=nrun)

nde.scad <- rep(0,nrun)
nde.scad.se <- rep(0,nrun)
nie.prod_scad <- rep(0,nrun)
nie.prod.scad.se <- rep(0,nrun)
nie.diff_scad <- rep(0,nrun)
cover.nde.scad <- rep(0,nrun)
cover.nie.prod.scad <- rep(0,nrun)

seednum <- as.integer(args[1])
	set.seed(1+seednum)
 	rho = 0.5
 	x <- rnorm(n,0,1)
 	p <- exp(-0.2+0.5*x)/(1 + exp(-0.2+0.5*x))
 	A <- rbinom(n,1,p)
	alpha10 <- 0
	alpha_1a <- -0.5
	alpha20 <- 0.2
	alpha_2a <- 0.3
	alpha_1b <- -0.4
	alpha_2b <- 0.5
	#M <- rnorm(n,mean=alpha10+alpha_1a*A,sd=1)
	#Sig <- matrix(c(1,0,0,1),ncol=2,nrow=2)
	Sig <- matrix(0,ncol=mm,nrow=mm)
	for(j in 1:mm){
		for(l in 1:mm){
			Sig[j,l] <- rho^abs(j-l)
		}
	}
	alpha1 <- alpha10 + alpha_1a*A + alpha_1b*x
	alpha2 <- alpha20 + alpha_2a*A + alpha_2b*x
	Mno.temp <- mvrnorm(n,mu = rep(0,mm),Sig=Sig)
	Mno <- cbind(Mno.temp[,1] + alpha1,Mno.temp[,2] + alpha2,Mno.temp[,3:mm])
	#M1 <- eps[,1] + alpha1; M2 <- eps[,2] + alpha2
	M1 <- Mno[,1]
	M2 <- Mno[,2]
	
	beta0 <- 0
	beta_a <- 1
	beta_m1 <- 0.5
	beta_m2 <- 0.3
	beta_x <- -0.2
	Time <- exp(rnorm(n,mean=beta0 + beta_a*A+beta_m1*M1 + beta_m2*M2 + beta_x*x,sd=1))
	Cen <- runif(n,min=0,max=5.8)
	obs <- pmin(Time,Cen)
	status <- as.numeric(Time <= Cen)
      cen.rate <- mean(status==0)
	Mmat <- cbind(1,A,Mno,x)
	p <- dim(Mmat)[2]
	# fit model
	fit.all <- survreg(Surv(obs,status) ~ -1 + Mmat, dist = "lognormal")
	source("lasso_code0622.R")
    tun_lasso <- c(tun_v[which.min(bic_v)])	 
    beta_lasso<-c(beta_m[which.min(bic_v),])
    sigma2_lasso<-c(sigma2_v[which.min(bic_v)])
    sel_lasso <-  which(abs(beta_lasso[3:(p-1)]) > 1e-3)
    
    Mmat_lasso <- Mno[,sel_lasso]
    #sel_lasso_m <- abs(beta_lasso[3:(p-1)]) > 1e-3
    sel_lasso_m[sel_lasso]<- ifelse(sel_lasso >=1, 1, 0)

    # reduced model
    exp.model <- survreg(Surv(obs,status) ~ A, dist = "lognormal")
    
    # Run mediator model
    med.model_lasso <- lm(Mmat_lasso ~ A + x)

    nde_lasso <- beta_lasso[2]
    
    sel_est_lasso <- sel_lasso + 2
    
    if(length(sel_lasso) > 1){
      nie.prod_lasso <- med.model_lasso$coefficients[2,]%*%beta_lasso[sel_est_lasso]
    }else{
      nie.prod_lasso <- med.model_lasso$coefficients[2] * beta_lasso[sel_est_lasso]
    }
    nie.diff_lasso <- exp.model$coefficients[2]-beta_lasso[2]
    beta_lasso_m <- beta_lasso
    sigma2_lasso_m <- sigma2_lasso

	source("alasso_code0622.R")
    iter<-iter_v[which.min(bic_v)]
    tun_alasso<-tun_v[which.min(bic_v)]
    beta_alasso<-c(beta_m[which.min(bic_v),])
    sigma2_alasso<-c(sigma2_v[which.min(bic_v)])
    SE_beta_alasso<-c(SE_beta_m[which.min(bic_v),])
    SE_sigma2_alasso<-c(SE_sigma2_v[which.min(bic_v)])
    alassocov <- cov.arr[,,which.min(bic_v)]
    bic<-min(bic_v)
    
  	# Est.
    beta_alasso_m<-beta_alasso  # beta_hat
    sigma2_alasso_m<-sigma2    # sigma2_hat
    
    	# SE
    se_beta_alasso<-SE_beta_alasso             # beta_hat
    se_sigma2_alasso<-SE_sigma2_alasso       

    sel_alasso <-  which(abs(beta_alasso[3:(p-1)]) > 1e-3)
    
    Mmat_alasso <- Mno[,sel_alasso]
 #   sel_alasso_m<- abs(beta_lasso[3:(p-1)]) > 1e-3
    
    sel_alasso_m[sel_alasso]<- ifelse(sel_alasso >=1, 1, 0)
    med.model_alasso <- lm(Mmat_alasso ~ A + x)

    nde.alasso <- beta_alasso[2]
    
    sel_est_alasso <- sel_alasso + 2
    if(length(sel_alasso) > 1){
      	nie.prod_alasso <- med.model_alasso$coefficients[2,]%*%beta_alasso[sel_est_alasso]
		SSCP.E <- crossprod(Mmat_alasso - med.model_alasso$fitted.values)
		MCP.E <- SSCP.E/(n-3)
		cov.mat.m <- kronecker(MCP.E,summary(med.model_alasso)[[1]]$cov.unscaled)
    	}else{
      	nie.prod_alasso <- med.model_alasso$coefficients[2] * beta_alasso[sel_est_alasso]
		SSCP.E <- sum((Mmat_alasso - med.model_alasso$fitted.values)^2)
		MCP.E <- SSCP.E/(n-3)
		cov.mat.m <- as.numeric(MCP.E)*summary(med.model_alasso)$cov.unscaled
    	}
    nie.diff_alasso <- exp.model$coefficients[2]-beta_alasso[2]
    
    
	if(length(sel_alasso) > 1){
		med.model.coef <- rep(0,dim(Mmat_alasso)[2])
		med.model.se <- rep(0,dim(Mmat_alasso)[2])
		for(j in 1:dim(Mmat_alasso)[2]){
			med.model.coef[j] <- as.numeric(summary(med.model_alasso)[[j]]$coefficients[2,1])
			med.model.se[j] <-   as.numeric(summary(med.model_alasso)[[j]]$coefficients[2,2])
		}
	}else{
		med.model.coef <- as.numeric(summary(med.model_alasso)$coefficients[2,1])
		med.model.se <-   as.numeric(summary(med.model_alasso)$coefficients[2,2])
	
	}
	alassoselcov <- alassocov[sel_est_alasso,sel_est_alasso]
	beta_alasso_sel <- beta_alasso[sel_est_alasso]
	len_b <- length(sel_alasso)
	if(length(sel_alasso) > 1){
		cov.mat.m.sel.all <- cov.mat.m[seq(2,dim(cov.mat.m)[1],3),seq(2,dim(cov.mat.m)[1],3)]
		cov.mat.m.sel <- cov.mat.m.sel.all[upper.tri(cov.mat.m.sel.all)]
		nde.alasso.se <- sqrt(alassocov[2,2])
		nie.prod.alasso.se <- sqrt(med.model.coef^2 %*% diag(alassoselcov) + beta_alasso_sel^2 %*% med.model.se^2 + 2*sum(combn(med.model.coef,m=2,prod)*alassoselcov[upper.tri(alassoselcov)])+ 
		2*sum(combn(beta_alasso_sel,m=2,prod)* cov.mat.m.sel))
      #nie.prod.se[i] <- as.numeric(sqrt(med.model1.coef^2 * summary(full.model)$var[3,3] + full.model$coefficients[3]^2*summary(med.model)[[1]]$coefficients[2,2]^2 + med.model2.coef^2 * summary(full.model)$var[4,4] + full.model$coefficients[4]^2*summary(med.model)[[2]]$coefficients[2,2]^2 + 2*med.model1.coef*med.model2.coef*summary(full.model)$var[3,4] + 2*full.model$coefficients[3]*full.model$coefficients[4]*cov.mat.m[2,4]))
		cover.nde.alasso <- as.numeric((nde.alasso - qnorm(0.975) * nde.alasso.se <= 1) & (1 <= nde.alasso + qnorm(0.975) * nde.alasso.se))
		cover.nie.prod.alasso <-  as.numeric((nie.prod_alasso - qnorm(0.975) * nie.prod.alasso.se <= -0.16) & (-0.16 <= nie.prod_alasso + qnorm(0.975) * nie.prod.alasso.se))
	}else{
		cov.mat.m.sel <- cov.mat.m[2,2]
		nde.alasso.se <- sqrt(alassocov[2,2])
		nie.prod.alasso.se <- sqrt(med.model.coef^2 * alassoselcov + beta_alasso_sel^2 %*% med.model.se^2)  
		cover.nde.alasso <- as.numeric((nde.alasso - qnorm(0.975) * nde.alasso.se <= 1) & (1 <= nde.alasso + qnorm(0.975) * nde.alasso.se))
		cover.nie.prod.alasso <-  as.numeric((nie.prod_alasso - qnorm(0.975) * nie.prod.alasso.se <= -0.16) & (-0.16 <= nie.prod_alasso + qnorm(0.975) * nie.prod.alasso.se))

	}
    
	source("scad_code0622.R")
    iter<-iter_v[which.min(bic_v)]
    tun_scad<-tun_v[which.min(bic_v)]
    beta_scad<-c(beta_m[which.min(bic_v),])
    sigma2_scad<-c(sigma2_v[which.min(bic_v)])
    SE_beta_scad<-c(SE_beta_m[which.min(bic_v),])
    SE_sigma2_scad<-c(SE_sigma2_v[which.min(bic_v)])
    scadcov <- cov.arr[,,which.min(bic_v)]
    bic<-min(bic_v)
    

    	# Est.
    beta_scad_m<-beta_scad  # beta_hat
    sigma2_scad_m<-sigma2    # sigma2_hat
    
    	# SE
    se_beta_scad<-SE_beta_scad             # beta_hat
    se_sigma2_scad<-SE_sigma2_scad       

    sel_scad <-  which(abs(beta_scad[3:(p-1)]) > 1e-3)
    
    Mmat_scad <- Mno[,sel_scad]
#    sel_scad_m<- abs(beta_scad[3:(p-1)]) > 1e-3
    

    sel_scad_m[sel_scad]<- ifelse(sel_scad >=1, 1, 0)
    med.model_scad <- lm(Mmat_scad ~ A + x)
    
    nde.scad <- beta_scad[2]
    
    sel_est_scad <- sel_scad + 2
    
    if(length(sel_scad) > 1){
      	nie.prod_scad <- med.model_scad$coefficients[2,]%*%beta_scad[sel_est_scad]
		SSCP.E <- crossprod(Mmat_scad - med.model_scad$fitted.values)
		MCP.E <- SSCP.E/(n-3)
		cov.mat.m <- kronecker(MCP.E,summary(med.model_scad)[[1]]$cov.unscaled)
    	}else{
      	nie.prod_scad <- med.model_scad$coefficients[2] * beta_scad[sel_est_scad]
		SSCP.E <- sum((Mmat_scad - med.model_scad$fitted.values)^2)
		MCP.E <- SSCP.E/(n-3)
		cov.mat.m <- as.numeric(MCP.E)*summary(med.model_scad)$cov.unscaled
    	}
    nie.diff_scad <- exp.model$coefficients[2]-beta_scad[2]
    
    

	if(length(sel_scad) > 1){
		med.model.coef <- rep(0,dim(Mmat_scad)[2])
		med.model.se <- rep(0,dim(Mmat_scad)[2])
		for(j in 1:dim(Mmat_scad)[2]){
			med.model.coef[j] <- as.numeric(summary(med.model_scad)[[j]]$coefficients[2,1])
			med.model.se[j] <-   as.numeric(summary(med.model_scad)[[j]]$coefficients[2,2])
		}
	}else{
		med.model.coef <- as.numeric(summary(med.model_scad)$coefficients[2,1])
		med.model.se <-   as.numeric(summary(med.model_scad)$coefficients[2,2])
	
	}
	scadselcov <- scadcov[sel_est_scad,sel_est_scad]
	beta_scad_sel <- beta_scad[sel_est_scad]
	if(length(sel_scad) > 1){
		cov.mat.m.sel.all <- cov.mat.m[seq(2,dim(cov.mat.m)[1],3),seq(2,dim(cov.mat.m)[1],3)]
		cov.mat.m.sel <- cov.mat.m.sel.all[upper.tri(cov.mat.m.sel.all)]
		nde.scad.se <- sqrt(scadcov[2,2])
		nie.prod.scad.se <- sqrt(med.model.coef^2 %*% diag(scadselcov) + beta_scad_sel^2 %*% med.model.se^2 + 2*sum(combn(med.model.coef,m=2,prod)*scadselcov[upper.tri(scadselcov)])+ 
		2*sum(combn(beta_scad_sel,m=2,prod)* cov.mat.m.sel))
      #nie.prod.se[i] <- as.numeric(sqrt(med.model1.coef^2 * summary(full.model)$var[3,3] + full.model$coefficients[3]^2*summary(med.model)[[1]]$coefficients[2,2]^2 + med.model2.coef^2 * summary(full.model)$var[4,4] + full.model$coefficients[4]^2*summary(med.model)[[2]]$coefficients[2,2]^2 + 2*med.model1.coef*med.model2.coef*summary(full.model)$var[3,4] + 2*full.model$coefficients[3]*full.model$coefficients[4]*cov.mat.m[2,4]))
		cover.nde.scad <- as.numeric((nde.scad - qnorm(0.975) * nde.scad.se <= 1) & (1 <= nde.scad + qnorm(0.975) * nde.scad.se))
		cover.nie.prod.scad <-  as.numeric((nie.prod_scad - qnorm(0.975) * nie.prod.scad.se <= -0.16) & (-0.16 <= nie.prod_scad + qnorm(0.975) * nie.prod.scad.se))
	}else{
		cov.mat.m.sel <- cov.mat.m[2,2]
		nde.scad.se <- sqrt(scadcov[2,2])
		nie.prod.scad.se <- sqrt(med.model.coef^2 * scadselcov + beta_scad_sel^2 %*% med.model.se^2)  
		cover.nde.scad <- as.numeric((nde.scad - qnorm(0.975) * nde.scad.se <= 1) & (1 <= nde.scad[i] + qnorm(0.975) * nde.scad.se))
		cover.nie.prod.scad <-  as.numeric((nie.prod_scad - qnorm(0.975) * nie.prod.scad.se <= -0.16) & (-0.16 <= nie.prod_scad + qnorm(0.975) * nie.prod.scad.se))

	}
	#print(i)

# bootstrap

source("boot_pen.R")


cover.nde.lasso.b <- as.numeric((nde_lasso - qnorm(0.975) * sd_meff_lasso_b[1] <= 1) & (1 <= nde_lasso + qnorm(0.975) * sd_meff_lasso_b[1]))
cover.nie.prod.lasso.b <-  as.numeric((nie.prod_lasso - qnorm(0.975) * sd_meff_lasso_b[2] <= -0.16) & (-0.16 <= nie.prod_lasso + qnorm(0.975) * sd_meff_lasso_b[2]))

#cover.nde.alasso.b <- as.numeric((nde.alasso - qnorm(0.975) * sd_meff_alasso_b[1] <= 1) & (1 <= nde.alasso + qnorm(0.975) * sd_meff_alasso_b[1]))
#cover.nie.prod.alasso.b <-  as.numeric((nie.prod_alasso - qnorm(0.975) * sd_meff_alasso_b[2] <= -0.16) & (-0.16 <= nie.prod_alasso + qnorm(0.975) * sd_meff_alasso_b[2]))

#cover.nde.scad.b <- as.numeric((nde.scad - qnorm(0.975) * sd_meff_scad_b[1] <= 1) & (1 <= nde.scad + qnorm(0.975) * sd_meff_scad_b[1]))
#cover.nie.prod.scad.b <-  as.numeric((nie.prod_scad - qnorm(0.975) * sd_meff_scad_b[2] <= -0.16) & (-0.16 <= nie.prod_scad + qnorm(0.975) * sd_meff_scad_b[2]))

est <- cbind(nde_lasso,nie.prod_lasso,nie.diff_lasso,nde.alasso,nie.prod_alasso,nie.diff_alasso,nde.scad,nie.prod_scad,nie.diff_scad)
se <- cbind(nde.alasso.se,nie.prod.alasso.se,nde.scad.se,nie.prod.scad.se)
#se_b <- cbind(sd_meff_lasso_b[1],sd_meff_lasso_b[2],sd_meff_alasso_b[1],sd_meff_alasso_b[2],sd_meff_scad_b[1],sd_meff_scad_b[2])
se_b <- cbind(sd_meff_lasso_b[1],sd_meff_lasso_b[2])
cover <- cbind(cover.nde.alasso,cover.nie.prod.alasso,cover.nde.scad,cover.nie.prod.scad,cen.rate)
#cover_b <- cbind(cover.nde.lasso.b,cover.nie.prod.lasso.b,cover.nde.alasso.b,cover.nie.prod.alasso.b,cover.nde.scad.b,cover.nie.prod.scad.b)
cover_b <- cbind(cover.nde.lasso.b,cover.nie.prod.lasso.b)
sel <- t(as.matrix(c(as.numeric(sel_lasso_m),as.numeric(sel_alasso_m),as.numeric(sel_scad_m))))

all <- cbind(est,se,se_b,cover,cover_b)

res.m.all.seed <- paste('all_m.',seednum,'.csv',sep='')
res.m.sel.seed <- paste('sel_m.',seednum,'.csv',sep='')

write.csv(all,file=res.m.all.seed)
write.csv(sel,file=res.m.sel.seed)
