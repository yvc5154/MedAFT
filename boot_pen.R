# LASSO with bootstrapped sample 
    
boot.pen.lassos <- function(fit.all.b,obs.b,status.b,Mmat.b,A.b,x.b,type){    	      
Mno.b <- Mmat.b[,3:(p-1)]
   sel_lasso_m.b <- rep(0, mm)
   sel_alasso_m.b <- rep(0, mm)

      # Newton-Rapshon 
   if(type == "lasso"){      
      betahat<-c(fit.all.b$coefficients)
      sigma2<-(fit.all.b$scale^2)
      iter<-1
      
      repeat{
        
        	eta<-Mmat.b%*%betahat
        	theta<-as.matrix(c(betahat,sigma2))
        	m<-(log(obs.b)-eta)/(sigma2^(1/2))
        	n_p<-dnorm(m,mean=0,sd=1)
        	n_c<-(1-pnorm(m,mean=0,sd=1))
        	V <- ifelse(n_c == 0, 0, n_p/n_c)  
        	
        	
        
        	WL=diag(c(0,0,c(tun_lasso/abs(betahat[3:(p-1)]+10^(-6))),0))
        	d_b=(1/(sigma2^(1/2)))*(t(Mmat.b)%*%(status.b*m+(1-status.b)*V))-n*WL%*%betahat
        	d_p=(1/(2*(sigma2)))*sum(status.b*(m^2-1)+(1-status.b)*V*m)
        	d1_v=rbind(d_b,d_p)
        
        	xi=V*(V-m)
        	W=diag(c((1/sigma2)*(status.b+(1-status.b)*xi)))
        
        	d2_b=t(Mmat.b)%*%W%*%Mmat.b+n*WL
        	d2_p=(1/(2*(sigma2^2)))*sum(status.b*(2*(m^2)-1)+(1-status.b)*m*((m*xi+3*V)/2))
        	d2_bp=(1/(2*(sigma2^(3/2))))*(t(Mmat.b)%*%(status.b*2*m+(1-status.b)*(m*xi+V)))
        
        	H=rbind(cbind(d2_b,d2_bp),cbind(t(d2_bp),d2_p))
        
        	#inv_H=solve(H)
        	if(!is.singular.matrix(H)){
        	  inv_H=solve(H)
        	}else{
        	  inv_H=ginv(H)
        	}
        	theta_a=as.numeric(theta+inv_H%*%d1_v)
        
        	beta_h<-theta_a[1:p]
        	sigma2_h<-theta_a[(p+1)]
        
        	if((max(abs(theta_a-theta))<10^(-6))|(iter>500)|(sigma2_h <= 0.001))(break)
        
        	iter<-iter+1
        	betahat<-beta_h
        	sigma2<-sigma2_h
        
      }
      val_lassos <- c(betahat, sigma2)
      sel_lasso.b <-  which(abs(betahat[3:(p-1)]) > 1e-3)
      Mmat_lasso.b <- Mno.b[,sel_lasso.b]
      sel_lasso_m.b[sel_lasso.b]<- ifelse(sel_lasso.b >=1, 1, 0)

      # Run mediator model
      med.model_lasso.b <- lm(Mmat_lasso.b ~ A.b + x.b)

       nde_lasso.b <- betahat[2]
    
    sel_est_lasso.b <- sel_lasso.b + 2
    
    if(length(sel_lasso.b) > 1){
      nie.prod_lasso.b <- med.model_lasso.b$coefficients[2,]%*%betahat[sel_est_lasso.b]
    }else{
      nie.prod_lasso.b <- med.model_lasso.b$coefficients[2] * betahat[sel_est_lasso.b]
    }

   val_lassos <- c(betahat,sigma2,nde_lasso.b,nie.prod_lasso.b)

  }else{
       beta_ini <- c(fit.all.b$coefficients[3:(p-1)])
	weight0=abs(1/beta_ini)
       # Newton-Rapshon 
      
      	betahat<-c(fit.all.b$coefficients)
      	sigma2<-(fit.all.b$scale^2)
      	iter<-1
      
      repeat{
        
        	eta<-Mmat.b%*%betahat
        	theta<-as.matrix(c(betahat,sigma2))
        	m<-(log(obs.b)-eta)/(sigma2^(1/2))
        	n_p<-dnorm(m,mean=0,sd=1)
        	n_c<-(1-pnorm(m,mean=0,sd=1))
        	V<-n_p/(n_c)
        
        	# 
        	WL=diag(c(0,0,c(tun_alasso/abs(betahat[3:(p-1)]+10^(-6)))*weight0,0))
        	d_b=(1/(sigma2^(1/2)))*(t(Mmat.b)%*%(status.b*m+(1-status.b)*V))-n*WL%*%betahat
        	d_p=(1/(2*(sigma2)))*sum(status.b*(m^2-1)+(1-status.b)*V*m)
        	d1_v=rbind(d_b,d_p)
        
        	xi=V*(V-m)
        	W=diag(c((1/sigma2)*(status.b+(1-status.b)*xi)))
        
        	d2_b=t(Mmat.b)%*%W%*%Mmat.b+n*WL
        	d2_p=(1/(2*(sigma2^2)))*sum(status.b*(2*(m^2)-1)+(1-status.b)*m*((m*xi+3*V)/2))
        	d2_bp=(1/(2*(sigma2^(3/2))))*(t(Mmat.b)%*%(status.b*2*m+(1-status.b)*(m*xi+V)))
        
        	H=rbind(cbind(d2_b,d2_bp),cbind(t(d2_bp),d2_p))
        
        	inv_H=solve(H)
        	theta_a=as.numeric(theta+inv_H%*%d1_v)
        
        	betahat_h<-theta_a[1:p]
        	sigma2_h<-theta_a[(p+1)]
        
        	if((max(abs(theta_a-theta))<10^(-6))|(iter>500))(break)
        
        	iter<-iter+1
        	betahat<-betahat_h
        	sigma2<-sigma2_h
        
      	}
   sel_alasso.b <-  which(abs(betahat[3:(p-1)]) > 1e-3)
    
    Mmat_alasso.b <- Mno.b[,sel_alasso.b]
    
    sel_alasso_m.b[sel_alasso.b]<- ifelse(sel_alasso.b >=1, 1, 0)
    med.model_alasso.b <- lm(Mmat_alasso.b ~ A.b + x.b)

    nde_alasso.b <- betahat[2]
    
    sel_est_alasso.b <- sel_alasso.b + 2
    if(length(sel_alasso.b) > 1){
      	nie.prod_alasso.b <- med.model_alasso.b$coefficients[2,]%*%betahat[sel_est_alasso.b]
    	}else{
      	nie.prod_alasso.b <- med.model_alasso.b$coefficients[2] * betahat[sel_est_alasso.b]
    	}
   val_lassos <- c(betahat,sigma2,nde_alasso.b,nie.prod_alasso.b)

  }
 return(val_lassos)
}
   
# SCAD 
boot.pen.scad <- function(fit.all.b,obs.b,status.b,Mmat.b,A.b,x.b,beta_lasso){
Mno.b <- Mmat.b[,3:(p-1)]
sel_scad_m.b <- rep(0, mm)

 # Newton-Rapshon 
      
      	betahat<-c(beta_lasso)
      	sigma2<-(fit.all.b$scale^2)
      	iter<-1
            tun<-c(rep(tun_scad,p-3))
      repeat{
        
        	eta<-Mmat.b%*%betahat
        	theta<-as.matrix(c(betahat,sigma2))
        	m<-(log(obs.b)-eta)/(sigma2^(1/2))
        	n_p<-dnorm(m,mean=0,sd=1)
        	n_c<-(1-pnorm(m,mean=0,sd=1))
        	V<-n_p/(n_c)
        
        	# penalty SCAD
        	a=3.7
        	nft<-vector()
        	for ( k in 3:(p-1)){
          		if(0<abs(betahat[k]) & abs(betahat[k])<=tun[k-2]){
            		nft[k-2]=tun[k-2]}else if(tun[k-2]<abs(betahat[k]) & abs(betahat[k])<a*tun[k-2]){
              		nft[k-2]=(a*tun[k-2]-abs(betahat[k]))/(a-1)} else {
                	nft[k-2]=0
            	}
        	}
        
        	#WL=diag(c((nft/abs(betahat+10^(-6)))))
        	WL=diag(c(0,0,c((nft/abs(betahat[3:(p-1)]+10^(-6)))),0))
        	d_b=(1/(sigma2^(1/2)))*(t(Mmat.b)%*%(status.b*m+(1-status.b)*V))-n*WL%*%betahat
        	d_p=(1/(2*(sigma2)))*sum(status.b*(m^2-1)+(1-status.b)*V*m)
        	d1_v=rbind(d_b,d_p)
        
        	xi=V*(V-m)
        	W=diag(c((1/sigma2)*(status.b+(1-status.b)*xi)))
        
        	d2_b=t(Mmat.b)%*%W%*%Mmat.b+n*WL
        	d2_p=(1/(2*(sigma2^2)))*sum(status.b*(2*(m^2)-1)+(1-status.b)*m*((m*xi+3*V)/2))
        	d2_bp=(1/(2*(sigma2^(3/2))))*(t(Mmat.b)%*%(status.b*2*m+(1-status.b)*(m*xi+V)))
        
        	H=rbind(cbind(d2_b,d2_bp),cbind(t(d2_bp),d2_p))
        
        	inv_H=solve(H)
        	theta_a=as.numeric(theta+inv_H%*%d1_v)
        
        	betahat_h<-theta_a[1:p]
        	sigma2_h<-theta_a[(p+1)]
        
        	if((max(abs(theta_a-theta))<10^(-6))|(iter>500)|(sigma2_h <= 0.001))(break)
        
        	iter<-iter+1
        	betahat<-betahat_h
        	sigma2<-sigma2_h
        
      }
       sel_scad.b <-  which(abs(betahat[3:(p-1)]) > 1e-3)
    
    Mmat_scad.b <- Mno[,sel_scad.b]
#    sel_scad_m<- abs(beta_scad[3:(p-1)]) > 1e-3
    

    sel_scad_m.b[sel_scad.b]<- ifelse(sel_scad.b >=1, 1, 0)
    med.model_scad.b <- lm(Mmat_scad.b ~ A.b + x.b)
    
    nde.scad.b <- betahat[2]
    
    sel_est_scad.b <- sel_scad.b + 2
    
    if(length(sel_scad.b) > 1){
      	nie.prod_scad.b <- med.model_scad.b$coefficients[2,]%*%betahat[sel_est_scad.b]
    	}else{
      	nie.prod_scad.b <- med.model_scad.b$coefficients[2] * betahat[sel_est_scad.b]
    	}
        val_scad <- c(betahat,sigma2,nde.scad.b, nie.prod_scad.b)

  return(val_scad)               
}

B <- 200
meff_lasso_b <- matrix(0,ncol=2,nrow=B)
meff_alasso_b <- matrix(0,ncol=2,nrow=B)
meff_scad_b <- matrix(0,ncol=2,nrow=B)
for(b in 1:B){
ind <- sample(1:n, size=n, replace=T)
obs.b <- obs[ind]
status.b <- status[ind]
Mmat.b <- Mmat[ind,]
A.b <- A[ind]
x.b <- x[ind]

fit.all.b <- survreg(Surv(obs.b,status.b) ~ -1 + Mmat.b, dist = "lognormal")

work.mis.hs <- try(boot.pen.lassos(fit.all.b,obs.b,status.b,Mmat.b,A.b,x.b,"lasso"))
if(class(work.mis.hs) == "try-error"){
  break
}else{
  theta_lasso_b <- work.mis.hs
}

#theta_lasso_b <- boot.pen.lassos(fit.all.b,obs.b,status.b,Mmat.b,A.b,x.b,"lasso")
#theta_alasso_b <- boot.pen.lassos(fit.all.b,obs.b,status.b,Mmat.b,A.b,x.b,"alasso")
beta_lasso_b <- theta_lasso_b[(1:p)]
#theta_scad_b <- boot.pen.scad(fit.all.b,obs.b,status.b,Mmat.b,A.b,x.b,beta_lasso_b)
meff_lasso_b[b,] <- theta_lasso_b[((p+2):(p+3))] 
#meff_alasso_b[b,] <- theta_alasso_b[((p+2):(p+3))]
#meff_scad_b[b,] <- theta_scad_b[((p+2):(p+3))] 
}

if(b == B){
  sd_meff_lasso_b <- apply(meff_lasso_b,2,sd)
}else{
  sd_meff_lasso_b <- apply(meff_lasso_b[1:(b-1),],2,sd)
}
#sd_meff_alasso_b <- apply(meff_alasso_b,2,sd)
#sd_meff_scad_b <- apply(meff_scad_b,2,sd)