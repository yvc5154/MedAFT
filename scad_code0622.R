 ### Start SCAD
    	# Range of tuning parameter 
    	start<-0
    	end<-0.3
    	by<-0.01
    	tun_range<-seq(start,end,by)
    	lth<-length(c(tun_range))
    
    	count<-1
    
    	iter_v<-vector()
    	bic_v<-vector() # BIC
    	tun_v<-vector() # Tuning parameter
    	beta_m<-matrix(nrow=lth,ncol=p)      # Est_beta_matrix
    	#beta_d1<-matrix(nrow=lth,ncol=p)
    	sigma2_v<-vector()
    	#sigma2_d1<-vector()
    	SE_beta_m<-matrix(nrow=lth,ncol=p)   # SE_beta_matrix
    	SE_sigma2_v<-vector()  # SE_sigma2_vector
      cov.arr <- array(0,c(p,p,lth))
    	for (j in tun_range){
      
      	tun<-c(rep(j,p-3))
      
      # Newton-Rapshon 
      
      	betahat<-c(beta_lasso)
      	sigma2<-(fit.all$scale^2)
      	iter<-1
      
      	repeat{
        
        	eta<-Mmat%*%betahat
        	theta<-as.matrix(c(betahat,sigma2))
        	m<-(log(obs)-eta)/(sigma2^(1/2))
        	n_p<-dnorm(m,mean=0,sd=1)
        	n_c<-(1-pnorm(m,mean=0,sd=1))
        	#V<-n_p/(n_c)
        	V <- ifelse(n_c == 0, 0, n_p/n_c)
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
        	d_b=(1/(sigma2^(1/2)))*(t(Mmat)%*%(status*m+(1-status)*V))-n*WL%*%betahat
        	d_p=(1/(2*(sigma2)))*sum(status*(m^2-1)+(1-status)*V*m)
        	d1_v=rbind(d_b,d_p)
        
        	xi=V*(V-m)
        	W=diag(c((1/sigma2)*(status+(1-status)*xi)))
        
        	d2_b=t(Mmat)%*%W%*%Mmat+n*WL
        	d2_p=(1/(2*(sigma2^2)))*sum(status*(2*(m^2)-1)+(1-status)*m*((m*xi+3*V)/2))
        	d2_bp=(1/(2*(sigma2^(3/2))))*(t(Mmat)%*%(status*2*m+(1-status)*(m*xi+V)))
        
        	H=rbind(cbind(d2_b,d2_bp),cbind(t(d2_bp),d2_p))
        
        	#inv_H=solve(H)
        	if(!is.singular.matrix(H)){
        	  inv_H=solve(H)
        	}else{
        	  inv_H=ginv(H)
        	}
        	theta_a=as.numeric(theta+inv_H%*%d1_v)
        
        	betahat_h<-theta_a[1:p]
        	sigma2_h<-theta_a[(p+1)]
        
        	if((max(abs(theta_a-theta))<10^(-6))|(iter>500))(break)
        
        	iter<-iter+1
        	betahat<-betahat_h
        	sigma2<-sigma2_h
        
      	}
      
      	# Log likelihood
      	log_L<-sum(status*(-(1/2)*log(sigma2)+log(n_p))+(1-status)*log(n_c))
      
      	# BIC
      	H_bb<-t(Mmat)%*%W%*%Mmat
      	#df<-sum(diag(solve(H_bb+n*WL)%*%H_bb))
      	#if(!is.singular.matrix(H_bb+n*WL)){
      	 # df<-sum(diag(solve(H_bb+n*WL)%*%H_bb))
      	#}else{
      	 # df<-sum(diag(ginv(H_bb+n*WL)%*%H_bb))
      	#}
      	
      	if(!is.singular.matrix(H_bb+n*WL)){
      	  df<-sum(diag(solve(H_bb+n*WL)%*%H_bb))
      	  # SE
      	  cov<-solve(H_bb+n*WL)%*%H_bb%*%solve(H_bb+n*WL)
      	  SE_beta_hat<-sqrt(diag(cov))                # beta_hat
      	  SE_sigma2_hat<-sqrt(diag(solve(d2_p)))         # sigma2_hat
      	  
      	}else{
      	  df<-sum(diag(ginv(H_bb+n*WL)%*%H_bb))
      	  # SE
      	  cov<-ginv(H_bb+n*WL)%*%H_bb%*%ginv(H_bb+n*WL)
      	  SE_beta_hat<-sqrt(diag(cov))                # beta_hat
      	  SE_sigma2_hat<-sqrt(diag(ginv(d2_p)))         # sigma2_hat
      	}
      	bic<--2*log_L+log(n)*df
      
      	# SE
      	#cov<-solve(H_bb+n*WL)%*%H_bb%*%solve(H_bb+n*WL)
      
      	#SE_beta_hat<-sqrt(diag(cov))                # beta_hat
      	#SE_sigma2_hat<-sqrt(diag(solve(d2_p)))         # sigma2_hat
      
      	# Tuning parameter
      	iter_v[count]<-iter
      	tun_v[count]<-tun[1]
      	bic_v[count]<-bic
      	beta_m[count,]<-betahat         # Est_beta_matrix
      	#beta_d1[count,]<-d_b
      	sigma2_v[count]<-sigma2
      	#sigma2_d1[count]<-d_p
      	SE_beta_m[count,]<-SE_beta_hat    # SE_beta_matrix
      	SE_sigma2_v[count]<-SE_sigma2_hat
		cov.arr[,,count] <- cov
      	count<-count+1
            
    		}