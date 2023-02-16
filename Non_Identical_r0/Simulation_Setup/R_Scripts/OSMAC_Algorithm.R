# Cordeiro Bias Estimation ----
Cordeiro<-function(XData,With_bias)
{
  p <- as.vector(exp(XData%*%as.vector(With_bias)))
  W <- diag(p)
  inverse_term <- solve(t(XData)%*%W%*%XData)
  
  Term1 <- inverse_term%*%t(XData)%*%W
  Term2 <- diag(diag(XData%*%(inverse_term)%*%t(XData))) %*% rep(-0.5,nrow(XData))
  
  bias <- as.vector(Term1%*%Term2)
  return(bias)
}

# Two step OSMAC ----
AlgTwoStp <- function(r1=r1, r2=r2,Y,X,Real_Data,n,Model) {
  if(Model=="Real_Model")
  {
    PI.prop <- rep(1/n, n)
    idx.prop <- sample(1:n, r1, T)
    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop]
    pinv.prop <- rep(n,r1)
    
    fit.prop<-glm(y.prop~x.prop-1,family = "poisson")
    beta.prop <-fit.prop$coefficients

    if (is.na(beta.prop[1]))
      return(list(opt=NA, msg="first stage not converge"))
    P.prop  <- exp(X%*% beta.prop)
    
    beta.mVc<-matrix(nrow = length(r2),ncol = ncol(X) )
    Utility_mVc<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mVc<-matrix(nrow = length(r2),ncol = ncol(X) )
    
    beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X) )
    Utility_mMSE<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mMSE<-matrix(nrow = length(r2),ncol = ncol(X) )
    
    #Sample.mMSE<-list()
    #Sample.mVc<-list()
    
    ## mVC
    for (i in 1:length(r2)) 
    {
      PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
      PI.mVc <- PI.mVc / sum(PI.mVc)
      
      idx.mVc <- sample(1:n, r2[i]-r1, T, PI.mVc)
      x.mVc <- X[c(idx.mVc, idx.prop),]
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
      fit.mVc <- glm(y.mVc~x.mVc-1,family="poisson",weights = pinv.mVc)
      #Sample.mVc[[i]]<-cbind(r2[i],y.mVc,1/pinv.mVc,x.mVc)
      
      ru <- length(y.mVc)
      beta.mVc[i,] <- fit.mVc$coefficients
      
      pi<- c(exp(x.mVc %*% beta.mVc[i,]))
      Mx<-solve(t(x.mVc) %*% (x.mVc*pi*pinv.mVc))
      V_Temp<-t(x.mVc) %*% (x.mVc*((as.vector(y.mVc)-pi)*pinv.mVc)^2)
      
      V_Final<-Mx %*% V_Temp %*% Mx
      
      Utility_mVc[i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final)))
      Bias_mVc[i,]<-Cordeiro(XData=x.mVc,With_bias = beta.mVc[i,])
      
      ## mMSE
      p.prop <- P.prop[idx.prop]
      w.prop <- p.prop
      W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))
      PI.mMSE <- sqrt((Y - P.prop)^2 * rowSums((X%*%W.prop)^2))
      PI.mMSE <- PI.mMSE / sum(PI.mMSE)
      idx.mMSE <- sample(1:n, r2[i]-r1, T, PI.mMSE)
      x.mMSE <- X[c(idx.mMSE, idx.prop),]
      y.mMSE <- Y[c(idx.mMSE, idx.prop)]
      pinv.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)
      fit.mMSE <- glm(y.mMSE~x.mMSE-1, weights=pinv.mMSE,family = "poisson")
      #Sample.mMSE[[i]]<-cbind(r2[i],y.mMSE,1/pinv.mMSE,x.mMSE)
      
      ru <- length(y.mMSE)
      beta.mMSE[i,] <- fit.mMSE$coefficients
      
      pi<-c(exp(x.mMSE %*% beta.mMSE[i,]))
      Mx<-solve(t(x.mMSE) %*% (x.mMSE*pi*pinv.mMSE))
      V_Temp<-t(x.mMSE) %*% (x.mMSE*((as.vector(y.mMSE)-pi)*pinv.mMSE)^2)
      
      V_Final<-Mx %*% V_Temp %*% Mx
      
      Utility_mMSE[i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final)))
      Bias_mMSE[i,]<-Cordeiro(XData=x.mMSE,With_bias = beta.mMSE[i,])
    }
    
    #Sample_mVc<-do.call(rbind,Sample.mVc)
    #Sample_mMSE<-do.call(rbind,Sample.mMSE)
    
    opt <- list("Est_Theta_mMSE"=cbind(r2,beta.mMSE),"Utility_mMSE"=Utility_mMSE,"Bias_mMSE"=cbind(r2,Bias_mMSE),
                #"Sample_mMSE"=Sample_mMSE,
                "Est_Theta_mVc"=cbind(r2,beta.mVc),"Utility_mVc"=Utility_mVc,"Bias_mVc"=cbind(r2,Bias_mVc)#,
                #"Sample_mVc"=Sample_mVc
                )
    msg <- c(fit.prop$message, fit.mMSE$message, fit.mVc$message)
    return(list(opt=opt, msg=msg)) 
  }
  else{
    x_Real<-Real_Data[,-1]
    
    PI.prop <- rep(1/n, n)
    idx.prop <- sample(1:n, r1, T)
    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop]
    x_Real.prop <- x_Real[idx.prop,] #
    pinv.prop <- rep(n,r1)
    
    fit.prop <- glm(y.prop~x.prop-1,family = "poisson")
    fit_Real.prop <- glm(y.prop~x_Real.prop-1,family="poisson") #
    beta.prop <- fit.prop$coefficients
    beta_Real.prop <- fit_Real.prop$coefficients #
    
    if (is.na(beta.prop[1]) & is.na(beta_Real.prop[1]))
      return(list(opt=NA, msg="first stage not converge"))
    P.prop  <- exp(X %*% beta.prop)
    P_Real.prop  <- exp(x_Real %*% beta_Real.prop) #
    
    beta.mVc<-matrix(nrow = length(r2),ncol = ncol(x_Real) )
    Utility_mVc<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mVc<-matrix(nrow = length(r2),ncol = ncol(x_Real) )
    
    beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(x_Real) )
    Utility_mMSE<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mMSE<-matrix(nrow = length(r2),ncol = ncol(x_Real) )
    
    #Sample.mMSE<-list()
    #Sample.mVc<-list()
    
    ## mVC
    for (i in 1:length(r2)) 
    {
      PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
      PI.mVc <- PI.mVc / sum(PI.mVc)
      
      idx.mVc <- sample(1:n, r2[i]-r1, T, PI.mVc)
      
      x_Real.mVc <- x_Real[c(idx.mVc, idx.prop),]
      y_Real.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv_Real.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
      fit.Real.mVc<-glm(y_Real.mVc~x_Real.mVc-1,family = "poisson",weights = pinv_Real.mVc)
      #Sample.mVc[[i]]<-cbind(r2[i],y_Real.mVc,1/pinv_Real.mVc,x_Real.mVc)
      
      ru <- length(y_Real.mVc)
      beta.mVc[i,] <- fit_Real.mVc$coefficients
      
      pi<-c(exp(x_Real.mVc %*% beta.mVc[i,]))
      Mx<-solve(t(x_Real.mVc) %*% (x_Real.mVc*pi*pinv_Real.mVc))
      V_Temp<-t(x_Real.mVc) %*% (x_Real.mVc*((as.vector(y_Real.mVc)-pi)*pinv_Real.mVc)^2)
      
      V_Final<-Mx %*% V_Temp %*% Mx
      
      Utility_mVc[i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final)))
      Bias_mVc[i,]<-Cordeiro(XData=x_Real.mVc,With_bias = beta.mVc[i,])
      
      ## mMSE
      p.prop <- P.prop[idx.prop]
      w.prop <- p.prop
      W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))
      PI.mMSE <- sqrt((Y - P.prop)^2 * rowSums((X%*%W.prop)^2))
      PI.mMSE <- PI.mMSE / sum(PI.mMSE)
      
      idx.mMSE <- sample(1:n, r2[i]-r1, T, PI.mMSE)
      x_Real.mMSE <- x_Real[c(idx.mMSE, idx.prop),]
      y_Real.mMSE <- Y[c(idx.mMSE, idx.prop)]
      pinv_Real.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)
      fit_Real.mMSE <-glm(y_Real.mMSE~x_Real.mMSE-1,weights = pinv_Real.mMSE,family = "poisson")
      #Sample.mMSE[[i]]<-cbind(r2[i],y_Real.mMSE,1/pinv_Real.mMSE,x_Real.mMSE)
      
      ru <- length(y_Real.mMSE)
      beta.mMSE[i,] <- fit_Real.mMSE$coefficients
      
      pi<-c(exp(x_Real.mMSE %*% beta.mMSE[i,]))
      Mx<-solve(t(x_Real.mMSE) %*% (x_Real.mMSE*pi*pinv_Real.mMSE))
      V_Temp<-t(x_Real.mMSE) %*% (x_Real.mMSE*((as.vector(y_Real.mMSE)-pi)*pinv_Real.mMSE)^2)
      
      V_Final<-Mx %*% V_Temp %*% Mx
      
      Utility_mMSE[i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final)))
      Bias_mMSE[i,]<-Cordeiro(XData=x_Real.mMSE,With_bias = beta.mMSE[i,])
    }
    
    if(anyNA(beta.mMSE) || anyNA(beta.mVc) || anyNA(Utility_mMSE) || anyNA(Utility_mVc) || anyNA(Bias_mMSE) || anyNA(Bias_mVc) )
    {
      stop("There are NA or NaN values")
    }
    
    #Sample_mVc<-do.call(rbind,Sample.mVc)
    #Sample_mMSE<-do.call(rbind,Sample.mMSE)
    
    opt <- list("Est_Theta_mMSE"=cbind(r2,beta.mMSE),"Utility_mMSE"=Utility_mMSE,"Bias_mMSE"=cbind(r2,Bias_mMSE),
                #"Sample_mMSE"=Sample_mMSE,
                "Est_Theta_mVc"=cbind(r2,beta.mVc),"Utility_mVc"=Utility_mVc,"Bias_mVc"=cbind(r2,Bias_mVc)#,
                #"Sample_mVc"=Sample_mVc
                )
    msg <- c(fit_Real.prop$message, fit_Real.mMSE$message, fit_Real.mVc$message)
    return(list(opt=opt, msg=msg)) 
  }
}
