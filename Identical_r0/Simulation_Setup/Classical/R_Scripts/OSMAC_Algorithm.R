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
AlgTwoStp <- function(r1=r1, r2=r2,Y,X,n,alpha,All_Combinations,All_Covariates){
    PI.prop <- rep(1/n, n)
    idx.prop <- sample(1:n, r1, T)
    
    x.prop<-lapply(1:length(All_Combinations),function(j){
      X[idx.prop,All_Covariates %in% All_Combinations[[j]]]
    })
    y.prop <- Y[idx.prop,]
    
    pinv.prop <- n
    pinv.prop <- 1/PI.prop[idx.prop]
    
    fit.prop <- lapply(1:length(All_Combinations), function(j){
      glm(y.prop~x.prop[[j]]-1,family = "poisson")
    })
    
    beta.prop<-list()
    for (j in 1:length(All_Combinations)) 
      {
      beta.prop[[j]] <- fit.prop[[j]]$coefficients
      if(anyNA(beta.prop[[j]]))
        {
        return(list(opt=NA, msg="first stage not converge"))
        }
      }
    
    P.prop  <- lapply(1:length(All_Combinations),function(j){
      exp(X[,All_Covariates %in% All_Combinations[[j]] ] %*% beta.prop[[j]])
    })
    
    beta.mVc_Single<-Utility_mVc_Single<-Bias_mVc_Single<-list() # Single Models
    beta.mMSE_Single<-Utility_mMSE_Single<-Bias_mMSE_Single<-list()
    #Sample.mMSE_Single<-list()
    #Sample.mVc_Single<-list()
    beta.mVc_MF<-Utility_mVc_MF<-Bias_mVc_MF<-list() # Model Free Models
    beta.mMSE_MF<-Utility_mMSE_MF<-Bias_mMSE_MF<-list()
    #Sample.mMSE_MF<-list()
    #Sample.mVc_MF<-list()
    
    # For the models, single and Model Free
    for (a in 1:length(All_Combinations)) 
    {
      beta.mVc_Single[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]]) ) # Single models
      Utility_mVc_Single[[a]]<-matrix(nrow = length(r2),ncol = 3 )
      Bias_mVc_Single[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]]) )

      beta.mMSE_Single[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]]) )
      Utility_mMSE_Single[[a]]<-matrix(nrow = length(r2),ncol = 3 )
      Bias_mMSE_Single[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]]) )
      
      #Sample.mMSE_Single[[a]]<-list()
      #Sample.mVc_Single[[a]]<-list()
      
      beta.mVc_MF[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]]) ) # Model Free models
      Utility_mVc_MF[[a]]<-matrix(nrow = length(r2),ncol = 3 )
      Bias_mVc_MF[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]]) )
    
      beta.mMSE_MF[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]]) )
      Utility_mMSE_MF[[a]]<-matrix(nrow = length(r2),ncol = 3 )
      Bias_mMSE_MF[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]]) )
    
      #Sample.mMSE_MF[[a]]<-list()
      #Sample.mVc_MF[[a]]<-list()
    }
    
    for (i in 1:length(r2)) 
    {
      ## mVC
      PI_Single.mVc <- lapply(1:length(All_Combinations), function(j){
        PI.mVc<-sqrt((Y - P.prop[[j]])^2 * rowSums(X[,All_Covariates %in% All_Combinations[[j]] ]^2)) # Single Model 
        return(PI.mVc/sum(PI.mVc))
      })
      PI_MF.mVc<-rowSums2(do.call(cbind,PI_Single.mVc)%*%diag(alpha)) # Model Free Model
        
      idx_Single.mVc <- lapply(1:length(All_Combinations), function(j){
        sample(1:n, r2[i]-r1, T, PI_Single.mVc[[j]]) # Single Model
      })
      idx_MF.mVc <- sample(1:n, r2[i]-r1, T, PI_MF.mVc) # Model Free Model
      
      x_Single.mVc <-lapply(1:length(All_Combinations),function(j){
        X[c(idx_Single.mVc[[j]], idx.prop),All_Covariates %in% All_Combinations[[j]] ] # Single Model
      }) 
      y_Single.mVc <- lapply(1:length(All_Combinations),function(j){
        Y[c(idx_Single.mVc[[j]], idx.prop)]  # Single Model
      })
        
      x_MF.mVc <-  lapply(1:length(All_Combinations),function(j){
        X[c(idx_MF.mVc, idx.prop),All_Covariates %in% All_Combinations[[j]] ] # Model Free Model
      })
      y_MF.mVc <- Y[c(idx_MF.mVc, idx.prop)] 
      pinv_MF.mVc<-c(1 / PI_MF.mVc[idx_MF.mVc], pinv.prop)
        
      fit_Single.mVc <-lapply(1:length(All_Combinations), function(j){
        pinv_Single.mVc<-list()
        pinv_Single.mVc[[j]]<-c(1 / PI_Single.mVc[[j]][idx_Single.mVc[[j]]], pinv.prop)
        glm(y_Single.mVc[[j]]~x_Single.mVc[[j]]-1, family = "poisson",
            weights=pinv_Single.mVc[[j]]) # Single Model
      })
      
      fit_MF.mVc <- lapply(1:length(All_Combinations),function(j){
        glm(y_MF.mVc~x_MF.mVc[[j]]-1,family="poisson",
            weights=pinv_MF.mVc) # Model Free Model
      })
        
      # for (j in 1:length(All_Combinations))
      # {
      #   Sample.mVc_Single[[j]][[i]]<-cbind(r2[i],y_Single.mVc[[j]],x_Single.mVc[[j]],
      #                                      c(PI_Single.mVc[[j]][idx_Single.mVc[[j]] ], 1 / pinv.prop))
      #   Sample.mVc_MF[[j]][[i]]<-cbind(r2[i],y_MF.mVc,x_MF.mVc[[j]],
      #                                  c(PI_MF.mVc[idx_MF.mVc], 1 / pinv.prop))
      # }

      if(anyNA(fit_Single.mVc$coefficients) || anyNA(fit_MF.mVc$coefficients) )
      {
        stop("There are NA or NaN values")
      }
      
      for (j in 1:length(All_Combinations)) 
      {
        beta.mVc_Single[[j]][i,] <- fit_Single.mVc[[j]]$coefficients 
        beta.mVc_MF[[j]][i,] <- fit_MF.mVc[[j]]$coefficients
      }
      
      # Single Model
      pi<- lapply(1:length(All_Combinations), function(j){
        c(exp(x_Single.mVc[[j]] %*% beta.mVc_Single[[j]][i,]))
      })
      Mx<-lapply(1:length(All_Combinations),function(j){
        pinv_Single.mVc<-list()
        pinv_Single.mVc[[j]]<-c(1 / PI_Single.mVc[[j]][idx_Single.mVc[[j]]], pinv.prop)
        solve(t(x_Single.mVc[[j]]) %*% (x_Single.mVc[[j]] * pi[[j]] * pinv_Single.mVc[[j]]))
      }) 
      V_Temp<-lapply(1:length(All_Combinations), function(j){
        pinv_Single.mVc<-list()
        pinv_Single.mVc[[j]]<-c(1 / PI_Single.mVc[[j]][idx_Single.mVc[[j]]], pinv.prop)
        t(x_Single.mVc[[j]]) %*% (x_Single.mVc[[j]]*((as.vector(y_Single.mVc[[j]])-pi[[j]])*pinv_Single.mVc[[j]])^2) 
      })
      V_Final<-lapply(1:length(All_Combinations),function(j){
        Mx[[j]] %*% V_Temp[[j]] %*% Mx[[j]]
      }) 
      
      for (j in 1:length(All_Combinations)) 
      {
        Utility_mVc_Single[[j]][i,]<-cbind(r2[i],tr(V_Final[[j]]),det(solve(V_Final[[j]])))
        Bias_mVc_Single[[j]][i,]<-Cordeiro(XData=x_Single.mVc[[j]],With_bias = beta.mVc_Single[[j]][i,])    
      }
      
      # Model Free
      pi<- lapply(1:length(All_Combinations), function(j){
        c(exp(x_MF.mVc[[j]] %*% beta.mVc_MF[[j]][i,]))
      })
      Mx<-lapply(1:length(All_Combinations),function(j){
        solve(t(x_MF.mVc[[j]]) %*% (x_MF.mVc[[j]]*pi[[j]]*pinv_MF.mVc))
      })
      V_Temp<-lapply(1:length(All_Combinations), function(j){
        t(x_MF.mVc[[j]]) %*% (x_MF.mVc[[j]]*((as.vector(y_MF.mVc)-pi[[j]])*pinv_MF.mVc)^2)
      })
      V_Final<-lapply(1:length(All_Combinations),function(j){
        Mx[[j]] %*% V_Temp[[j]] %*% Mx[[j]]
      }) 
      
      for (j in 1:length(All_Combinations)) 
      {
        Utility_mVc_MF[[j]][i,]<-cbind(r2[i],tr(V_Final[[j]]),det(solve(V_Final[[j]])))
        Bias_mVc_MF[[j]][i,]<-Cordeiro(XData=x_MF.mVc[[j]],With_bias = beta.mVc_MF[[j]][i,])    
      }
      
      ## mMSE
      p_Single.prop <- lapply(1:length(All_Combinations),function(j){
        P.prop[[j]][idx.prop] # Single Model
      })
      w_Single.prop <- p_Single.prop
      W_Single.prop <- lapply(1:length(All_Combinations),function(j){
        solve(t(x.prop[[j]]) %*% (x.prop[[j]] * w_Single.prop[[j]] * pinv.prop)) # Single Model
      })
      PI_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
        Pi_mMSE<-sqrt((Y - P.prop[[j]])^2 * rowSums((X[,All_Covariates %in% All_Combinations[[j]] ]%*%
                                                     W_Single.prop[[j]])^2)) # Single Model
        return(Pi_mMSE/sum(Pi_mMSE))
      })
      PI_MF.mMSE<-rowSums2(do.call(cbind,PI_Single.mMSE)%*%diag(alpha))  # Model Free
      
      idx_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
        sample(1:n, r2[i]-r1, T, PI_Single.mMSE[[j]]) # Single Model
      }) 
      idx_MF.mMSE <- sample(1:n, r2[i]-r1, T, PI_MF.mMSE) # Model Free
      
      x_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
        X[c(idx_Single.mMSE[[j]], idx.prop),All_Covariates %in% All_Combinations[[j]] ] # Single Model
      })
      y_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
        Y[c(idx_Single.mMSE[[j]], idx.prop)] # Single Model
      })
      
      x_MF.mMSE <- lapply(1:length(All_Combinations),function(j){
        X[c(idx_MF.mMSE, idx.prop),All_Covariates %in% All_Combinations[[j]] ] # Model Free
      })  
      y_MF.mMSE <- Y[c(idx_MF.mMSE, idx.prop)] # Model Free
      pinv_MF.mMSE<-c(1 / PI_MF.mMSE[idx_MF.mMSE], pinv.prop)
      
      fit_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
        pinv_Single.mMSE<-list()
        pinv_Single.mMSE[[j]]<-c(1 / PI_Single.mMSE[[j]][idx_Single.mMSE[[j]]], pinv.prop)
        glm(y_Single.mMSE[[j]]~x_Single.mMSE[[j]]-1, family = "poisson", 
            weights=pinv_Single.mMSE[[j]]) # Single Model
      })
        
      fit_MF.mMSE <- lapply(1:length(All_Combinations), function(j){
        glm(y_MF.mMSE~x_MF.mMSE[[j]]-1, family = "poisson",
            weights=pinv_MF.mMSE) # Model Free
      })
        
      # for (j in 1:length(All_Combinations))
      # {
      #   Sample.mMSE_Single[[j]][[i]]<-cbind(r2[i],y_Single.mMSE[[j]],x_Single.mMSE[[j]],
      #                                            c(PI_Single.mMSE[[j]][idx_Single.mMSE[[j]]], 1 / pinv.prop)) # Single Model
      #   
      #   Sample.mMSE_MF[[i]]<-cbind(r2[i],y_MF.mMSE,x_MF.mMSE[[j]],
      #                              c(PI_MF.mMSE[idx_MF.mMSE], 1 / pinv.prop)) # Model Free
      # }
      
      if(anyNA(fit_Single.mMSE$coefficients) || anyNA(fit_MF.mMSE$coefficients))
      {
        stop("There are NA or NaN values")
      }
      
      for (j in 1:length(All_Combinations)) 
      {
        beta.mMSE_Single[[j]][i,] <-fit_Single.mMSE[[j]]$coefficients # Single Model
        beta.mMSE_MF[[j]][i,] <-fit_MF.mMSE[[j]]$coefficients # Model Free
      }
      
      # Single Model
      pi<-lapply(1:length(All_Combinations),function(j){
        c(exp(x_Single.mMSE[[j]] %*% beta.mMSE_Single[[j]][i,]))
      }) 
      Mx<-lapply(1:length(All_Combinations),function(j){
        pinv_Single.mMSE<-list()
        pinv_Single.mMSE[[j]]<-c(1 / PI_Single.mMSE[[j]][idx_Single.mMSE[[j]]], pinv.prop)
        solve(t(x_Single.mMSE[[j]]) %*% (x_Single.mMSE[[j]]*pi[[j]]*pinv_Single.mMSE[[j]]))
      })
      V_Temp<-lapply(1:length(All_Combinations),function(j){
        pinv_Single.mMSE<-list()
        pinv_Single.mMSE[[j]]<-c(1 / PI_Single.mMSE[[j]][idx_Single.mMSE[[j]]], pinv.prop)
        t(x_Single.mMSE[[j]]) %*% (x_Single.mMSE[[j]]*((as.vector(y_Single.mMSE[[j]])-pi[[j]])*pinv_Single.mMSE[[j]])^2)
        }) 
      V_Final<-lapply(1:length(All_Combinations),function(j){
        Mx[[j]] %*% V_Temp[[j]] %*% Mx[[j]]
      }) 
      
      for (j in 1:length(All_Combinations)) 
      {
        Utility_mMSE_Single[[j]][i,]<-cbind(r2[i],tr(V_Final[[j]]),det(solve(V_Final[[j]])))
        Bias_mMSE_Single[[j]][i,]<-Cordeiro(XData=x_Single.mMSE[[j]],With_bias = beta.mMSE_Single[[j]][i,])    
      }
      
      # Model Free
      pi<-lapply(1:length(All_Combinations),function(j){
        c(exp(x_MF.mMSE[[j]] %*% beta.mMSE_MF[[j]][i,]))
      })
      Mx<-lapply(1:length(All_Combinations),function(j){
        solve(t(x_MF.mMSE[[j]]) %*% (x_MF.mMSE[[j]]*pi[[j]]*pinv_MF.mMSE))
      })
      V_Temp<-lapply(1:length(All_Combinations),function(j){
        t(x_MF.mMSE[[j]]) %*% (x_MF.mMSE[[j]]*((as.vector(y_MF.mMSE)-pi[[j]])*pinv_MF.mMSE)^2)
      }) 
      V_Final<-lapply(1:length(All_Combinations),function(j){
        Mx[[j]] %*% V_Temp[[j]] %*% Mx[[j]]
      }) 
      
      for (j in 1:length(All_Combinations)) 
      {
        Utility_mMSE_MF[[j]][i,]<-cbind(r2[i],tr(V_Final[[j]]),det(solve(V_Final[[j]])))
        Bias_mMSE_MF[[j]][i,]<-Cordeiro(XData=x_MF.mMSE[[j]],With_bias = beta.mMSE_MF[[j]][i,])    
      }
    }
    
    # Full_SP_Single<-cbind(X,do.call(cbind,PI_Single.mMSE),do.call(cbind,PI_Single.mVc))
    # Full_SP_MF<-cbind(X,PI_MF.mMSE,PI_MF.mVc)
    # 
    # for (j in 1:length(All_Combinations)) 
    # {
    #   assign(paste0("Sample_mVc_Single_",j),do.call(rbind,Sample.mVc_Single[[j]]))
    #   assign(paste0("Sample_mVc_MF_",j),do.call(rbind,Sample.mVc_MF[[j]]))
    # 
    #   assign(paste0("Sample_mMSE_Single_",j),do.call(rbind,Sample.mMSE_Single[[j]]))
    #   assign(paste0("Sample_mMSE_MF_",j),do.call(rbind,Sample.mMSE_MF[[j]]))
    # }
    
    for(j in 1:length(All_Combinations))
    {
      assign(paste0("opt_Single_",j),
             list("Est_Theta_mMSE"=cbind(r2,beta.mMSE_Single[[j]]),"Utility_mMSE"=Utility_mMSE_Single[[j]],
                  "Bias_mMSE"=cbind(r2,Bias_mMSE_Single[[j]]),
                  "Est_Theta_mVc"=cbind(r2,beta.mVc_Single[[j]]),"Utility_mVc"=Utility_mVc_Single[[j]],
                  "Bias_mVc"=cbind(r2,Bias_mVc_Single[[j]]))
             )
      assign(paste0("opt_MF_",j),
             list("Est_Theta_mMSE"=cbind(r2,beta.mMSE_MF[[j]]),"Utility_mMSE"=Utility_mMSE_MF[[j]],
                  "Bias_mMSE"=cbind(r2,Bias_mMSE_MF[[j]]),
                  "Est_Theta_mVc"=cbind(r2,beta.mVc_MF[[j]]),"Utility_mVc"=Utility_mVc_MF[[j]],
                  "Bias_mVc"=cbind(r2,Bias_mVc_MF[[j]]))
      )
    }
    
    msg <- NULL
    return(list(opt=list("opt_Single"=mget(paste0("opt_Single_",1:length(All_Combinations))),
                         "opt_MF"=mget(paste0("opt_MF_",1:length(All_Combinations)))), 
                msg=msg)) 
}
