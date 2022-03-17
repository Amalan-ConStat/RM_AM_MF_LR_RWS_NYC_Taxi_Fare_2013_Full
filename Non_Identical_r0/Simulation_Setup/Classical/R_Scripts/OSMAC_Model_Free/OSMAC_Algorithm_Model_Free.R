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

# Two Step OSMAC ----
OSMAC_MF <- function(r1,r2,Y,X,n,alpha,combs,All_Covariates) {
  PI.prop <- rep(1/n, n)
  idx.prop <- sample(1:n, r1, T, PI.prop)
  
  y.prop <- Y[idx.prop]
  pinv.prop <- rep(n,r1)
  
  beta.prop<-list()
  beta.prop<-lapply(1:length(combs),function(i) {
    glm(y.prop~X[idx.prop,colnames(X) %in% combs[[i]]]-1,family = "poisson")})
  
  for (i in 1:length(combs)) 
  {
    if (anyNA(beta.prop[[i]]$coefficients))
    {
      stop("There are NA or NaN values")
    }
  }
  
  P.prop<-list()
  P.prop<-lapply(1:length(combs), function(i)
    exp(X[,colnames(X) %in% combs[[i]]] %*% beta.prop[[i]]$coefficients))
  
  beta.mVc<-Utility_mVc<-Bias_mVc<-list()
  beta.mMSE<-Utility_mMSE<-Bias_mMSE<-list()
  
  for (a in 1:length(combs)) 
  {
    beta.mVc[[a]]<-matrix(nrow = length(r2),ncol = length(combs[[a]])) #
    Utility_mVc[[a]]<-matrix(nrow = length(r2),ncol = 3 ) #
    Bias_mVc[[a]]<-matrix(nrow = length(r2),ncol = length(combs[[a]])) #
    
    beta.mMSE[[a]]<-matrix(nrow = length(r2),ncol = length(combs[[a]])) #
    Utility_mMSE[[a]]<-matrix(nrow = length(r2),ncol = 3) #
    Bias_mMSE[[a]]<-matrix(nrow = length(r2),ncol = length(combs[[a]]))  #
  }
  
  #Sample.mMSE<-list()
  #Sample.mVc<-list()
  
  for (i in 1:length(r2)) 
  {
    ## mVC
    PI.mVc<-list()
    PI.mVc<-lapply(1:length(combs),function(j){
      PI.mVc <- sqrt((Y - P.prop[[j]])^2 * rowSums(X[,All_Covariates %in% combs[[j]] ]^2)) #
      return(PI.mVc / sum(PI.mVc))})
    
    PI.mVc<-matrix(unlist(PI.mVc),nrow = n,byrow = FALSE)
    pjoin.mVc<-rowWeightedMeans(PI.mVc,w=alpha)
    
    idx.mVc <- sample(1:n, r2[i]-r1, T, pjoin.mVc) # ##
    
    fit.mVc<-list()
    fit.mVc<-lapply(1:length(combs),function(i){
      x.mVc <- X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[i]] ] #
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv.mVc <- c(1 / pjoin.mVc[idx.mVc], pinv.prop)
      return(glm(y.mVc~x.mVc-1,weights=pinv.mVc,family = "poisson"))
    })
    
    # Sample.mVc[[i]]<-cbind(r2=r2[i],Y=Y[c(idx.mVc, idx.prop)],
    #                        SP=1/c(1 / pjoin.mVc[idx.mVc], pinv.prop),
    #                        X[c(idx.mVc, idx.prop),])
    
    for (j in 1:length(combs))
    {
      ru <- length(Y[c(idx.mVc, idx.prop)])
      beta.mVc[[j]][i,] <- fit.mVc[[j]]$coefficients
      Bias_mVc[[j]][i,]<-Cordeiro(XData=X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ],
                                  With_bias = beta.mVc[[j]][i,])
      
      pi <- c(exp(X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ] %*% beta.mVc[[j]][i,])) #
      pinv_join.mVc<-c(1 / pjoin.mVc[idx.mVc], pinv.prop)
      Mx <- solve(t(X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ]) %*% 
                   (X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ]*pi*pinv_join.mVc))
      V_Temp<- t(X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ]) %*% 
                (X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ]*
                   ((as.vector(Y[c(idx.mVc, idx.prop)])-pi)*pinv_join.mVc)^2)#
      V_Final<-Mx %*% V_Temp %*% Mx #
      
      Utility_mVc[[j]][i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final))) #
    }
    
    ## mMSE
    p.prop <-lapply(1:length(combs),function(i) P.prop[[i]][idx.prop])  #
    w.prop <-p.prop  #
    
    W.prop <-lapply(1:length(combs),function(i) {
      solve(t(X[idx.prop,All_Covariates %in% combs[[i]] ]) %*% 
             (X[idx.prop,All_Covariates %in% combs[[i]] ] * w.prop[[i]] * pinv.prop)) })   #
    
    PI.mMSE <-lapply(1:length(combs), function(i){
      PI.mMSE<-sqrt((Y - P.prop[[i]])^2 *rowSums((X[,All_Covariates %in% combs[[i]] ]%*%W.prop[[i]])^2)) 
      return(PI.mMSE / sum(PI.mMSE)) } )   #
    
    PI.mMSE<-matrix(unlist(PI.mMSE),nrow = n,byrow = FALSE)
    pjoin.mMSE<-rowWeightedMeans(PI.mMSE,w=alpha)
    
    idx.mMSE <- sample(1:n, r2[i]-r1, T, pjoin.mMSE)
    
    fit.mMSE<-list()
    fit.mMSE<-lapply(1:length(combs),function(i){
      x.mMSE <- X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[i]] ] #
      y.mMSE <- Y[c(idx.mMSE, idx.prop)]
      pinv.mMSE <- c(1 / pjoin.mMSE[idx.mMSE], pinv.prop)
      return(glm(y.mMSE~x.mMSE-1, weights=pinv.mMSE,family = "poisson"))
    })
    
    # Sample.mMSE[[i]]<-cbind(r2=r2[i],Y=Y[c(idx.mMSE, idx.prop)],
    #                         SP=1/c(1 / pjoin.mMSE[idx.mMSE], pinv.prop),
    #                         X[c(idx.mMSE, idx.prop),])
    
    for (j in 1:length(combs)) 
    {
      ru <- length(Y[c(idx.mMSE, idx.prop)])
      beta.mMSE[[j]][i,] <-fit.mMSE[[j]]$coefficients
      Bias_mMSE[[j]][i,]<-Cordeiro(XData=X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ],
                                   With_bias = beta.mMSE[[j]][i,])
      
      pi <-c(exp(X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ] %*% beta.mMSE[[j]][i,])) #
      pinv_join.mMSE<-c(1 / pjoin.mMSE[idx.mMSE], pinv.prop)
      Mx <- solve(t(X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ]) %*% 
                   (X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ]*pi*pinv_join.mMSE)) 
      V_Temp<- t(X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ]) %*% 
                (X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ]*
                   ((as.vector(Y[c(idx.mMSE, idx.prop)])-pi)*pinv_join.mMSE)^2) #
      V_Final<-Mx %*% V_Temp %*% Mx #
      
      Utility_mMSE[[j]][i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final))) #
    }
  }
  
  # Sample_mVc<-do.call(rbind,Sample.mVc) #
  # Sample_mMSE<-do.call(rbind,Sample.mMSE) #
  opt<-list()
  
  for (i in 1:length(combs)) 
  {
    opt[[i]] <- list("Est_Theta_mMSE"=cbind(r2,beta.mMSE[[i]]),"Utility_mMSE"=Utility_mMSE[[i]],
                     "Bias_mMSE"=cbind(r2,Bias_mMSE[[i]]),
                     "Est_Theta_mVc"=cbind(r2,beta.mVc[[i]]),"Utility_mVc"=Utility_mVc[[i]],
                     "Bias_mVc"=cbind(r2,Bias_mVc[[i]]))
    if(anyNA(opt[[i]]$Est_Theta_mMSE) || anyNA(opt[[i]]$Est_Theta_mVc) )
    {
      stop("There are NA or NaN values")
    }
    
    if(anyNA(opt[[i]]$Utility_mMSE) || anyNA(opt[[i]]$Utility_mVc) )
    {
      stop("There are NA or NaN values")
    }
    
    if(anyNA(opt[[i]]$Bias_mMSE) || anyNA(opt[[i]]$Bias_mVc) )
    {
      stop("There are NA or NaN values")
    }
  }
  # if(anyNA(opt$Sample_mMSE) || anyNA(opt$Sample_mVc) )
  # {
  #   stop("There are NA or NaN values")
  # }
  return(list("opt"=opt#,"Sample_mMSE"=Sample_mMSE,"Sample_mVc"=Sample_mVc
              ))
}
