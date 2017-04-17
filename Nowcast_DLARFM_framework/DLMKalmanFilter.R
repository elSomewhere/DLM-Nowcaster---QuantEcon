runKF <- function(y,A,C,Q,R,Z_0,V_0){
  
  #note that y is transposed - time is accross columns. I'll correct this later
  
  #C = Z = lambda = loadings in observation matrix = design matrix = OBSERVATION EQUATION
  #A = T = transitions matrix = STATE  EQUATION
  #Q = covariance matrix from VAR and AR = state disturbance cov mx. = STATE EQUATION
  #R = observation disturbance covariance mx. (very small) = OBSERVATION EQUATION
  #initZ = initial state vector - prior state (alpha0) mean at t=0
  #initV = initial variance of state vector - prior state (alpha0) cov at t=0
  #there is no element for the state (alpha) itself!
  S <- KFRun(y,A,C,Q,R,Z_0,V_0)
  S <- KSRun(y,C,R,A,Q,S)
  
  xSmooth <- S$alpha_t_T
  VSmooth <- S$cov_alpha_t_T
  VVSmooth <- S$cov_alpha_t_T1
  loglik <- S$logLik
  res <- list(xSmooth,VSmooth,VVSmooth,loglik)
  names(res) <- c("xSmooth","VSmooth","VVSmooth","loglik")
  return(res)
}
KFRun <- function(y,A,C,Q,R,Z_0,V_0){
  
  #note: in some papers:
  #A = T
  #C = Z
  #R = H
  #Q = Q
  #Z_0 <- alpha0
  #V_0 <- v0
  
  
  n <- dim(C)[1]
  m <- dim(C)[2]
  numObs <- dim(y)[2]
  
  #create empty output matrices
  alpha_t_tm1 <- matrix(NA,m,numObs)#Predicted state vector
  cov_alpha_t_tm1 <- array(rep(NA,m*m*numObs),dim=c(m,m,numObs))#Predicted state covariance
  alpha_t_t <- matrix(NA,m,numObs+1)#Filtered  state vector
  cov_alpha_t_t <- array(rep(NA,m*m*numObs),dim=c(m,m,numObs+1))#Filtered  state covariance - has an addiitonal element at t0
  logLik = 0
  
  
  #Alpha = state
  #P = state convariance
  alpha_temp <- Z_0 #init state mean
  cov_alpha_temp <- V_0 #init state covariance
  
  #set first observation as the initial state (at time t0...)
  alpha_t_t[,1] <- alpha_temp
  cov_alpha_t_t[,,1] <- cov_alpha_temp 
  
  
  
  for(t in 1:numObs){
    #print(t)
    alpha <- A%*%alpha_temp
    cov_alpha <- A%*%cov_alpha_temp%*%t(A) + Q
    cov_alpha <- 0.5 * (cov_alpha + t(cov_alpha)) #just double check diagnoality... cus this amtrix is diagonal, the whole operation results in itself
    
    select <- missingData(y[,t,drop=FALSE],C,R)
    y_t <- select$y_t
    C_t <- select$C
    R_t <- select$R
    
    if(getDim(y_t)==0){
      alpha_temp <- alpha
      cov_alpha_temp <- cov_alpha
    }else{
      eq_1 <- cov_alpha%*%t(C_t)
      eq_2 <- ginv(C_t%*%eq_1+R_t)
      eq_3 <- eq_1%*%eq_2
      
      #V= forecast eror
      V <- y_t-C_t%*%alpha
      
      alpha_temp <- alpha+eq_3%*%V
      cov_alpha_temp <- cov_alpha-eq_3%*%t(eq_1)
      cov_alpha_temp <- 0.5*(cov_alpha_temp+t(cov_alpha_temp))
      logLik <- logLik+0.5*(log(det(eq_2))-t(V)%*%eq_2%*%V)
    }
    alpha_t_tm1[,t] <- alpha
    cov_alpha_t_tm1[,,t] <- cov_alpha
    
    alpha_t_t[,t+1] <- alpha_temp
    cov_alpha_t_t[,,t+1] <- cov_alpha_temp
    
  }
  if(getDim(y_t)==0){
    kalmanGain <- matrix(0,m,m)
  }else{
    kalmanGain <- eq_3%*%C_t
  }
  
  res <- list(alpha_t_tm1,cov_alpha_t_tm1,alpha_t_t,cov_alpha_t_t,kalmanGain,logLik)
  names(res) <- c("alpha_t_tm1","cov_alpha_t_tm1","alpha_t_t","cov_alpha_t_t","kalmanGain","logLik")
  return(res)
}


#KFRunResult <- KFRun(y,A,C,Q,R,Z_0,V_0)
KSRun <- function(y,C,R,A,Q,KFRunResult){
  n <- dim(C)[1]
  m <- dim(C)[2]
  numObs <- dim(y)[2]
  
  #create empty output matrices
  alpha_t_T <- matrix(0,m,numObs+1)#Filtered  state vector
  cov_alpha_t_T <- array(rep(0,m*m*numObs),dim=c(m,m,numObs+1))#Filtered  state covariance - has an addiitonal element at t0
  alpha_t_T[,numObs+1] <- drop(KFRunResult$alpha_t_t[,numObs+1])
  cov_alpha_t_T[,,numObs+1] <- drop(KFRunResult$cov_alpha_t_t[,,numObs+1])
  
  cov_alpha_t_T1 <- array(rep(0,m*m*numObs),dim=c(m,m,numObs))
  cov_alpha_t_T1[,,numObs] <- (diag(m)-KFRunResult$kalmanGain)%*%A%*%drop(KFRunResult$cov_alpha_t_t[,,numObs])
  
  
  
  #here, we should explore the possibility of replacing ginv simply with solve. Why do we actually need ginv?
  
  J_2 <- drop(KFRunResult$cov_alpha_t_t[,,numObs])%*%t(A)%*%ginv(KFRunResult$cov_alpha_t_tm1[,,numObs])
  #J_2 <- drop(KFRunResult$cov_alpha_t_t[,,numObs])%*%t(A)%*%solve(KFRunResult$cov_alpha_t_tm1[,,numObs])
  
  for(t in numObs:1){
    cov_alpha_t_t_temp <- drop(KFRunResult$cov_alpha_t_t[,,t])
    cov_alpha_t_tm1_temp <- drop(KFRunResult$cov_alpha_t_tm1[,,t])
    cov_alpha_t_T_temp <- drop(cov_alpha_t_T[,,t+1])
    cov_alpha_t_T1_temp <- drop(cov_alpha_t_T1[,,t])
    
    J_1 <- J_2
    
    alpha_t_T[,t] <- KFRunResult$alpha_t_t[,t] + J_1 %*% (alpha_t_T[,t+1]-A%*%KFRunResult$alpha_t_t[,t])
    cov_alpha_t_T[,,t] <-  cov_alpha_t_t_temp + J_1 %*% (cov_alpha_t_T_temp-cov_alpha_t_tm1_temp) %*% t(J_1)
    
    if(t>1){
      J_2 <- KFRunResult$cov_alpha_t_t[,,t-1]%*%t(A)%*%ginv(drop(KFRunResult$cov_alpha_t_tm1[,,t-1]))
      #J_2 <- KFRunResult$cov_alpha_t_t[,,t-1]%*%t(A)%*%solve(drop(KFRunResult$cov_alpha_t_tm1[,,t-1]))
      cov_alpha_t_T1[,,t-1] <- cov_alpha_t_t_temp%*%t(J_2)+J_1%*%(cov_alpha_t_T1_temp-(A%*%cov_alpha_t_t_temp))%*%t(J_2)
    }
  }
  res <- KFRunResult
  res[[7]] <- alpha_t_T
  res[[8]] <- cov_alpha_t_T
  res[[9]] <- cov_alpha_t_T1
  names(res)[c(7,8,9)] <- c("alpha_t_T","cov_alpha_t_T","cov_alpha_t_T1")
  return(res)
}
missingData <- function(y_t,C,R){
  ix <- !is.na(y_t)
  e <- diag(getDim(ix)[1])
  L <- e[,ix,drop=FALSE]
  y_t <- y_t[ix,drop=FALSE]
  C <- C[ix,,drop=FALSE]
  R <- R[ix,ix,drop=FALSE]
  res <- list(y_t,C,R,L)
  names(res) <- c("y_t","C","R","L")
  return(res)
}