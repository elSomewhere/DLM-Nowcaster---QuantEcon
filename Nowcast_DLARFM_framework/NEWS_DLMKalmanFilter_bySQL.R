#this re-runs the kalman filter but by augmenting the input variables with the carry-over of lags
#I believe this is essentially the same as running the KF on a lag of the state vector
#the difference is essentially only that we operate via lists instead of a 3d-array for covariance stuff

augmentLags <- function(y,A,C,Q,R,Z_0,V_0,lags){
  n <- dim(C)[1]
  r <- dim(C)[2]
  if(lags>0){
    C <- cbind(C,matrix(0,n,lags*r))
    A <- fillDiag(A,matrix(0,lags*r,lags*r),0)
    A[(r+1):dim(A)[1],1:(lags*r)] <- diag(lags*r)
    Q <- fillDiag(Q,matrix(0,lags*r,lags*r),0)
    Z_0 <- rbind(Z_0,matrix(0,lags*r,1))
    V_0 <- fillDiag(V_0,matrix(0,lags*r,lags*r),0)
  }
  res <- list(y,A,C,Q,R,Z_0,V_0)
  names(res) <- c("y","A","C","Q","R","Z_0","V_0")
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
KFKSRun_SQL <- function(y,A,C,Q,R,Z_0,V_0,location=NULL){
  if(is.null(location)){
    location <- getwd()
  }
  setwd(location)
  drivee <- SQLite()
  
  if(file.exists("kfdm.sqlite")){
    print("deleting existing KFDM file")
    file.remove("kfdm.sqlite")
  }
  
  
  con <- dbConnect(drivee, dbname="kfdm.sqlite")
  
  
  n <- dim(C)[1]
  m <- dim(C)[2]
  numObs <- dim(y)[2]
  
  #create empty output matrices
  #alpha_t_tm1 <- matrix(NA,m,numObs)#Predicted state vector
  dbGetQuery(con, "CREATE TABLE IF NOT EXISTS alpha_t_tm1 (_id integer PRIMARY KEY, data blob)")
  for(i in 1: numObs){
    df <- data.frame(t=i,d=I(list(serialize(matrix(NA,m),NULL))))
    dbGetPreparedQuery(con, 'insert into alpha_t_tm1 (_id, data) values (:t,:d)', bind.data=df)
  }
  
  #
  dbGetQuery(con, "CREATE TABLE IF NOT EXISTS cov_alpha_t_tm1 (_id integer PRIMARY KEY, data blob)")
  for(i in 1: numObs){
    df <- data.frame(t=i,d=I(list(serialize(matrix(NA,m,m),NULL))))
    dbGetPreparedQuery(con, 'insert into cov_alpha_t_tm1 (_id, data) values (:t,:d)', bind.data=df)
  }
  
  
  alpha_t_t <- matrix(NA,m,numObs+1)#Filtered  state vector
  dbGetQuery(con, "CREATE TABLE IF NOT EXISTS alpha_t_t (_id integer PRIMARY KEY, data blob)")
  for(i in 1: (numObs+1)){
    df <- data.frame(t=i,d=I(list(serialize(matrix(NA,m),NULL))))
    dbGetPreparedQuery(con, 'insert into alpha_t_t (_id, data) values (:t,:d)', bind.data=df)
  }
  
  dbGetQuery(con, "CREATE TABLE IF NOT EXISTS cov_alpha_t_t (_id integer PRIMARY KEY, data blob)")
  for(i in 1: (numObs+1)){
    df <- data.frame(t=i,d=I(list(serialize(matrix(NA,m,m),NULL))))
    dbGetPreparedQuery(con, 'insert into cov_alpha_t_t (_id, data) values (:t,:d)', bind.data=df)
  }
  
  logLik = 0
  
  
  #Alpha = state
  #P = state convariance
  alpha_temp <- Z_0 #init state mean
  cov_alpha_temp <- V_0 #init state covariance
  
  #set first observation as the initial state (at time t0...)
  #alpha_t_t[,1] <- alpha_temp
  df <- data.frame(t=1,d=I(list(serialize(alpha_temp,NULL))))
  dbGetPreparedQuery(con, 'UPDATE alpha_t_t SET data = (:d) WHERE _id = (:t)', bind.data=df)
  
  df <- data.frame(t=1,d=I(list(serialize(cov_alpha_temp,NULL))))
  dbGetPreparedQuery(con, 'UPDATE cov_alpha_t_t SET data = (:d) WHERE _id = (:t)', bind.data=df)
  
  
  for(t in 1:numObs){
    print(paste0("doing kalman filter run, iteration", t))
    #flush.console()
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
      eq_2 <- solve(C_t%*%eq_1+R_t)
      eq_3 <- eq_1%*%eq_2
      
      #V= forecast eror
      V <- y_t-C_t%*%alpha
      
      alpha_temp <- alpha+eq_3%*%V
      cov_alpha_temp <- cov_alpha-eq_3%*%t(eq_1)
      cov_alpha_temp <- 0.5*(cov_alpha_temp+t(cov_alpha_temp))
      logLik <- logLik+0.5*(log(det(eq_2))-t(V)%*%eq_2%*%V)
    }
    #alpha_t_tm1[,t] <- alpha
    df <- data.frame(t=t,d=I(list(serialize(alpha,NULL))))
    dbGetPreparedQuery(con, 'UPDATE alpha_t_tm1 SET data = (:d) WHERE _id = (:t)', bind.data=df)
    
    #cov_alpha_t_tm1[[t]] <- cov_alpha
    df <- data.frame(t=t,d=I(list(serialize(cov_alpha,NULL))))
    dbGetPreparedQuery(con, 'UPDATE cov_alpha_t_tm1 SET data = (:d) WHERE _id = (:t)', bind.data=df)
    
    #alpha_t_t[,t+1] <- alpha_temp
    df <- data.frame(t=t+1,d=I(list(serialize(alpha_temp,NULL))))
    dbGetPreparedQuery(con, 'UPDATE alpha_t_t SET data = (:d) WHERE _id = (:t)', bind.data=df)
    #cov_alpha_t_t[[t+1]] <- cov_alpha_temp
    df <- data.frame(t=t+1,d=I(list(serialize(cov_alpha_temp,NULL))))
    dbGetPreparedQuery(con, 'UPDATE cov_alpha_t_t SET data = (:d) WHERE _id = (:t)', bind.data=df)
    
  }
  if(getDim(y_t)==0){
    kalmanGain <- matrix(0,m,m)
  }else{
    kalmanGain <- eq_3%*%C_t
  }
  
  
  #create empty output matrices
  #alpha_t_T <- matrix(0,m,numObs+1)#Filtered  state vector
  dbGetQuery(con, "CREATE TABLE IF NOT EXISTS alpha_t_Tbig (_id integer PRIMARY KEY, data blob)")
  for(i in 1: (numObs+1)){
    df <- data.frame(t=i,d=I(list(serialize(matrix(0,m),NULL))))
    dbGetPreparedQuery(con, 'insert into alpha_t_Tbig (_id, data) values (:t,:d)', bind.data=df)
  }
  
  
  dbGetQuery(con, "CREATE TABLE IF NOT EXISTS cov_alpha_t_tbig (_id integer PRIMARY KEY, data blob)")
  for(i in 1: numObs+1){
    df <- data.frame(t=i,d=I(list(serialize(matrix(NA,m,m),NULL))))
    dbGetPreparedQuery(con, 'insert into cov_alpha_t_tbig (_id, data) values (:t,:d)', bind.data=df)
  }
  
  
  query <- paste0("SELECT * FROM alpha_t_t WHERE _id = ",numObs+1)
  pp <- dbGetQuery(con, query)
  pp <- unserialize(pp$data[[1]])
  df <- data.frame(t=numObs+1,d=I(list(serialize(pp,NULL))))
  dbGetPreparedQuery(con, 'UPDATE alpha_t_Tbig SET data = (:d) WHERE _id = (:t)', bind.data=df)
  
  
  
  query <- paste0("SELECT * FROM cov_alpha_t_t WHERE _id = ",numObs+1)
  pp <- dbGetQuery(con, query)
  pp <- unserialize(pp$data[[1]])
  
  df <- data.frame(t=numObs+1,d=I(list(serialize(pp,NULL))))
  dbGetPreparedQuery(con, 'UPDATE cov_alpha_t_tbig SET data = (:d) WHERE _id = (:t)', bind.data=df)
  
  #cov_alpha_t_T[[numObs+1]] <- KFRunResult$cov_alpha_t_t[[numObs+1]]
  
  
  
  
  
  
  
  
  
  
  
  dbGetQuery(con, "CREATE TABLE IF NOT EXISTS cov_alpha_t_T1 (_id integer PRIMARY KEY, data blob)")
  for(i in 1: numObs){
    df <- data.frame(t=i,d=I(list(serialize(matrix(NA,m,m),NULL))))
    dbGetPreparedQuery(con, 'insert into cov_alpha_t_T1 (_id, data) values (:t,:d)', bind.data=df)
  }
  query <- paste0("SELECT * FROM cov_alpha_t_t WHERE _id = ",numObs)
  pp1 <- dbGetQuery(con, query)
  pp1 <- unserialize(pp1$data[[1]])
  forin1 <-   (diag(m)-kalmanGain)%*%A%*%(pp1)
  df <- data.frame(t=numObs,d=I(list(serialize(forin1,NULL))))
  dbGetPreparedQuery(con, 'UPDATE cov_alpha_t_T1 SET data = (:d) WHERE _id = (:t)', bind.data=df)
  
  
  
  query <- paste0("SELECT * FROM cov_alpha_t_tm1 WHERE _id = ",numObs)
  pp2 <- dbGetQuery(con, query)
  pp2 <- unserialize(pp2$data[[1]])
  forin2 <-   pp1%*%t(A)%*%ginv(pp2)
  J_2 <- forin2
  
  
  for(t in numObs:1){
    print(paste0("doing augmented kalman run, currently at smoothing iteration ",t))
    
    
    
    
    query <- paste0("SELECT * FROM cov_alpha_t_t WHERE _id = ",t)
    cov_alpha_t_t_temp <- dbGetQuery(con, query)
    cov_alpha_t_t_temp <- unserialize(cov_alpha_t_t_temp$data[[1]])
    
    query <- paste0("SELECT * FROM cov_alpha_t_tm1 WHERE _id = ",t)
    cov_alpha_t_tm1_temp <- dbGetQuery(con, query)
    cov_alpha_t_tm1_temp <- unserialize(cov_alpha_t_tm1_temp$data[[1]])
    
    query <- paste0("SELECT * FROM cov_alpha_t_tbig WHERE _id = ",t+1)
    cov_alpha_t_T_temp <- dbGetQuery(con, query)
    cov_alpha_t_T_temp <- unserialize(cov_alpha_t_T_temp$data[[1]])
    
    query <- paste0("SELECT * FROM cov_alpha_t_T1 WHERE _id = ",t)
    cov_alpha_t_T1_temp <- dbGetQuery(con, query)
    cov_alpha_t_T1_temp <- unserialize(cov_alpha_t_T1_temp$data[[1]])
    
    
    J_1 <- J_2
    
    
    query <- paste0("SELECT * FROM alpha_t_t WHERE _id = ",t)
    aaa <- dbGetQuery(con, query)
    aaa <- unserialize(aaa$data[[1]])
    
    
    query <- paste0("SELECT * FROM alpha_t_Tbig WHERE _id = ",t+1)
    bbb <- dbGetQuery(con, query)
    bbb <- unserialize(bbb$data[[1]])
    
    ccc <- aaa + J_1 %*% (bbb-A%*%aaa)
    
    df <- data.frame(t=t,d=I(list(serialize(ccc,NULL))))
    dbGetPreparedQuery(con, 'UPDATE alpha_t_Tbig SET data = (:d) WHERE _id = (:t)', bind.data=df)
    
    forin <-  cov_alpha_t_t_temp + J_1 %*% (cov_alpha_t_T_temp-cov_alpha_t_tm1_temp) %*% t(J_1)
    df <- data.frame(t=t,d=I(list(serialize(forin,NULL))))
    dbGetPreparedQuery(con, 'UPDATE cov_alpha_t_tbig SET data = (:d) WHERE _id = (:t)', bind.data=df)
    
    
    
    if(t>1){
      
      query <- paste0("SELECT * FROM cov_alpha_t_t WHERE _id = ",t-1)
      a <- dbGetQuery(con, query)
      a <- unserialize(a$data[[1]])
      
      query <- paste0("SELECT * FROM cov_alpha_t_tm1 WHERE _id = ",t-1)
      b <- dbGetQuery(con, query)
      b <- unserialize(b$data[[1]])
      
      J_2 <- a%*%t(A)%*%ginv(b)
      
      pp <-  cov_alpha_t_t_temp%*%t(J_2)+J_1%*%(cov_alpha_t_T1_temp-(A%*%cov_alpha_t_t_temp))%*%t(J_2)
      df <- data.frame(t=t-1,d=I(list(serialize(pp,NULL))))
      dbGetPreparedQuery(con, 'UPDATE cov_alpha_t_T1 SET data = (:d) WHERE _id = (:t)', bind.data=df)
      
      
      forin <-  cov_alpha_t_t_temp%*%t(J_2)+J_1%*%(cov_alpha_t_T1_temp-(A%*%cov_alpha_t_t_temp))%*%t(J_2)
      df <- data.frame(t=t-1,d=I(list(serialize(forin,NULL))))
      dbGetPreparedQuery(con, 'UPDATE cov_alpha_t_T1 SET data = (:d) WHERE _id = (:t)', bind.data=df)
    }
  }
  
  
  
  #query <- paste0("SELECT * FROM alpha_t_Tbig")
  #alpha_t_T <- dbGetQuery(con, query)
  #alpha_t_T <- lapply(alpha_t_T$data,function(x){
   # res <- unserialize(x)
  #})
  #alpha_t_T <- do.call(cbind,alpha_t_T)
  
  res <- list(con,getwd())
  
  return(res)
  
  
  
  
  
  
}


runKF_SQL <- function(y,A,C,Q,R,Z_0,V_0){
  library(RSQLite)
  #note that y is transposed - time is accross columns. I'll correct this later
  
  #C = Z = lambda = loadings in observation matrix = design matrix = OBSERVATION EQUATION
  #A = T = transitions matrix = STATE  EQUATION
  #Q = covariance matrix from VAR and AR = state disturbance cov mx. = STATE EQUATION
  #R = observation disturbance covariance mx. (very small) = OBSERVATION EQUATION
  #initZ = initial state vector - prior state (alpha0) mean at t=0
  #initV = initial variance of state vector - prior state (alpha0) cov at t=0
  #there is no element for the state (alpha) itself!
  S <- KFKSRun_SQL(y,A,C,Q,R,Z_0,V_0)
 
  return(S)
 
}


