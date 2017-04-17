initialize_theta <- function(data_scaled,params_init,preInitR,preInitSteps){
  #now, the spanish nowcaster simply assumes:
  #identity matrix for preliminary state cov V
  #0s for preliminary state mean Z
  #identity matrix for R
  
  
  data_scaled=data_scaled
  r=params_init[[1]]
  p=params_init[[2]]
  g=params_init[[3]]
  index_idio_m=params_init[[4]]
  index_idio_q=params_init[[5]]
  restrictions_matrix=params_init[[6]]
  restrictions_equalizer=params_init[[7]]
  
  #spline data - later, we should also regress fill in data backwards to have a ore complete balanced set
  #we use the hard timmer here - because the current spline function has problem with leading and trailing NAs
  #correct this! use some kind of smart interpolation or so
  
  #data_scaled_trimmed <- trimmer_hard(data_scaled)
  
  data_scaled_trimmed <- data_scaled
  #data_scaled_trimmed <- trimmer_soft(data_scaled,1)
  #res_spliner <- spliner(data_scaled_trimmed)
  res_spliner <- spliner_staticModel(data_scaled_trimmed,r=1,max_iter=preInitSteps,params_init)
  
  
  data_scaled_splined_trimmed_filled <- res_spliner[[1]]
  
  
  #some vars
  periods_data_scaled_splined <- dim(data_scaled_splined_trimmed_filled)[1]
  nQM <- dim(data_scaled_splined_trimmed_filled)[2]
  nM <- dim(index_idio_m)[2]
  nQ <- dim(index_idio_q)[2]
  size_f_normal <- p+1
  size_f_restriction <- dim(restrictions_matrix)[2]
  if(!is.null(nQ)){
    size_f_final <- max(size_f_normal,size_f_restriction)
  }else{
    size_f_final <- size_f_normal
  }
  
  
  
  
  #########################
  #monthly factor loadings#
  #########################
  #fill in loadings for the monthly part - simply the eigenvectors of the covariance matrix prescaled data - this is the definition of a pca!
  #the eigenvectors are actually only retrieved from the monthly series
  decompose <- eigen(cov(data_scaled_splined_trimmed_filled[,index_idio_m]))
  eigenvectors <- decompose$vectors[,1:r,drop=FALSE]
  #eigenvectors <- eigenvectors[,ncol(eigenvectors):1] #!!!!!!!!!!!Check how the output of eigenvectors really is.... do we actually have to reverse?
  C_cf_m <- eigenvectors
  
  ##############################################
  #Create factors matrix, incl. the lagged ones#
  ##############################################
  #we multiply original data with the loadings to get the commonfactors
  cf <- as.matrix(data_scaled_splined_trimmed_filled[,index_idio_m,drop=FALSE])%*%eigenvectors
  #this next series is shorter due to the lags - since we are only using it to create loadings, this is not much of an issue
  #this is the lagged factors matrix - we need this to fill up quarterly loadings, 
  #because we will be creating quarterly loadings by regressing to the factors extracted from before
  #this next series is shorter due to the lags - since we are only using it to create loadings, this is not much of an issue
  if(!is.null(nQ)){
    shortenCfBy <- size_f_restriction-1
    cf_lagged <- NULL
    for(i in 1:(size_f_restriction-1)){
      cf_lagged <- cbind(cf_lagged,takeLag(cf,(-1)*i,NA,FALSE)[-seq(shortenCfBy),])
    }
    cf_todayAndLagged <- cbind(cf[-seq(shortenCfBy),],cf_lagged)
    cf_lagged <- NULL
    for(i in 1:(size_f_restriction-1)){
      cf_lagged <- cbind(cf_lagged,takeLag(cf,(-1)*i,NA,FALSE)[-seq(shortenCfBy),])
    }
    cf_todayAndLagged <- cbind(cf[-seq(shortenCfBy),],cf_lagged)
  }else{
    shortenCfBy <- size_f_normal-1
    cf_lagged <- NULL
    for(i in 1:(size_f_normal-1)){
      cf_lagged <- cbind(cf_lagged,takeLag(cf,(-1)*i,NA,FALSE)[-seq(shortenCfBy),])
    }
    cf_todayAndLagged <- cbind(cf[-seq(shortenCfBy),],cf_lagged)
    cf_lagged <- NULL
    for(i in 1:(size_f_normal-1)){
      cf_lagged <- cbind(cf_lagged,takeLag(cf,(-1)*i,NA,FALSE)[-seq(shortenCfBy),])
    }
    cf_todayAndLagged <- cbind(cf[-seq(shortenCfBy),],cf_lagged)
  }
  
  
  
  ###########################
  #quarterly factor loadings#
  ###########################
  #we extend the restrictionstrix by the amount of quarterly values
  
  if(!is.null(nQ)){
    restrictions_matrix_extended_by_r <- kronecker(restrictions_matrix,diag(r))
    restrictions_equalizer_extended_by_r <- kronecker(restrictions_equalizer,matrix(0,r,1))
    C_cf_q <- matrix(NA,nQ,size_f_restriction*r)
    for(j in 1:nQ){
      data_for_q <- data_scaled_splined_trimmed_filled[-seq(shortenCfBy),index_idio_q[j],drop=FALSE]
      a <- solve(t(cf_todayAndLagged)%*%cf_todayAndLagged)
      b <- restrictions_matrix_extended_by_r
      c <- simpleMVRegression(data_for_q,cf_todayAndLagged)$beta
      d <- restrictions_equalizer_extended_by_r
      C_cf_q[j,] <- t(c-(a%*%t(b)%*%solve(b%*%a%*%t(b))%*%(b%*%c-d)))
    }
  }
  
  
  
  ##########
  #finish C#
  ##########
  
  if(!is.null(nQ)){
    #we need C  this - we temporarly construct it. In the end, we let the utility functions do this for consistency reasons. the resuls the same
    C_cf <- rbind(cbind(C_cf_m,matrix(0,nM,((size_f_final-1)*r))),C_cf_q)
    #the final part of C with the restrictions
    C_m <- rbind(kronecker(diag(nM),matrix(c(1,rep(0,g-1)),1,g)),matrix(0,nQ,g*nM))
    C_q <- rbind(matrix(0,nM,5*nQ),kronecker(diag(nQ),matrix(c(1,2,3,2,1),1,5)))
    C <- cbind(C_cf,C_m,C_q)
  }else{
    #we need C  this - we temporarly construct it. In the end, we let the utility functions do this for consistency reasons. the resuls the same
    C_cf <- rbind(cbind(C_cf_m,matrix(0,nM,((size_f_final-1)*r))))
    #the final part of C with the restrictions
    C_m <- rbind(kronecker(diag(nM),matrix(c(1,rep(0,g-1)),1,g)))
    C <- cbind(C_cf,C_m)
  }
  
  
  
  
  
  ################
  #Covar matrix R#
  ################
  if(!is.null(nQ)){
    #create the residuals, originaldata-commonfactors*loadings
    residuals_scaled_splined_trimmed_filled <- data_scaled_splined_trimmed_filled[-seq(shortenCfBy),]-cf_todayAndLagged%*%t(C_cf)
    #covar matrix of residuals (data minus the loadings*commonfactors), use the ones with holes to not introduce distortion
    R <- matrix(0,nQM,nQM)
    diag(R) <- diag(cov(residuals_scaled_splined_trimmed_filled))
    #I think this part is discarded again in the end... cus we set R very very small
    diag(R) <- 1e-04
  }else{
    #create the residuals, originaldata-commonfactors*loadings
    residuals_scaled_splined_trimmed_filled <- data_scaled_splined_trimmed_filled[-seq(shortenCfBy),]-cf_todayAndLagged%*%t(C_cf)
    #covar matrix of residuals (data minus the loadings*commonfactors), use the ones with holes to not introduce distortion
    R <- matrix(0,nQM,nQM)
    diag(R) <- diag(cov(residuals_scaled_splined_trimmed_filled))
    #I think this part is discarded again in the end... cus we set R very very small
    diag(R) <- 1e-04
  }
  
  
  
  ########################################
  #VAR betas and sigmas for commonfactors#
  ########################################
  cf_var <- simpleVAR(cf,p,constant=FALSE)
  A_cf <- rbind(cbind(cf_var$B_hat_UR,matrix(0,r,size_f_final*r-r*p)),cbind(diag(r*(size_f_final-1)),matrix(0,(size_f_final-1)*r,r)))
  Q_cf <- fillDiag(cf_var$sigma_UR,matrix(0,(size_f_final-1)*r,(size_f_final-1)*r),fill=0)
  
  
  ##################
  #alphas for idios#
  ##################
  alphas_m <- matrix(0,nM*g,nM*g)
  cov_idios_m <- matrix(0,nM*g,nM*g)
  pers <- dim(residuals_scaled_splined_trimmed_filled)[1]
  for (i in 1:nM){
    res_m <- matrix(residuals_scaled_splined_trimmed_filled[,i,drop=FALSE])
    
    a_t0 <- res_m[(1+g):pers,]
    a_tml <- matrix(NA,pers-g,g)
    for(j in 1:g){
      a_tml[,j] <- res_m[(1+g-j):(pers-j),]
    }
    
    
    alpha_res_i <- solve(t(a_tml)%*%a_tml)%*%t(a_tml)%*%a_t0
    alphas_m[((i-1)*g+1),((i-1)*g+1):(i*g)] <- t(alpha_res_i)
    if(g==2){
      alphas_m[((i-1)*g+g),((i-1)*g+1)] <- 1 
    }else if(g>2){
      diag(alphas_m[((i-1)*g+2):((i-1)*g+g),((i-1)*g+1):(i*g)]) <- 1  
    }
    cov_res_i <- cov(a_t0-a_tml%*%alpha_res_i)
    cov_idios_m[((i-1)*g+1),((i-1)*g+1)] <- cov_res_i
  }
  
  
  if(!is.null(nQ)){
    alphas_q <- matrix(0,nQ*5,nQ*5)
    cov_idios_q <- matrix(0,nQ*5,nQ*5)
    for (i in 1:nQ){
      res_m <- matrix(residuals_scaled_splined_trimmed_filled[,i,drop=FALSE])
      
      a_t0 <- res_m[(1+g):pers,]
      a_tml <- matrix(NA,pers-g,g)
      for(j in 1:g){
        a_tml[,j] <- res_m[(1+g-j):(pers-j),]
      }
      
      
      alpha_res_i <- solve(t(a_tml)%*%a_tml)%*%t(a_tml)%*%a_t0
      alphas_q[((i-1)*5+1),((i-1)*5+1):((i-1)*5+g)] <- t(alpha_res_i)
      diag(alphas_q[((i-1)*5+2):((i-1)*5+5),((i-1)*5+1):((i-1)*5+4)]) <- 1  
      cov_res_i <- cov(a_t0-a_tml%*%alpha_res_i)
      cov_idios_q[((i-1)*5+1),((i-1)*5+1)] <- cov_res_i
    }
  }
  
  if(!is.null(nQ)){
    A <- fillDiag(A_cf,fillDiag(alphas_m,alphas_q,fill=0),fill=0)
    Q <- fillDiag(Q_cf,fillDiag(cov_idios_m,cov_idios_q,fill=0),fill=0)  
  }else{
    A <- fillDiag(A_cf,alphas_m,fill=0)
    Q <- fillDiag(Q_cf,cov_idios_m,fill=0)  
  }
  
  
  
  
  ####################################################
  #V matrix - I think for the the transition equation#
  ####################################################
  #V_cf <- matrix(solve(diag((r*size_f_final)^2)-kronecker(A_cf,A_cf))%*%matrix(Q_cf,dim(Q_cf)[1]*dim(Q_cf)[2],1),r*size_f_final,r*size_f_final)
  
  #V1
  #V_m <- diag(1/diag(diag(nM*g)-(alphas_m^2)))*cov_idios_m
  #V_q <- (solve(diag((5*nQ)^2)-kronecker(A_restrictionspart,A_restrictionspart))%*%vecFun(Q_restrictionspart),5*nQ,5*nQ)
  
  #V2
  #V_m <- matrix(solve(diag((g*nM)^2)-kronecker(alphas_m,alphas_m))%*%vecFun(cov_idios_m),g*nM,g*nM)
  #V_q <- matrix(solve(diag((5*nQ)^2)-kronecker(alphas_q,alphas_q))%*%vecFun(cov_idios_q),5*nQ,5*nQ)
  
  #V <- fillDiag(V_cf,(fillDiag(V_m,V_q,fill=0)),fill=0)
  if(!is.null(nQ)){
    V <- diag(size_f_final*r+nM*g+nQ*5)
  }else{
    V <- diag(size_f_final*r+nM*g)
  }
  
  
  
  #################################
  #Z matrix - preliminarytate mean#
  #################################
  Z <- matrix(0,dim(A)[1],1)
  
  
  
  
  
  
  res <- list(A,C,Q,R,Z,V)
  
  
  #C = Z = lambda = loadings in observation matrix = design matrix = OBSERVATION EQUATION
  #A = T = transitions matrix = STATE  EQUATION
  #Q = covariance matrix from VAR and AR = state disturbance cov mx. = STATE EQUATION
  #R = observation disturbance covariance mx. (very small) = OBSERVATION EQUATION
  
  #initZ = initial state vector - prior state (alpha0) mean at t=0
  #initV = initial variance of state vector - prior state (alpha0) cov at t=0
  
  #there is no element for the state (alpha) itself!
  names(res) <- cbind("A","C","Q","R","initZ","initV")
  return(res)
}