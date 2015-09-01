
cal_beta_loglik = function(para,Z.aug,Y,subject.n,time.n,
                               Z.test.coeff.index,
                               quad.n=30){
  Y.tem <- Y
  s2 <- para[1]
  v  <- para[2]
  #### Z.test.coeff.index == TRUE means we want to test this parameter
  #### Therefore we want to keep it as 0 for LRT
  beta <- rep(NA,length(Z.test.coeff.index))
  beta[Z.test.coeff.index]  <- 0 ## initialized as 0
  beta[!Z.test.coeff.index] <- para[-(1:2)][1:sum(!Z.test.coeff.index)]
  beta <- as.matrix(beta)
  ####
  gherm <- generate_gaussian_quad_points(quad.n = quad.n)
  gh.weights <- matrix(rep(gherm$weights,subject.n),nrow = subject.n,byrow = TRUE)
  gh.nodes <- matrix(rep(gherm$nodes,subject.n * time.n),
                     nrow = subject.n * time.n,byrow = TRUE)
  #browser()
  u <- 1 / (1 + exp(-(Z.aug %*% beta[,rep(1,quad.n)] + gh.nodes * s2 * sqrt(2))))
  #### replace Y==0 with NA, so don't use them in the likelihood calculation
  Y.tem[Y == 0] <- NA
  #### log likelihood
  logL <- sum(log(rowSums(
    gh.weights / sqrt(pi) *
      apply(u,2, ## for each quad point
            function(u_m) {
              apply(matrix(gamma(v)/(gamma(u_m*v)*gamma((1-u_m)*v))*Y.tem^(u_m*v-1)*(1-Y.tem)^((1-u_m)*v-1),nrow = subject.n,ncol = time.n,byrow = TRUE),1,prod,na.rm=TRUE)
            }#function
      ) # apply
    ,na.rm=TRUE)),na.rm=TRUE)
  return(-logL)
}


#######################################
fit_beta_random_effect = function(Z=Z,Y=Y,
                                      subject.ind=subject.ind,time.ind=time.ind,
                                      quad.n=30,verbose=FALSE){

  ######
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)
  Z.aug <- cbind(intersept = 1, Z)
  ######
  est.table <- matrix(NA,ncol=2,nrow=ncol(Z.aug),
                      dimnames=list(colnames(Z.aug),c('Estimate','Pvalue')))
  #############
  subject.n <- length(unique(subject.ind))
  time.n   <- length(unique(time.ind))
  #### re-order Z,Y so that values belongs to the same subject are together
  #### need the values in this format for loglikelihood calculation
  gind <- sort(subject.ind,index.return=TRUE)$ix
  Y <- Y[gind,]
  Z <- Z[gind,]
  Z.aug <- Z.aug[gind,]
  ##### H1: estimate all parameters
  Z.test.coeff.index <- rep(FALSE,ncol(Z.aug))
  opt.H1 <- nlminb(start=c(1,2,rep(0,sum(!Z.test.coeff.index))), ## s2,v,alpha
                   objective=cal_beta_loglik,
                   lower = c(0.00001,0.00001,
                             rep(-Inf,sum(!Z.test.coeff.index))),
                   upper = Inf,
                   Z.test.coeff.index = Z.test.coeff.index,
                   Y=Y,Z.aug=Z.aug,time.n=time.n,subject.n=subject.n,
                   quad.n=quad.n,
                   control=list(trace=ifelse(verbose,2,0))
  )
  #### save the estimatd results
  s2.est   <- opt.H1$par[1]
  v.est    <- opt.H1$par[2]
  beta.est <- opt.H1$par[-(1:2)]
  est.table[,'Estimate'] <- beta.est
  ########################
  ####### H0
  for (test.i in 1:ncol(Z.aug)){
    Z.test.coeff.index <- rep(FALSE,ncol(Z.aug))
    Z.test.coeff.index[test.i] <- TRUE
    opt.H0 <- nlminb(start=c(1,2,rep(0,sum(!Z.test.coeff.index))), ## s2,alpha
                     objective=cal_beta_loglik,
                     lower = c(0.00001,0.00001,
                               rep(-Inf,sum(!Z.test.coeff.index))),
                     upper = c(Inf),
                     Z.test.coeff.index = Z.test.coeff.index,
                     Y=Y,Z.aug=Z.aug,time.n=time.n,subject.n=subject.n,
                     quad.n=quad.n,
                     control=list(trace=ifelse(verbose,2,0))
    )
    #### likelihood ratio test
    likelihodd.ratio <- -2*(-opt.H0$objective-(-opt.H1$objective))
    LRT.p <- 1-pchisq(likelihodd.ratio,df=1)
    est.table[test.i,'Pvalue'] <- LRT.p
  }

  return(list(est.table=est.table,s2.est=s2.est,v.est=v.est))
}
