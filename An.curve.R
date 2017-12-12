an.com.H0 <- function(dat,tri=2:9){
  
  
  pheno1 <- dat$pheno1[,tri]/10
  pheno2 <- dat$pheno2[,tri]/10
  times <- dat$sample.time[tri]
  
  
  #parin2 <- c(1.721163,0.1969936,1.676018,0.1931441,0.544902,
              #4.247237,7.763093,14.14129,-3.682609,-0.2948541,
              #1.502972,7.579344,5.573393,-1.005777,-0.9886827)
  parin2 <- c(2.0469798,0.1274323,2.0159414,0.1214398,0.1323994,
              4.6961261,7.7807460,17.3646514,-4.1523418,-0.2554463,
              1.5665796,15.3576427,13.8041671,-1.1094356,-0.8260141)
  mpheno1 <- as.numeric(colMeans(pheno1,na.rm=T))
  mpheno2 <- as.numeric(colMeans(pheno2,na.rm=T))
  
  loop_k <- 1
  max_iter <- 100
  epsi <- 10^-4
  max_err <- 1
  
  while(loop_k<max_iter && max_err>epsi){
    
    
    oldpar <-c(parin2);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin2[6:15])
      AA <- curve.mle(nnpar,y=cbind(pheno1,pheno2),time.std =times,x0=mpheno1[1],y0=mpheno2[1])
      AA
    }
    r1.covar <- optim(parin2[1:5],mle.covar1,method = "BFGS",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle.1 <- function(npar){
      
      nnpar <- c(new.covar1,npar)
      AA <- curve.mle(nnpar,y=cbind(pheno1,pheno2),time.std =times,x0=mpheno1[1],y0=mpheno2[1])
      AA
    }
    r1 <- optim(c(parin2[6:15]),mle.1,method = "BFGS",control=list(maxit=32000))    
    new1 <- r1$par
    
    #r2<- optim(c(new.covar1,new1),curve.mle,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[2],
    #method = "BFGS",control=list(maxit=32000,trace=T))    
    #cat("new1:",unlist( new1), "\n");
    
    nparin <- c(new.covar1,new1)
    
    newpar <- c(nparin)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin2 <- nparin
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  LL <- curve.mle(parin2,cbind(pheno1,pheno2),time.std =times,x0=mpheno1[1],y0=mpheno2[1])
  return(c(LL,parin2))
}


curve.mle <-function( par,y,time.std,x0,y0)
{
  len.cov <- 5
  par.covar <- par[1:len.cov]
  n  <- length(y[,1])
  #sig.inv3 <- SAD3.get_inv_mat(par.covar,time.std, 2)
  cov.mat <- SAD3.get_mat(par.covar,time.std, 2)#solve( sig.inv3 )
  
  curve.par <- par[(len.cov+1):(len.cov+ 10)]
  mu <- com.get_mu(curve.par,time.std,x0=x0,y0=y0)
  
  if (any(is.na(mu) ))
    return (NaN);
  
  Y.delt <- y- matrix(rep(mu,n),nrow=n,byrow=T)
  
  nna.vec <- get_non_na_number(Y.delt)
  
  pv <- dmvnorm_fast( Y.delt, rep(0, NCOL(cov.mat)), cov.mat, nna.vec=nna.vec, log=T)
  
  A <- -sum( pv );
  return(A)
}





get_non_na_number<-function(y.resd)
{
  nna.mat <- array(1, dim=c(NROW(y.resd), NCOL(y.resd)))
  nna.mat[ which(is.na(y.resd)) ] <- 0;
  
  nna.vec <- rep(0, NROW(nna.mat));
  for(i in 1:NCOL(nna.mat))
    nna.vec <- nna.vec*2 + nna.mat[,i];
  
  return(nna.vec);
}


dmvnorm_fast<-function( y.resd, mu, cov.mat, nna.vec=NULL, log=T)
{
  if(is.null(nna.vec))
    nna.vec = get_non_na_number(y.resd);
  
  pv <- rep(NA, NROW(nna.vec));
  for( nna in unique(nna.vec) )
  {
    if(nna>0)
    {
      nna.idx <- which(nna.vec==nna);
      col.idx <- !is.na(y.resd[nna.idx[1],]);
      pv[nna.idx] <- dmvnorm( y.resd[nna.idx, col.idx, drop=F ],
                              mu[col.idx] ,
                              cov.mat[col.idx, col.idx, drop=F], log=log);
    }
  }
  
  return(pv);
}

