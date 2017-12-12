An.est1 <- function(dat,interval=c(1,10),tri=2:9){
  
  y1 <- as.matrix( dat$pheno1)/10
  y2 <- as.matrix( dat$pheno2)/10
  
  times <- dat$sample.time[tri]
  geno_table <- dat$snp
  nm <- dat$nm
  n1 <- interval[1]
  if(length(interval)==1)
    n2 <- interval[1]
  else
    n2 <- interval[2]
  if(n2 >=nm)
    n2 <- nm
  res <- matrix(NA,nrow=length(c(n1:n2)),ncol=100)
  for(i in n1:n2){
    SNP <- geno_table[i,]
    NSNP <- as.character(apply(matrix(SNP,1),2,function(x){c(rep(x,10))}))
    missing <- which(is.na(NSNP))
    if ( length(missing) > 0)
    {
      SNP1 <- NSNP[-(missing)]
      y11 <- y1[ -(missing), ]
      y22 <- y2[ -(missing), ]
    }else{
      SNP1 <- NSNP
      y11 <- y1
      y22 <- y2
    }
    
    ndat <- dat
    ndat$pheno1 <- y11*10
    ndat$pheno2 <- y22*10
    ndat$n <- dim(y11)[1]
    h01 <- try(an.com.H0(ndat,tri),TRUE)
    if (class(h01) == "try-error") 
      h01 <- NA
    #h02 <- try(an.com.H1(y11/10,y22/10,SNP1,init.par=h01,times),TRUE)
    ny11 <- y11[,tri];ny22 <- y22[,tri]
    h02 <- try(an.com.H11(y11=ny11,y22=ny22,SNP1,init.par=h01,times),TRUE)
    #if (class(h02) == "try-error"){
    #  h02 <- try(an.com.H111(y11=ny11,y22=ny22,SNP1,init.par=h01,times),TRUE)
    #}
    if (class(h02) == "try-error")
      h02 <- NA

    LR <- 2*(h01[1]-h02[1])
    if(is.na(h01)||is.na(h02)){
      allpar <- c(LR,rep(NA,25))
    }else{
      allpar <- c(LR,h02)
    }
    
    #if(i==309){
    #  allpar[1] <- NA
    #}
    cat("snp", i, "=", allpar, "\n");
    res[(i-(n1-1)),(1:length(allpar))] <- allpar
  }
  return(res)
}


an.com.H1 <- function(y11,y22,SNP1,init.par,times){
  
  index <- table(SNP1)
  snp.type <- names(index)
  
  g.par <- c()
  SNP.index <- list()
  for(i in 1:length(snp.type)){
    
    SNP.n <- which(SNP1==snp.type[i])
    #m1 <- as.numeric(colMeans(y11[SNP.n,],na.rm=T))
    #m2 <- #as.numeric(colMeans(y22[SNP.n,],na.rm=T))
    na.i <- unique(c(which(is.na(m1)),which(is.na(m2))))
    if(length(na.i)>0){
      ntimes <- times[-na.i]
      nm1 <- m1[-na.i]
      nm2 <- m2[-na.i]
    }else{
      ntimes <- times
      nm1 <- m1
      nm2 <- m2
    }
    
    r1 <- optim(init.par[7:16],s.mle,s.y=c(nm1,nm2),s.t=ntimes,x0=nm1[1],y0=nm2[1],
                method="BFGS",control=list(maxit=32000))
    par <- r1$par
    g.par <- c(g.par,par)
    SNP.index[[i]] <- SNP.n
  }
  
  loop_k <- 1;
  max_iter <- 100;
  epsi <- 10^-5;
  max_err <- 1;
  
  #parin <- c(init.par[2:16],init.par[7:16],init.par[7:16],init.par[7:16])
  parin <- c(init.par[2:6],g.par)
  while(loop_k<max_iter && max_err>epsi){
    
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[-(1:5)])
      AA <- mle.fun(nnpar,y11=y11,y22=y22,times=times,SNP.index,snp.type)
      AA
    }
    r1.covar <- optim(parin[1:5],mle.covar1,method = "BFGS",control=list(maxit=2000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle1.g <- function(npar){
      nnpar <- c(new.covar1,npar)
      AA <- mle.fun(nnpar,y11=y11,y22=y22,times=times,SNP.index,snp.type)
      AA
    }
    r1.g <- optim(c(parin[-(1:5)]),mle1.g,method = "BFGS",control=list(maxit=32000))
    #cat("r1.g:",unlist(r1.g$par), "\n");
    
    newpar <- c(new.covar1,r1.g$par)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin <- newpar
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  return(c(r1.g$value,newpar))
}

an.com.H11 <- function(y11,y22,SNP1,init.par,times){
  
  index <- table(SNP1)
  snp.type <- names(index)
  
  #g.par <- c()
  SNP.index <- list()
  for(i in 1:length(snp.type)){
    
    SNP.n <- which(SNP1==snp.type[i])
    #m1 <- as.numeric(colMeans(y11[SNP.n,],na.rm=T))
    #m2 <- as.numeric(colMeans(y22[SNP.n,],na.rm=T))
    #na.i <- unique(c(which(is.na(m1)),which(is.na(m2))))
    #if(length(na.i)>0){
    #  ntimes <- times[-na.i]
    #  nm1 <- m1[-na.i]
    #  nm2 <- m2[-na.i]
    #}else{
    #  ntimes <- times
    #  nm1 <- m1
    #  nm2 <- m2
   # }
    
    #r1 <- optim(init.par[7:16],s.mle,s.y=c(nm1,nm2),s.t=ntimes,x0=nm1[1],y0=nm2[1],
               # method="BFGS",control=list(maxit=32000))
    #par <- r1$par
    #g.par <- c(g.par,par)
    SNP.index[[i]] <- SNP.n
  }
  
  loop_k <- 1;
  max_iter <- 100;
  epsi <- 10^-5;
  max_err <- 1;
  
  parin <- c(init.par[2:16],init.par[7:16],init.par[7:16],init.par[7:16])
  #parin <- c(init.par[2:6],g.par)
  while(loop_k<max_iter && max_err>epsi){
    
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[-(1:5)])
      AA <- mle.fun(nnpar,y11=y11,y22=y22,times=times,SNP.index,snp.type)
      AA
    }
    r1.covar <- optim(parin[1:5],mle.covar1,method = "BFGS",control=list(maxit=2000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle1.g <- function(npar){
      nnpar <- c(new.covar1,npar)
      AA <- mle.fun(nnpar,y11=y11,y22=y22,times=times,SNP.index,snp.type)
      AA
    }
    r1.g <- optim(c(parin[-(1:5)]),mle1.g,method = "BFGS",control=list(maxit=32000))
    #cat("r1.g:",unlist(r1.g$par), "\n");
    
    newpar <- c(new.covar1,r1.g$par)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin <- newpar
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  return(c(r1.g$value,newpar))
}
mle.fun <- function(par,y11=y11,y22=y22,times=times,
                    SNP.index=SNP.index,snp.type=snp.type){
  
  Y1 <- cbind(y11,y22)
  len.cov <- 5
  par.covar <- par[1:len.cov]
  
  cov.mat <- SAD3.get_mat(par.covar,times, 2)

  M1 <- as.numeric(colMeans(y11,na.rm=T))
  M2 <- as.numeric(colMeans(y22,na.rm=T))
  len.gen <- 10
  len <- 0
  A1 <- c()
  for(i in 1:length(snp.type)){
    mu.g <- par[(len.cov+len+1):(len.cov+len+len.gen)]
    #M1 <- as.numeric(colMeans(y11[SNP.index[[i]],],na.rm=T))
    #M2 <- as.numeric(colMeans(y22[SNP.index[[i]],],na.rm=T))
    mu <- com.get_mu(mu.g,times,x0=M1[1],y0=M2[1])
    ns <- length(SNP.index[[i]])
    if (any(is.na(mu) ))
      return (NaN);
    y <- Y1[SNP.index[[i]],]
    Y.delt <- y- matrix(rep(mu,ns),nrow=ns,byrow=T)
    
    nna.vec <- get_non_na_number(Y.delt)
    
    pv <- dmvnorm_fast( Y.delt, rep(0, NCOL(cov.mat)), cov.mat, nna.vec=nna.vec, log=T)
    
    A1 <- c(A1,-sum( pv ))
    
    #M1 <- as.numeric(colMeans(y11[SNP.index[[i]],],na.rm=T))
    #M2 <- as.numeric(colMeans(y22[SNP.index[[i]],],na.rm=T))
    #mu <- com.get_mu(mu.g,times,x0=0.663,y0=0.633)
    
    #for ( ii in 1:dim(yy1)[2] )
    #{
    #  yy1.miss <- which( is.na(yy1[,ii]));
    #  nyy1[yy1.miss,ii] <- mu[ii];
    #}
    
    #fy1 <- dmvnorm( nyy1, mu, sigma)
    #fy1[which(fy1<=.Machine$double.eps)] <- .Machine$double.eps
    len <- len + len.gen
  }
  A <- sum(A1)
  if ((is.na(A) ))
    return (NaN);
  #cat("LL=",A,"\n")
  return (A);
}

