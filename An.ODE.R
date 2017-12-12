require("deSolve")


get_con_param<-function(parm.id)
{
        for (e in commandArgs())
        {
                ta = strsplit(e,"=", fixed=TRUE);
                if(! is.na( ta[[1]][2]))
                {
                        temp = ta[[1]][2];
                        if( ta[[1]][1] == parm.id) {
                                return (as.character(temp));
                        }
                }
        }

        return(NA);
}













s.mle <- function(s.par,s.y,s.t,x0,y0){
  A <- sum((s.y - com.get_mu(s.par,s.t,x0,y0))^2 )
  A
}

com.get_mu <- function(par, times, x0,y0)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      a1 = par[1],
      k1 = par[2],
      r1 = par[3],
      beta12 = par[4],
      s12 = par[5],
      a2 = par[6],
      k2 = par[7],
      r2 = par[8],
      beta21 = par[9],
      s21 = par[10]);
  }
  
  state0 <- c(X=x0, Y=y0);
  y <- COMP.f( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


COMP.f <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- a1*(1-(X/k1)^r1) + beta12*Y^s12
            dY <- a2*(1-(Y/k2)^r2) + beta21*X^s21
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}

com.get_mu.ind <- function(par, times, x0,y0)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      a1 = par[1],
      k1 = par[2],
      r1 = par[3],
      a2 = par[4],
      k2 = par[5],
      r2 = par[6]);
  }
  
  state0 <- c(X=x0, Y=y0);
  y <- COMP.f.ind( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}




COMP.f.ind <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- a1*(1-(X/k1)^r1) 
            dY <- a2*(1-(Y/k2)^r2) 
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}



Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL )
{
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1] + ti*mu[2];
  if (np.order>=2)
    L <- L + 0.5*(3*ti*ti-1)* mu[3] ;
  if (np.order>=3)
    L <- L + 0.5*(5*ti^3-3*ti)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125*(35*ti^4-30*ti^2+3)* mu[5];
  if (np.order>=5)
    L <- L + 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6];
  if (np.order>=6)
    L <- L + (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)* mu[7];
  if (np.order>=7)
    L <- L + (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)* mu[8];
  if (np.order>=8)
    L <- L + (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)* mu[9];
  if (np.order>=9)
    L <- L + (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)* mu[10];
  if (np.order>=10)
    L <- L + (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)* mu[11];
  if (np.order>=11)
  {
    for(r in 11:(np.order))
    {
      kk <- ifelse(r%%2==0, r/2, (r-1)/2);
      for (k in c(0:kk) )
      {
        L <- L + (-1)^k*factorial(2*r-2*k)/factorial(k)/factorial(r-k)/factorial(r-2*k)/(2^r)*ti^(r-2*k)*mu[r+1];
      }
    }
  }
  return(L);
}



ODEABC.sum <- function(para1,A,B,D){
  
  
  sumA <- sum((D-ODEABC(para1,A,B))^2)
  sumA
}





ODEABC <- function(para1,A,B){
  
  para1[1]*(1-A/ para1[2])^ para1[3] + para1[4]*B^para1[5]
  
}




ode.sovle1 <- function(para1,para2,times,inter1,inter2,nstep=100){
  
  stp <- (max(times)-min(times))/nstep
  nt1 <- seq(min(times),max(times),stp)
  LX1 <- predict(inter1,nt1)$y
  nt2 <- nt1 + stp/2
  LX2 <- predict(inter1,nt2)$y
  nt3 <- nt2 + stp/2
  LX3 <- predict(inter1,nt3)$y
  
  LY1 <- predict(inter2,nt1)$y
  LY2 <- predict(inter2,nt2)$y
  LY3 <- predict(inter2,nt3)$y
  
  NG <- matrix(NA,nrow=nstep+1,ncol=2)
  NG[1,] <- c(LX1[1],LY1[1])
  for(j in 1:nstep){
    
    tg1 <- ODE.G(para1[1:3],LX1[j])
    tg2 <- ODE.G(para1[1:3],LX2[j])
    tg3 <- ODE.G(para1[1:3],LX2[j])
    tg4 <- ODE.G(para1[1:3],LX3[j])
    NG[j+1,1] <- NG[j,1]+stp*(tg1+2*tg2+2*tg3+tg4)/6 
    
    tg1 <- ODE.G(para2[1:3],LY1[j])
    tg2 <- ODE.G(para2[1:3],LY2[j])
    tg3 <- ODE.G(para2[1:3],LY2[j])
    tg4 <- ODE.G(para2[1:3],LY3[j])
    NG[j+1,2] <- NG[j,2]+stp*(tg1+2*tg2+2*tg3+tg4)/6 
  
  }
  list(NG=NG,nt=nt1)
}


ode.sovle <- function(para1,para2,times,inter1,inter2,nstep=100){
  
  stp <- (max(times)-min(times))/nstep
  nt1 <- seq(min(times),max(times),stp)
  LX1 <- predict(inter1,nt1)$y
  nt2 <- nt1 + stp/2
  LX2 <- predict(inter1,nt2)$y
  nt3 <- nt2 + stp/2
  LX3 <- predict(inter1,nt3)$y
  
  LY1 <- predict(inter2,nt1)$y
  LY2 <- predict(inter2,nt2)$y
  LY3 <- predict(inter2,nt3)$y
  
  NG <- matrix(NA,nrow=nstep+1,ncol=2)
  NG[1,] <- c(LX1[1],LY1[1])
  for(j in 1:nstep){
    
    tg1 <- ODEABC(para1,LX1[j],LY1[j])
    tg2 <- ODEABC(para1,LX2[j],LY2[j])
    tg3 <- ODEABC(para1,LX2[j],LY2[j])
    tg4 <- ODEABC(para1,LX3[j],LY3[j])
    NG[j+1,1] <- NG[j,1]+stp*(tg1+2*tg2+2*tg3+tg4)/6 
    
    #tg1 <- ODEABC(para2,LY1[j],LX1[j],LZ1[j])
    #tg2 <- ODEABC(para2,LY2[j],LX2[j],LZ2[j])
    #tg3 <- ODEABC(para2,LY2[j],LX2[j],LZ2[j])
    #tg4 <- ODEABC(para2,LY3[j],LX3[j],LZ3[j])
    #NG[j+1,2] <- NG[j,2]+stp*(tg1+2*tg2+2*tg3+tg4)/6 
   
  }
  list(NG=NG,nt=nt1)
}



#parin <- c(0.2,2,0.3,1,0,1)
ODE.G <- function(parin,X){
  parin[1]*(1-X/ parin[2])^ parin[3]
}




smL <- function(times,para,DS1){
  
  sum((DS1-Legendre.model(t=times,mu=para))^2)
}