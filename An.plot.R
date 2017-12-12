An.plot.mean <- function(dat,par0,tri=2:9){
  
  require(ggplot2)
  
  p1 <- dat$pheno1[,tri]/10
  p2 <- dat$pheno2[,tri]/10
  times <- dat$sample.time[tri]
  
  mpheno1 <- as.numeric(colMeans(p1,na.rm=T))
  mpheno2 <- as.numeric(colMeans(p2,na.rm=T))
  init.par <- par0[7:16]
  tmin <- min(times);tmax <- max(times)
  s.t <- seq(tmin,tmax,0.1)
  
  NY <- com.get_mu(init.par,times=s.t,x0=mpheno1[1],y0=mpheno2[1])
  IND <- com.get_mu.ind(c(init.par[1:3],init.par[6:8]),times=s.t,x0=mpheno1[1],y0=mpheno2[1])
  INTER <- NY-com.get_mu.ind(c(init.par[1:3],init.par[6:8]),times=s.t,x0=mpheno1[1],y0=mpheno2[1])
  
  n <- dim(dat$pheno1)[1]
  tp1 <- c()
  for(i in 1:n){
    
    nd <- as.numeric(p1[i,])
    nd.1 <- cbind(rep(i,length(nd)),times,nd)
    tp1 <- rbind(tp1,nd.1)
  }
  
  tp2 <- c()
  for(i in 1:n){
    
    ndbh <- as.numeric(p2[i,])
    ndbh.1 <- cbind(rep(i,length(ndbh)),times,ndbh)
    tp2 <- rbind(tp2,ndbh.1)
  }
  
  colnames(tp1) <- c("index","time","pheno")
  tp1 <- as.data.frame(tp1)
  colnames(tp2) <- c("index","time","pheno")
  tp2 <- as.data.frame(tp2)
  
  fit.p <- data.frame(index=rep(1,length(s.t)),time=s.t,pheno=NY[1:length(s.t)])
  fit.ind <- data.frame(index=rep(1,length(s.t)),time=s.t,pheno=IND[1:length(s.t)])
  fit.inter <- data.frame(index=rep(1,length(s.t)),time=s.t,pheno=INTER[1:length(s.t)])
  
  g1 <- ggplot(tp1)
  g1 <- g1 + geom_line(aes(x=time,y=pheno,group=index),colour="#666666",size=0.5,alpha=0.5)
  g1 <- g1 + geom_line(data=fit.p,aes(x=time,y=pheno,group=index),colour="red",size=1)
  g1 <- g1 + geom_line(data=fit.ind,aes(x=time,y=pheno,group=index),colour="#4F94CD",size=1,linetype=2)
  g1 <- g1 + geom_line(data=fit.inter,aes(x=time,y=pheno,group=index),colour="#A020F0",size=1,linetype=3)
  g1 <- g1 + scale_x_continuous(limits=c(tmin,tmax),breaks=seq(tmin,tmax,1),labels=seq(tmin,tmax,1))
  g1 <- g1 + geom_hline(yintercept=0,linetype=1,size=0.3)
  g1 <- g1 + scale_y_continuous(limits=c(-10,32))
  g1 <- g1 +  annotate("text", x = 2.5,y =32*0.95, label = "A",size=8)
  g1 <- g1 + xlab("Time (week)")+ylab("Height (cm)") + theme_zg()
  
  fit.pp <- data.frame(index=rep(1,length(s.t)),time=s.t,pheno=NY[(1+length(s.t)):(2*length(s.t))])
  fit.IND1 <- data.frame(index=rep(1,length(s.t)),time=s.t,pheno=IND[(1+length(s.t)):(2*length(s.t))])
  fit.INTER1 <- data.frame(index=rep(1,length(s.t)),time=s.t,pheno=INTER[(1+length(s.t)):(2*length(s.t))])

  g2 <- ggplot(tp2)
  g2 <- g2 + geom_line(aes(x=time,y=pheno,group=index),colour="#666666",size=0.5,alpha=0.5)
  g2 <- g2 + geom_line(data=fit.pp,aes(x=time,y=pheno,group=index),colour="red",size=1)
  g2 <- g2 + geom_line(data=fit.IND1,aes(x=time,y=pheno,group=index),colour="#4F94CD",size=1,linetype=2)
  g2 <- g2 + geom_line(data=fit.INTER1,aes(x=time,y=pheno,group=index),colour="#A020F0",size=1,linetype=3)
  g2 <- g2 + scale_x_continuous(limits=c(tmin,tmax),breaks=seq(tmin,tmax,1),labels=seq(tmin,tmax,1))
  g2 <- g2 + geom_hline(yintercept=0,linetype=1,size=0.3)
  g2 <- g2 + scale_y_continuous(limits=c(-10,32))
  g2 <- g2 +  annotate("text", x = 2.5, y =32*0.95, label = "B",size=8)
  g2 <- g2 + xlab("Time (week)")+ylab("Height (cm)") + theme_zg()
  
  
  pdf("Figure-1.pdf",height=6,width=14)
  multiplot(g1,g2,cols=2)
  dev.off()
  
  
  
  
}

an.genotype.plot <- function(dat,allpar,tri=2:9,s.i=1){
  
  
  y1 <- as.matrix( dat$pheno1[,tri])/10
  y2 <- as.matrix( dat$pheno2[,tri])/10
  
  times <- dat$sample.time[tri]
  geno_table <- dat$snp
  SNP <- geno_table[s.i,]
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
  
  index <- table(SNP1)
  snp.type <- names(index)
  parin <- matrix(allpar[s.i,8:47],nrow=4,byrow=T)
  nt <- seq(min(times),max(times),0.1)
  NY <- matrix(NA,nrow=4,ncol=length(nt)*2)
  IND <- matrix(NA,nrow=4,ncol=length(nt)*2)
  INTER <- matrix(NA,nrow=4,ncol=length(nt)*2)
  SNP.index <- list()
  m1 <- as.numeric(colMeans(y11,na.rm=T))
  m2 <- as.numeric(colMeans(y22,na.rm=T))
  for(i in 1:length(snp.type)){
    
    SNP.n <- which(SNP1==snp.type[i])

    NY[i,] <- com.get_mu(c(parin[i,]),times=nt,x0=m1[1],y0=m2[1])
    IND[i,] <- com.get_mu.ind(c(parin[i,1:3],parin[i,6:8]),times=nt,x0=m1[1],y0=m2[1])
    INTER[i,] <- NY[i,]-com.get_mu.ind(c(parin[i,1:3],parin[i,6:8]),times=nt,x0=m1[1],y0=m2[1])
    SNP.index[[i]] <- SNP.n
  }
  col.i <- c("red","green","blue","purple")
  file <- paste("M-",s.i,".pdf",sep="")
  pdf(file,height=8,width=10)
  par(mfcol=c(2,4),mar=c(0.1,0.1, 0.1, 0.1) ,mai=c(0.1,0.3,0.1,0.1),oma=c(6,6,3,3))
  for(i in 1:length(SNP.index)){
    
    plot(0,0,xlim=c(2,9),ylim=c(min(INTER[,1:length(nt)]),max(dat$pheno1/10,na.rm=T)),
         xlab="",ylab="",cex.axis=2,cex.lab=2,xaxt="n",yaxt="n")
    #axis(2,seq(2,9,2),seq(2,9,2),lwd=2,cex.axis=1.2)
    ind <- SNP.index[[i]]
    for(j in 1:length(SNP.index[[i]])){
      
      lines(2:9,y11[ind[j],],col="#B8B8B8",lwd=1)
    }
    abline(h=0,v = 0, col = "black",lwd=1)
    lines(nt,NY[i,1:length(nt)],col=col.i[i],lwd=2)
    lines(nt,IND[i,1:length(nt)],col=col.i[i],lwd=2,lty=2)
    lines(nt,INTER[i,1:length(nt)],col=col.i[i],lwd=2,lty=3)
    text(2.8,max(dat$pheno1/10,na.rm=T)*0.97,labels=snp.type[i],cex=1.5)
    if(i==1){
      #y.l <- "log(Height of Arabidopsis 1 (mm))"
      mtext("Height of arabidopsis 1 (cm)",side=2, line=3,cex=1.2,adj=0.5)
      axis(2, seq(-10,30,5), seq(-10,30,5),lwd=1,cex.axis=1.5)
    }
    plot(0,0,xlim=c(2,9),ylim=c(min(INTER[,(length(nt)+1):(2*length(nt))]),max(dat$pheno2/10,na.rm=T)),
         xlab=" ",ylab="",cex.axis=2,cex.lab=2,xaxt="n",yaxt="n")
    ind <- SNP.index[[i]]
    for(j in 1:length(SNP.index[[i]])){
      
      lines(2:9,y22[ind[j],],col="#B8B8B8",lwd=1)
    }
    abline(h=0,v = 0, col = "black",lwd=1)
    lines(nt,NY[i,(length(nt)+1):(2*length(nt))],col=col.i[i],lwd=2)
    lines(nt,IND[i,(length(nt)+1):(2*length(nt))],col=col.i[i],lwd=2,lty=2)
    lines(nt,INTER[i,(length(nt)+1):(2*length(nt))],col=col.i[i],lwd=2,lty=3)
    text(2.8,max(dat$pheno2,na.rm=T)*0.97,labels=snp.type[i],cex=1.5)
    if(i==1){
      #y.l <- "log(Height of Arabidopsis 1 (mm))"
      mtext("Height of arabidopsis 2 (cm)",side=2, line=3,cex=1.2,adj=0.5)
      axis(2, seq(-10,30,5), seq(-10,30,5),lwd=1,cex.axis=1.5)
    }
    mtext("Time (week)",side=1, line=3,cex=1.2,adj=0.5)
    axis(1, seq(2,9,1),seq(2,9,1),lwd=1,cex.axis=1.5)
  }
  
  dev.off()
  
}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(base_family="Times") +
    theme(...,rect = element_rect(fill=bg),
          plot.margin = unit(c(0.5,0.5,1,1), 'lines'),
          panel.background = element_rect(fill='transparent', color='black'),
          panel.border = element_rect(fill='transparent', color='transparent'),
          panel.grid = element_blank(),
          axis.title = element_text(color='black',size=rel(3)),
          axis.title.x = element_text(margin = unit(c(0.5, 0, 0, 0), "cm")),
          axis.title.y = element_text(margin = unit(c(0, 0.5, 0, 0), "cm")),
          axis.ticks.length = unit(-.25, "cm"),
          axis.ticks = element_line(color='black',size = 0.7),
          axis.text.x= element_text(color='black',margin = unit(c(0.5, 0, 0, 0), "cm"),size=rel(3)),
          axis.text.y= element_text(color='black',angle = 90,margin = unit(c(0, 0.6, 0, 0), "cm"),size=rel(3),hjust=0.5),
          legend.position = 'none',
          legend.title = element_blank(),
          legend.key = element_rect(fill='transparent', color='transparent'),
          strip.background = element_rect(fill='transparent', color='transparent'),
          strip.text.x = element_blank(),
          text = element_text(debug=FALSE),
          panel.spacing = unit(0.4, "lines"),
          #strip.placement = "top",
          plot.title = element_text(size = rel(5))
    )
}


draw.line <- function (xc, yc, w, l1, l2, col=col, lwd=lwd, lend=1) {
  w  <- (w/360)*2*pi;
  x1 <- xc+l1*cos(w);
  y1 <- yc-l1*sin(w);
  x2 <- xc+l2*cos(w);
  y2 <- yc-l2*sin(w);
  segments(x1, y1, x2, y2, col=col, lwd=lwd, lend=lend);
}



draw.arc.s <- function (xc, yc, r, w1, w2, col="lightblue", lwd=1, lend=1){
  
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 5;
  if (pix.n < 2){
    pix.n <- 2;
  }
  
  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.i.x <- xc + cos(ang.seq) * r;
  fan.i.y <- yc - sin(ang.seq) * r;
  lines(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend);
}

draw.arc.st <- function (xc, yc, r, w1, w2, col="lightblue", lwd=1, lend=1){
  
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 5;
  if (pix.n < 2){
    pix.n <- 2;
  }
  
  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.i.x <- xc + cos(ang.seq) * r;
  fan.i.y <- yc - sin(ang.seq) * r;
  lines(fan.i.x, fan.i.y, col=col, lwd=lwd, lty=2);
}
draw.text.rt <- function(xc, yc, r, w, n, col="black", cex=1, side="out",LG){
  
  w     <- w%%360;
  the.o <- w;
  
  the.w <- 360-w;
  w     <- w/360*2*pi;
  x     <- xc+r*cos(w);
  y     <- yc-r*sin(w);
  
  
  num2  <- 26;
  
  if (side=="out"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w ;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.w <- the.w +180;
      the.pos <- 4;
    }
    
    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  } 
  if(LG==1)
    text(x, y, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  if(LG==2)
    text(x, y, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  
  if(LG==3)
    text(x, y, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  if(LG==4)
    text(x, y, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  if(LG==5)
    text(x, y, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  
  if(LG==6)
    text(x+25, y-8, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  if(LG==7)
    text(x+15, y-20, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  if(LG==8)
    text(x, y-25, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  
  if(LG==9)
    text(x+5.5, y+10, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  
  if(LG==10)
    text(x+20, y+18, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  
  
  if(LG==11)
    text(x+18, y+5, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex);
  
  
}


draw.text.rt.axes <- function(xc, yc, r, w, n, col="black", cex=1, side="out"){
  
  w     <- w%%360;
  the.o <- w;
  
  the.w <- 360-w;
  w     <- w/360*2*pi;
  x     <- xc+r*cos(w);
  y     <- yc-r*sin(w);
  
  
  num2  <- 30;
  
  if (side=="out"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w +180;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.w <- the.w +180+180;
      the.pos <- 4;
    }
    
    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  } 
  
  text(x, y, adj=0, offset=1, labels=n, srt=the.w, 
       pos=the.pos, col=col, cex=cex);
}

do.scale.cir <-function (xc=xc, yc=yc, the.r=the.r, total.num=total.num, 
                         col="blue", lwd=0.001,
                         V1=V1, V2=V2, V3=V3, V4=V4){
  
  sum.po  <- as.numeric(total.num);         # total number
  scale.w <- 10#sum.po/150;                         # one scale size
  scale.l <- nchar(as.integer(scale.w))-1;       # length of the number
  scale.i <- as.integer(scale.w/(10^scale.l));   # the first digital
  scale.d <- scale.i * 10^scale.l;               # the first digital, then 0
  scale.m <- 1 * 10^scale.l;                     # the unit of the scale   
  #draw.arc.s(xc, yc, w1=V1, w2=V2, r=the.r, col=col, lwd=lwd);
  
  start.p  <- 0;
  w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);  
  scale.n  <- as.integer(as.numeric(V4)/scale.d+1);
  
  for (j in 1:scale.n){
    po.s <- as.integer(start.p/scale.m);
    
    draw.line(xc,    yc,   w, the.r-10, the.r+3, col=col, lwd=lwd);
    draw.text.rt.axes(xc, yc,  the.r+6, w, po.s*10, col=1, cex=0.5);
    
    start.p <- start.p + scale.d/2;
    w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);
    
    if (w <= V2){  
      draw.line(xc, yc, w, the.r-10, the.r+1, col=col, lwd=lwd);
    }
    
    start.p <- start.p + scale.d/2;
    w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);
  }
  
}


draw.line3 <- function (xc, yc, w1, w2, r1, r2, col=col, lwd=lwd,lty=1){
  theangle1 <- w1;
  theangle2 <- w2;
  l1        <- r1;
  l2        <- r2;
  
  theangle1 <- (theangle1/360)*2*pi;
  x1        <- xc+l1*cos(theangle1);
  y1        <- yc-l1*sin(theangle1);
  
  theangle2 <- (theangle2/360)*2*pi;
  x2        <- xc+l2*cos(theangle2);
  y2        <- yc-l2*sin(theangle2);
  
  segments(x1, y1, x2, y2, col=col, lwd=lwd,lty=lty,lend="butt");
}

scale.v <- function(v, a, b, min.v, max.v) {
  v <- v-min.v; 
  v <- v/(max.v-min.v); 
  v <- v*(b-a);  
  v+a
}


draw.arc.pg <- function (xc, yc, 
                         w1, w2, r1, r2, col="lightblue", border="lightblue", lwd=0.01
){
  
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.i.x <- xc + cos(ang.seq) * r1;
  fan.i.y <- yc - sin(ang.seq) * r1;
  
  
  fan.o.x <- xc + cos(ang.seq) * r2;
  fan.o.y <- yc - sin(ang.seq) * r2;
  
  polygon(c(rev(fan.i.x), fan.o.x ), c(rev(fan.i.y), fan.o.y), 
          fillOddEven=F, border=border, col=col, lwd=lwd, lend=1)
  
}



draw.line4 <- function (xc, yc, w1, w2, r1, r2, col=col, lwd=lwd,lty=1){
  theangle1 <- w1;
  theangle2 <- w2;
  l1        <- r1;
  l2        <- r2;
  
  theangle1 <- (theangle1/360)*2*pi;
  x1        <- xc+l1*cos(theangle1);
  y1        <- yc-l1*sin(theangle1);
  
  theangle2 <- (theangle2/360)*2*pi;
  x2        <- xc+l2*cos(theangle2);
  y2        <- yc-l2*sin(theangle2);
  
  segments(x1, y1, x2, y2, col=col, lwd=lwd);
}


draw.link <- function(xc, yc, r, w1, w2, col=col, lwd=lwd) {
  # for translocation
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  points(bezierCurve(x,y,60), type="l", col=col, lwd=lwd, lend="butt")
}

draw.link.pg <- function(xc, yc, r, w1.1, w1.2, w2.1, w2.2, col=col, lwd=lwd) {
  #####################################################
  w1 <- w1.1;
  w2 <- w2.2;
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  bc1 <- bezierCurve(x,y,60);
  
  ang.d <- abs(w1.1-w1.2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w1.1,w1.2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.1.x <- xc + cos(ang.seq) * r;
  fan.1.y <- yc - sin(ang.seq) * r;
  
  ######################################################
  w1 <- w1.2;
  w2 <- w2.1;
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  bc2 <- bezierCurve(x,y,60);
  
  ang.d <- abs(w2.1-w2.2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w2.1,w2.2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.2.x <- xc + cos(ang.seq) * r;
  fan.2.y <- yc - sin(ang.seq) * r;
  
  polygon(c(bc1$x, fan.2.x, rev(bc2$x), rev(fan.1.x)), 
          c(bc1$y, fan.2.y, rev(bc2$y), rev(fan.1.y)), 
          fillOddEven=T, border=col, col=col); 
}



bezierCurve <-function(x, y, n=10)  {  
  outx <- NULL	
  outy <- NULL 	
  i <- 1	
  for (t in seq(0, 1, length.out=n))		{		
    b <- bez(x, y, t)		
    outx[i] <- b$x		
    outy[i] <- b$y 		
    i <- i+1		
  } 	
  return (list(x=outx, y=outy))	
}


bez <-function(x, y, t)	{	
  outx <- 0	
  outy <- 0	
  n <- length(x)-1	
  for (i in 0:n)		{		
    outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]		
    outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]		
  } 	
  return (list(x=outx, y=outy))	
}

do.line <- function (xc=xc, yc=yc, the.r=the.r,
                     col="blue", lwd=0.001,
                     V1=V1, V2=V2, V3=V3, V4=V4,res){
  
  for (j in 1:length(res)){
    
    w <- scale.v(as.numeric(res[j]), V1, V2, V3, V4);  
    draw.line(xc,    yc,   w, the.r-30, the.r, col=col, lwd=lwd);
    #draw.arc.pg(xc, yc, w1=w, w2=w, r1=the.r-20, r2=the.r,
    #col="#00000050", border=NA, lwd=0.01
    #)
    
  }
  
}


draw.point.w <- function(xc, yc, r, w, col=col, cex=cex,pch=20){
  w <- w/360*2*pi;
  x <- xc+r*cos(w);
  y <- yc-r*sin(w);
  points(x, y, pch=pch, col=col, cex=cex);
}



ciro.dat.g <- function(dat,LR){
  
  seg.num    <- length(unique(dat$snp.i[,1]))
  seg.po <- c()
  for(i in 1:seg.num){
    
    index <- which(dat$snp.i[,1]==i)
    seg.po <- c(seg.po,dat$snp.i[max(index),2])
  }
  
  gap.angle.size <- 2
  seg.angle.from <- 270
  
  seg.full.l  <- sum(as.numeric(seg.po))
  cir.angle.r <- (360-seg.num*gap.angle.size)/seg.full.l
  
  out.s     <- c()
  l.old     <- 0
  gap.angle <- 0
  
  for (i in 1:seg.num){
    seg.n <- paste("Chr",i,sep="")
    len   <- cumsum(seg.po)[i]
    w1    <- cir.angle.r*l.old + gap.angle
    w2    <- cir.angle.r*len   + gap.angle
    out.s     <- rbind(out.s, c(seg.n, w1+seg.angle.from, w2+seg.angle.from, l.old, len, 0, seg.po[i]))
    gap.angle <- gap.angle + gap.angle.size
    l.old     <- len;
  }
  MAF <- rep(NA,dim(dat$snp)[1])
  for(i in 1:dim(dat$snp)[1]){
    mafs <- table(dat$snp[i,])
    MAF[i] <-(2*mafs[1]+mafs[2]+mafs[3])/(sum(mafs)*2)
  }
  logp <- -log10(pchisq(LR,df=30,lower.tail = F))
  mapping <- cbind(dat$snp.i[,1:2],MAF,logp)
  
  colnames(mapping) <- c("seg.name", "seg.po","MAF","logp")
  
  ciro <- list(chr.po=out.s,mapping=mapping)
  return(ciro)
}

ciro.plot <- function(ciro,xc=400,yc=400,r=350,xlab=NULL,ylab=NULL,
                      col.out=1,print.chr.lab=T,scale=T,filename="1.pdf"){
  
  pdf(filename, 6, 6)
  par(mar=c(0, 0, 0, 0))
  plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab=xlab, ylab=ylab, main="")
  redtrans <- rgb(0, 0, 0, 4, maxColorValue=255) 
  col.out=c("steelblue","salmon", "palevioletred","steelblue","darkgreen")
  chr.po <- ciro$chr.po
  chr.num     <- nrow(chr.po)
  mapping <- ciro$mapping
  for (chr.i in c(1:chr.num)){
    w1 <- as.numeric(chr.po[chr.i,2])
    w2 <- as.numeric(chr.po[chr.i,3])
    draw.arc.s(xc, yc, r+20, w1, w2, col=col.out[chr.i], lwd=7)
    draw.arc.s(xc, yc, r=r-110, w1, w2, col=col.out[chr.i], lwd=5)
    v1 <- as.numeric(chr.po[chr.i,2])
    v2 <- as.numeric(chr.po[chr.i,3])
    v3 <- as.numeric(chr.po[chr.i,6])
    v4 <- as.numeric(chr.po[chr.i,7])
    chr.n <- gsub("Chr","",chr.po[chr.i,1])
    res <- as.numeric((mapping[mapping[,1]==chr.n,2]))
    do.line(xc=xc, yc=yc, the.r=r+6,col=redtrans , lwd=0.01,
            V1=v1, V2=v2, V3=v3, V4=v4,res)
    draw.arc.s(xc, yc, r=r+10, w1, w2, col="grey", lwd=1)
    if (print.chr.lab){
      w.m <- (w1+w2)/2
      r.j <- r/20+30
      chr.t <- chr.po[chr.i,1]
      draw.text.rt(xc, yc, r+r.j, w.m, chr.t, cex=1.0,LG=chr.i)
    } 
    
    draw.arc.s(xc, yc, r=r-118, w1, w2, col="grey", lwd=1)
    draw.arc.s(xc, yc, r=r-240, w1, w2, col="black", lwd=1)
    #draw.line3(xc, yc, w1, w2, r1=r-10, r2=r-10, col="blue", lwd=lwd)
  }
  allcol <- c("#6495ED","#76EE00","#6495ED","#76EE00","#6495ED")
  mapping <- ciro$mapping
  for(col.v in 1:1){
    dat.min <- min(as.numeric(mapping[,col.v+2]), na.rm=T)
    dat.max <- max(as.numeric(mapping[,col.v+2]), na.rm=T)
    for (chr.i in 1:chr.num){
      
      w1 <- as.numeric(chr.po[chr.i,2])
      w2 <- as.numeric(chr.po[chr.i,3])
      chr.s <- chr.po[chr.i,1]   
      chr.s <- gsub("Chr","",chr.s)
      dat   <- subset(mapping, mapping[,1]==chr.s)
      dat   <- dat[order(as.numeric(dat[,2])),]
      
      v1 <- as.numeric(chr.po[chr.i,2])
      v2 <- as.numeric(chr.po[chr.i,3])
      v3 <- as.numeric(chr.po[chr.i,6])
      v4 <- as.numeric(chr.po[chr.i,7])
      
      my.R1 <- r -10
      my.R2 <- r -100
      my.v    <- as.numeric(dat[1, col.v+2])     
      v.old   <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max);
      po      <- as.numeric(dat[1,2]);
      w.from  <- scale.v(po, v1, v2, v3, v4);
      for (i in 2:nrow(dat)){
        my.v  <- as.numeric(dat[i, col.v+2]);
        
        if (is.na(my.v)){
          next;
        }  
        my.R1 <- r -25;
        my.R2 <- r -100;
        
        v     <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max);
        po    <- as.numeric(dat[i,2]);
        w.to  <- scale.v(po, v1, v2, v3, v4);
        if (v.old > 0){            
          draw.line3(xc, yc, w.from, w.to, v.old, v, col=col.out[chr.i], lwd=1)
        } 
        #draw.point.w(xc, yc, vt, w.to, col=coln[col.v], cex=0.2)
        #draw.line4(xc, yc, w.from, w.to, vt.old, vt, col=coln[col.v], lwd=1)
        v.old <- v;
        w.from <- w.to;
      }
    }
  }
  
  for(col.v in 2:2){
    dat.min <- min(as.numeric(mapping[,col.v+2]), na.rm=T)
    dat.max <- max(as.numeric(mapping[,col.v+2]), na.rm=T)
    for (chr.i in 1:chr.num){
      
      w1 <- as.numeric(chr.po[chr.i,2])
      w2 <- as.numeric(chr.po[chr.i,3])
      chr.s <- chr.po[chr.i,1]   
      chr.s <- gsub("Chr","",chr.s)
      dat   <- subset(mapping, mapping[,1]==chr.s)
      dat   <- dat[order(as.numeric(dat[,2])),]
      
      v1 <- as.numeric(chr.po[chr.i,2])
      v2 <- as.numeric(chr.po[chr.i,3])
      v3 <- as.numeric(chr.po[chr.i,6])
      v4 <- as.numeric(chr.po[chr.i,7])
      
      my.thred <- -log10(1e-20/400000)
      my.R1 <- r -118
      my.R2 <- r -208
      my.v    <- as.numeric(dat[1, col.v+2])     
      v.old   <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max);
      vt.old <- scale.v(my.thred, my.R1, my.R2, dat.min, dat.max);
      po      <- as.numeric(dat[1,2]);
      w.from  <- scale.v(po, v1, v2, v3, v4);
      for (i in 2:nrow(dat)){
        my.v  <- as.numeric(dat[i, col.v+2]);
        
        if (is.na(my.v)){
          next;
        }  
        my.R1 <- r -118;
        my.R2 <- r -228;
        
        v     <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max);
        vt     <- scale.v(my.thred, my.R1, my.R2, dat.min, dat.max);
        po    <- as.numeric(dat[i,2]);
        w.to  <- scale.v(po, v1, v2, v3, v4);
        if (v.old > 0){      
          draw.point.w(xc, yc, v, w.to, col=col.out[chr.i], cex=0.3,pch=16)
          #draw.line3(xc, yc, w.from, w.to, v.old, v, col="#76EE00", lwd=1)
        } 
        #draw.point.w(xc, yc, vt, w.to, col=coln[col.v], cex=0.2)
        #draw.line4(xc, yc, w.from, w.to, vt.old, vt, col=coln[col.v], lwd=1)
        v.old <- v;
        w.from <- w.to;
      }
      draw.arc.st(xc, yc, vt, w1, w2, col="red", lwd=1)
    }
  }
  
  redtrans <- rgb(0, 0, 0, 80, maxColorValue=255) 
  draw.arc.pg(xc, yc,492.16, 496.57, r-118, r-240, col=redtrans , border=NA , lwd=0.01)
  
  segments(340, 450, 390, 450, col="black", lwd=1.5)
  text(420,450,"MAF",cex=0.7)
  points(365, 430, col="black", cex=1,pch=16)
  text(425,427,"Height",cex=0.7)
  segments(340, 405, 390, 405, col="red", lwd=1.5,lty=2)
  text(425,405,"Threshold",cex=0.7)
  text(400,385,expression(italic(P)<1e-20),cex=0.7)
  rect(345,350,365,370,col=redtrans,border=NA)
  text(420,360,"QTL regions",cex=0.7)
  dev.off()
  
}

select.m.plot <- function(dat,res1){
  
  ll <- -log10(pchisq(res1[,1],df=30,lower.tail = F))
  thred <- -log10(1e-30/400000)
  index <- which(ll>thred)
  dlr <- ll[index]
  d <-(dat$snp.i[index,2])
  pdf("chr4.pdf",height=5,width=6.5)
  par(mar=c(3.2,4,2.2,1))
  plot(0,0,pch=16,type="n",col="#757575",xlab="SNP Position (Mb)",xlim=c(d[1],d[length(d)]),ylim=c(35,max(dlr)*1.05),
       ylab=expression(-log10(p-value)),cex.lab=1.5,mgp = c(2.3, 1, 0),xaxt="n",yaxt="n")
  axis(1,round(seq(0.002*10^6,0.759*10^6,length=4),3),round(seq(0.002,0.759,length=4),3),lwd=2,cex.axis=1.2)
  axis(2,seq(35,48,4),seq(35,48,4),lwd=2,cex.axis=1.2)
  points(d,dlr,pch=6,col="blue",cex=0.7)
  text(126500+254000,49.8,"Chr4",cex=2)
  dev.off()
}