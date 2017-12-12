
an.load <- function(phenotype="Height-pairs.txt",
                    genotype="../data/final-snp-n.txt",
                    pheno.info="Pairs.txt",
                    genotype.info="../data/snp-info.txt"){
  
  r.pheno <- read.table(phenotype,header=T)
  r.pheno.i <- read.table(pheno.info,header=T)
  r.snp <- read.table(genotype,header=T)
  r.snp.i <- read.table(genotype.info,header=T)
  r.name <- colnames(r.snp)
  r.n1 <- c()
  for(i in 1:length(r.name)){
    r.n1 <- c(r.n1,as.numeric(strsplit(r.name[i],"X")[[1]][2]))
  }
  
  all.p <- c()
  for(j in 1:dim(r.pheno.i)[1]){
    
    if(any(r.pheno.i[j,1]==r.n1)&&any(r.pheno.i[j,2]==r.n1))
      all.p <- c(all.p,j)
  }
  
  all.pp <- unique(as.numeric(unlist(c(r.pheno.i[-all.p,]))))
  
  pname <- rownames(r.pheno)
  p.i <- c()
  for(i in 1:length(pname)){
    
    p.i <- c(p.i,strsplit(pname[i],split="[.]")[[1]][2])
  }
  
  del.i <- c()
  for(i in 1:length(all.pp)){
    del.i <- c(del.i,which(all.pp[i]==p.i))
  }
  
  pheno <- r.pheno[-del.i,]
  pheno.i <- r.pheno.i[all.p,]
  s1 <- seq(1,dim(pheno)[1],20)
  s2 <- seq(11,dim(pheno)[1],20)
  
  ss1 <- c();ss2 <- c()
  for(i in 1:length(s1)){
    ss1 <- c(ss1,s1[i]:(s1[i]+9))
    ss2 <- c(ss2,s2[i]:(s2[i]+9))
  }
  
  pheno.n <- (p.i[-del.i])
  
  pheno.n1 <- pheno.n[ss1]
  pheno.n2 <- pheno.n[ss1]
  pheno1 <- pheno[ss1,]
  pheno2 <- pheno[ss2,]
  
  p1 <- c()
  for( p in 1:length(pheno.i[,1])){
    p1 <- c(p1,which(pheno.i[p,1]==r.n1))
  }
  p2 <- c()
  for( p in 1:length(pheno.i[,2])){
    p2 <- c(p2,which(pheno.i[p,2]==r.n1))
  }
  
  snp1 <- r.snp[,unique(p1)]
  snp2 <- r.snp[,unique(p2)]
  snp.ct <- c()
  for(k in 1:length(unique(p1))){
    tmp.c <- (paste(snp1[,k],snp2[,k],sep=""))
    snp.ct <- cbind(snp.ct,tmp.c)
  }
  rownames(snp.ct) <- rownames(snp1)
  colnames(snp.ct) <- paste("P",1:length(unique(p1)),sep="")
  
  for(i in 1:length(unique(p1))){
    snp.ct[which(snp.ct[,i]=="99"),i] <- NA
    snp.ct[which(snp.ct[,i]=="19"),i] <- NA
    snp.ct[which(snp.ct[,i]=="09"),i] <- NA
    snp.ct[which(snp.ct[,i]=="91"),i] <- NA
    snp.ct[which(snp.ct[,i]=="90"),i] <- NA
  }
  del.snp <- c()
  for(kk in 1:dim(snp.ct)[1])
  if(any(table(snp.ct[kk,])<5))
    del.snp <- c(del.snp,kk)
  
  snp <- snp.ct[-del.snp,]
  snp.i <- r.snp.i[-del.snp,]
  s.time <- 1:dim(pheno1)[2]
  nm <- dim(snp)[1]
  n <- dim(snp)[2]
  res <- list(pheno1=pheno1,pheno2=pheno2,pheno1L=log(pheno1),pheno2L=log(pheno2),
              pheno.i=pheno.i,snp=snp,snp.i=snp.i,
              sample.time=s.time,nm=nm,n=n,pheno.n1=pheno.n1,pheno.n2=pheno.n2)
  res
}