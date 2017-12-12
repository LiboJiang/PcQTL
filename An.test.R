setwd("/home/lbjiang/system mapping/WLN-NEW")
library(mvtnorm)
source("An.load.R")
source("An.ODE.R")
source("an.covar.R")
source("An.curve.R")
source("An.est.R")
source("An.plot.R")


dat <- an.load(phenotype="Height-pairs.txt",genotype="final-snp-example.txt",
               pheno.info="Pairs.txt",genotype.info="snp-info-example.txt")

H0par <- an.com.H0(dat,tri=2:9)
An.plot.mean(dat,par0=H0par,tri=2:9)


res <- An.est1(dat,interval=c(1,10),tri=2:9)


