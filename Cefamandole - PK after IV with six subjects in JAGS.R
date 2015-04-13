install.packages(c("MEMSS", "R2jags"))
library(MEMSS)
library(ggplot2)
library(R2jags)
png('Cefamandole.png')
qplot(y=conc,x=Time,col=Subject,data=Cefamandole) +
    scale_y_log10('cefamandole (mcg/ml)')+
    geom_line()+
    theme(legend.position='bottom')
dev.off()

library(R2jags)
datain <-  list(
    time=Cefamandole$Time,
    conc=Cefamandole$conc,
    n=nrow(Cefamandole),
    subject =c(1:nlevels(Cefamandole$Subject))[Cefamandole$Subject],
    nsubject = nlevels(Cefamandole$Subject)
)

model1 <- function() {
    tau <- 1/pow(sigma,2)
    sigma ~ dunif(0,100)
    for (i in 1:n) {
        pred[i] <- (ke1[subject[i]]<ke2[subject[i]])*
            (c1star[subject[i]]*exp(-ke1[subject[i]]*time[i])+
                 c2star[subject[i]]*exp(-ke2[subject[i]]*time[i]))
        conc[i] ~ dnorm(pred[i],tau)
    }
    
    for(i in 1:nsubject) {
        ke1[i] ~ dlnorm(ke[1],kemtau[1])
        ke2[i] ~ dlnorm(ke[2],kemtau[2]) 
        c1star[i] ~ dlnorm(cm[1],ctau[1])
        c2star[i] ~ dlnorm(cm[2],ctau[2])
    }
    for (i in 1:2) {
        kem[i] ~ dunif(-10,10) 
        kemtau[i] <- 1/pow(kesd[i],2) 
        kesd[i] ~ dunif(0,10)
        Ke[i] <- exp(ke[i])
        cm[i] ~ dunif(-10,20)
        ctau[i] <- 1/pow(csd[i],2)
        csd[i] ~ dunif(0,100)
        C[i] <- exp(cm[i])    
    }
    ke <- sort(kem)
} 

parameters <- c('Ke','ke1','ke2','c1star','c2star','C')
inits <- function() {
    list(ke1 = rnorm(6,0.13,.03),
         ke2=  rnorm(6,0.0180,.005), 
         kem=  log(c(rnorm(1,0.13,.03),rnorm(1,0.0180,.005))),
         kesd =runif(2,.1,.4),
         cm = log(rnorm(2,25,5)),
         c1star=rnorm(6,25,3),
         c2star=rnorm(6,25,3)
    )}
jagsfit1 <- jags(datain, model=model1, 
                 inits=inits,
                 parameters=parameters,progress.bar="gui",
                 n.iter=10000,
                 n.chains=4,n.thin=10)
jagsfit1
plot(jagsfit1)
traceplot(jagsfit1)

Time <- seq(0,max(Cefamandole$Time))
la1 <- sapply(1:nrow(jagsfit1$BUGSoutput$sims.matrix),
              function(i) {
                  C1 <- jagsfit1$BUGSoutput$sims.matrix[i,'C[1]']
                  C2 <- jagsfit1$BUGSoutput$sims.matrix[i,'C[2]']
                  k1 <- jagsfit1$BUGSoutput$sims.matrix[i,'Ke[1]']
                  k2 <- jagsfit1$BUGSoutput$sims.matrix[i,'Ke[2]']
                  y=C1*exp(-k1*Time)+C2*exp(-k2*Time)
              })
res1 <- data.frame(
    Time=Time,
    Conc=apply(la1,1,mean),
    C05=apply(la1,1,FUN=function(x) quantile(x,0.05)),
    C95=apply(la1,1,FUN=function(x) quantile(x,0.95)))
png('cef1.png')
ggplot(res1, aes(x=Time)) +
    geom_line(aes(y=Conc))+
    scale_y_log10('cefamandole (mcg/ml)')+
    geom_line(aes(y=conc,x=Time,col=Subject),alpha=1,data=Cefamandole)+
    geom_ribbon(aes(ymin=C05, ymax=C95),alpha=.2)+
    theme(legend.position='bottom')
dev.off()
#########
#
# model 2: proportional error
#
#########
model2 <- function() {
    tau <- 1/pow(sigma,2)
    sigma ~ dunif(0,100)
    for (i in 1:n) {
        pred[i] <- log(
            (ke1[subject[i]]<ke2[subject[i]])*
                (c1star[subject[i]]*exp(-ke1[subject[i]]*time[i])+
                     c2star[subject[i]]*exp(-ke2[subject[i]]*time[i]))
            +.001*(ke1[subject[i]]>ke2[subject[i]]))
        conc[i] ~ dlnorm(pred[i],tau)
    }
    
    
    for(i in 1:nsubject) {
        ke1[i] ~ dlnorm(ke[1],kemtau[1])
        ke2[i] ~ dlnorm(ke[2],kemtau[2]) 
        c1star[i] ~ dlnorm(cm[1],ctau[1])
        c2star[i] ~ dlnorm(cm[2],ctau[2])
    }
    for (i in 1:2) {
        kem[i] ~ dunif(-10,10) 
        kemtau[i] <- 1/pow(kesd[i],2) 
        kesd[i] ~ dunif(0,10)
        Ke[i] <- exp(ke[i])
        cm[i] ~ dunif(-10,20)
        ctau[i] <- 1/pow(csd[i],2)
        csd[i] ~ dunif(0,100)
        C[i] <- exp(cm[i])    
    }
    ke <- sort(kem)
    
} 

parameters <- c('Ke','ke1','ke2','c1star','c2star','C')
inits <- function() {
    list(ke1 = rnorm(6,0.13,.03),
         ke2=  rnorm(6,0.0180,.005), 
         kem=  log(c(rnorm(1,0.13,.03),rnorm(1,0.0180,.005))),
         kesd =runif(2,.1,.4),
         cm = log(rnorm(2,25,5)),
         c1star=rnorm(6,25,3),
         c2star=rnorm(6,25,3)
    )}
jagsfit2 <- jags(datain, model=model2, 
                 inits=inits,
                 parameters=parameters,progress.bar="gui",
                 n.iter=10000,
                 n.chains=4,n.thin=10)
jagsfit2
la2 <- sapply(1:nrow(jagsfit2$BUGSoutput$sims.matrix),
              function(i) {
                  C1 <- jagsfit2$BUGSoutput$sims.matrix[i,'C[1]']
                  C2 <- jagsfit2$BUGSoutput$sims.matrix[i,'C[2]']
                  k1 <- jagsfit2$BUGSoutput$sims.matrix[i,'Ke[1]']
                  k2 <- jagsfit2$BUGSoutput$sims.matrix[i,'Ke[2]']
                  y=C1*exp(-k1*Time)+C2*exp(-k2*Time)
              })


res2 <- data.frame(
    Time=Time,
    Conc=apply(la2,1,mean),
    C05=apply(la2,1,FUN=function(x) quantile(x,0.05)),
    C95=apply(la2,1,FUN=function(x) quantile(x,0.95)))

png('cef2.png')
ggplot(res2, aes(x=Time)) +
    geom_line(aes(y=Conc))+
    scale_y_log10('cefamandole (mcg/ml)')+
    geom_line(aes(y=conc,x=Time,col=Subject),alpha=1,data=Cefamandole)+
    geom_ribbon(aes(ymin=C05, ymax=C95),alpha=.2)+
    theme(legend.position='bottom')
dev.off()
