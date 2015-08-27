simulate_beta_data = function(){
##########################
##########################
#### Beta regression
sim.seed <- 1
time.n <- 20
subject.n <- 100
s2 <- 0.05
beta <- as.matrix(c(-0.5,0.5,0.5))
######
set.seed(sim.seed+10)
X <- as.matrix(data.frame(log.Time=as.matrix(log(rep(seq(1,time.n),subject.n))),
                          Treatment=as.matrix(c(rep(0,subject.n*time.n/2),rep(1,subject.n*time.n/2)))))
#X <- as.matrix(runif(subject.n*time.n,-1,1))
#X  <- as.matrix(c(rep(0,N*time.n/2),rep(1,N*time.n/2)))
Z <- X
set.seed(sim.seed+1)
c <- as.matrix(rnorm(subject.n,mean=0,sd=s2))
c.rep <- as.matrix(as.vector(matrix(c,nrow=time.n,ncol=length(c),byrow=TRUE)))
#####
subject.ind <- as.vector(matrix(paste('Subject_',seq(1,subject.n),sep=''),nrow=time.n,ncol=subject.n,byrow=TRUE))
time.ind  <- rep(seq(1,time.n),subject.n)
######
Z.aug <- cbind(intersept=1,Z)
#browser()
logit.u  <- Z.aug %*% beta + c.rep
######
u <- 1 / (1 + exp(-logit.u))
## v is the phi
v <- 2
######
set.seed(sim.seed+4)
Y <- rbeta(subject.n*time.n, shape1 = u*v, shape2=(1-u)*v)
if(any(Y>1-10^(-6))){Y[Y>1-10^(-6)] <- 1-10^(-6)}

#### For test purpose
#### betareg can not fit random effect model
#### so set the s2 to a small value (small random effect)
#library(betareg)
#tdata <- data.frame(Y=Y,Z,SID=subject.ind)
#gy <- betareg(Y ~ log.Time + as.factor(Treatment) , data = tdata,type='ML')
#summary(gy)
#gy <- betareg(Y ~ log.Time + as.factor(Treatment) , data = tdata,type='BC')
#summary(gy)
#gy <- betareg(Y ~ log.Time + as.factor(Treatment) , data = tdata,type='BR')
#summary(gy)
#
#
#fit_beta_random_effect(Z=Z,Y=Y,subject.ind=subject.ind,time.ind=time.ind,
#                      quad.n=30,verbose=FALSE)
}
