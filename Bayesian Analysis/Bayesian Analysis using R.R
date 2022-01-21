
library(rjags)
setwd("/Users/agboo/OneDrive - University of Northern Colorado/Fall 2019/SRM 636/Assignments/TakeHome2")
########### 1a
# The model specification
model_string <- "model {
for (i in 1:n){
y[i]~dnorm(mu[i],tau)
mu[i]<-b0+b[1]*x1[i]+b[2]*x2[i]+b[3]*x3[i]+b[4]*x4[i]+b[5]*x5[i]+b[6]*x6[i]+b[7]*x7[i]+b[8]*x8[i]+b[9]*x9[i]

x4[i]<-x1[i]*x2[i]
x5[i]<-x1[i]*x3[i]
x6[i]<-x2[i]*x3[i]
x7[i]<-pow(x1[i],2)
x8[i]<-pow(x2[i],2)
x9[i]<-pow(x3[i],2)
}

# Priors
b0~dunif(-10000, 10000)
 

for( j in 1:k){
b[j]~dunif(-10000, 10000)

}

tau~dgamma(0.01,0.01)
sigma2<-1/tau
}"

# Running the model
model <- jags.model(textConnection(model_string), data = list(n=26,k=9, y = c(0.22,0.40,0.42,0.44,0.43,0.47,0.44,0.38,0.49,0.46,0.45,0.11,0.43,0.10,0.23,0.31,0.09,0.12,0.08,0.44,0.09,0.12,0.07,0.04,0.25,0.00),
                                                              x1=c(7.30,8.70,8.80,8.10,9.00,8.70,9.30,7.60,10.00,8.40,9.30,7.70,9.80,7.30,8.50,9.50,7.40,7.80,7.70,10.30,7.80,7.10,7.70,7.40,7.30,7.60) ,
                                                              x2 = c(0.00,0.00,0.70,4.00,0.50,1.50,2.10,5.10,0.00,3.70,3.60,2.80,4.20,2.50,2.00,2.50,2.80,2.80,3.00,1.70,3.30,3.90,4.30,6.00,2.00,7.80),
                                                              x3 = c(0.00,0.30,1.00,0.20,1.00,2.80,1.00,3.40,0.30,4.10,2.00,7.10,2.00,6.80,6.60,5.00,7.80,7.70,8.00,4.20,8.50,6.60,9.50,10.90,5.20,20.70)), n.chains = 5, n.adapt= 10000)
update(model, 10000); # Burning for 10000 samples
mcmc_samples <- coda.samples(model, variable.names=c("b","b0"), n.iter=2000000, thin = 1000)


#### PLOT POSTERIOR STATISTICS
plot(mcmc_samples)


#### ERGODIC MEAN
dev.new(width=16,height=12);par(mfrow=c(2,2))
M<-as.matrix(mcmc_samples)
plot(cumsum(M[,10])/cumsum(rep(1,NROW(M[,1]))),type = "l",xlab = "Iterations",ylab = "mu",main = "b11")


#### SHOW POSTERIOR STATISTICS
summary(mcmc_samples)

out <- capture.output(summary(mcmc_samples))
cat("THstat1a", out, file="stats1a.txt", sep="\n", append=TRUE)

#pdf(file="stat1a.pdf")
#plot(...)
#dev.off


#### HPD CREDIBLE INTERVAL
library(HDInterval)
out <- capture.output(hdi(mcmc_samples))
cat("HDI 1a", out, file="stats1a.txt", sep="\n", append=TRUE)



######### 1b
model_string1b <- "model {
  for (i in 1:n){
    y[i]~dnorm(mu[i],tau)
    mu[i]<-b0+b[1]*x1[i]+b[2]*x2[i]+b[3]*x3[i]+b[4]*x4[i]+b[5]*x5[i]+b[6]*x6[i]+b[7]*x7[i]+b[8]*x8[i]+b[9]*x9[i]
    
    x4[i]<-pow(x1[i],2)
    x5[i]<-pow(x2[i],2)
    x6[i]<-pow(x3[i],2)
    x7[i]<-x1[i]*x2[i]
    x8[i]<-x1[i]*x3[i]
    x9[i]<-x3[i]*x2[i]
  }
  
  # Priors
  b0~dunif(-10000, 10000)
  
  for( j in 1:k){
    b[j]~dt(0,1,6)
  }
  
  tau~dgamma(0.1,0.1)
  sigma2<-1/tau
}"

# Running the model
model1b <- jags.model(textConnection(model_string1b), data = list(n=26,k=9, y = c(0.22,0.40,0.42,0.44,0.43,0.47,0.44,0.38,0.49,0.46,0.45,0.11,0.43,0.10,0.23,0.31,0.09,0.12,0.08,0.44,0.09,0.12,0.07,0.04,0.25,0.00),
                                                              x1=c(7.30,8.70,8.80,8.10,9.00,8.70,9.30,7.60,10.00,8.40,9.30,7.70,9.80,7.30,8.50,9.50,7.40,7.80,7.70,10.30,7.80,7.10,7.70,7.40,7.30,7.60) ,
                                                              x2 = c(0.00,0.00,0.70,4.00,0.50,1.50,2.10,5.10,0.00,3.70,3.60,2.80,4.20,2.50,2.00,2.50,2.80,2.80,3.00,1.70,3.30,3.90,4.30,6.00,2.00,7.80),
                                                              x3 = c(0.00,0.30,1.00,0.20,1.00,2.80,1.00,3.40,0.30,4.10,2.00,7.10,2.00,6.80,6.60,5.00,7.80,7.70,8.00,4.20,8.50,6.60,9.50,10.90,5.20,20.70)), n.chains = 3, n.adapt= 10000)
update(model, 1000); # Burning for 1000 samples
mcmc_samples1b <- coda.samples(model1b, variable.names=c("b","b0"), n.iter=3000000, thin = 1000)



#### PLOT POSTERIOR STATISTICS
plot(mcmc_samples1b)


#### ERGODIC MEAN
dev.new(width=16,height=12);par(mfrow=c(2,2))
M1b<-as.matrix(mcmc_samples1b)
plot(cumsum(M1b[,9])/cumsum(rep(1,NROW(M1b[,9]))),type = "l",xlab = "Iterations",ylab = "mu",main = "b9")


#### SHOW POSTERIOR STATISTICS
summary(mcmc_samples1b)

out <- capture.output(summary(mcmc_samples1b))
cat("THstat1b", out, file="stats1b.txt", sep="\n", append=TRUE)

#pdf(file="stat1a.pdf")
#plot(...)
#dev.off


#### HPD CREDIBLE INTERVAL
library(HDInterval)
out <- capture.output(hdi(mcmc_samples1b))
cat("HDI 1b", out, file="stats1b.txt", sep="\n", append=TRUE)

##### Question 5
setwd("/Users/agboo/OneDrive - University of Northern Colorado/Fall 2019/SRM 636/Assignments/TakeHome2")
data<-read.csv("Uneployment54-2018.csv", header = TRUE)
use<-cbind(data[,3],data[,6],data[,9],data[,12])

# Model A: Shifts Transformed Series
model_string5 <- "model{
  for (t in 1:T+1) {u[t] <- mean(U[C[t]+1:C[t+1]])} # calculate six-monthly averages
  for (t in 1:T) {y[t] <- 100*log(1+u[t+1]/100)-100*log(1+u[t]/100)
  LL[t] <- 0.5*log(tau[t]/6.28)-0.5*tau[t]*pow(y[t]-mu[t],2)
  G[t] <- 1/exp(LL[t])}
  y[1] ~ dnorm(mu[1],tau.first); 
  mu[1] <- beta0
  for (t in 2:T) {y[t] ~ dnorm(mu[t],tau[t])
    mu[t] <- gamma*y[t-1]+beta0*(1-gamma)+d1[t]*nu[t]
    
    # shift in means mechanism
    d1[t] ~ dbern(eta1)
    nu[t] ~ dnorm(0,tau.nu)
    
    # shift in variance mechanism
    tau[t] <- tau[t-1]*pow(omega[t],d2[t])
    d2[t] ~ dbern(eta2)
    omega[t] ~ dexp(1)}
    
  # priors
  tau.nu <- 0.1    gamma ~ dnorm(0,0.1);  beta0 ~ dnorm(0,1)
  tau[1] <- tau.1;   tau.first <- (1-gamma*gamma)*tau.1;
  tau.1 ~ dgamma(1,0.001)
  # update density of shift probabilities
  s1 <- sum(d1[2:T])+1;    s2 <- (T-1)-s1+19
  t1 <- sum(d2[2:T])+1;    t2 <- (T-1)-t1+19
  eta1 ~ dbeta(s1,s2);      eta2 ~ dbeta(t1,t2)
  }"


library(coda)
library(rjags)
library(R2jags)

model_string5 <-"model  { 

mu[1] ~ dnorm(0,1)
eps[1] ~ dnorm(0,1)
V[1] ~ dgamma(1,0.001) 

for (t in 1:T) { y[t] ~ dnorm(mu[t]+eps[t],V[t])
}      
      
for (t in 2:T) {
                d1[t] ~ dbern(eta1)
                d2[t] ~ dbern(eta2)
                nu[t] ~ dnorm(0,tau.nu)
                u[t] ~ dnorm(0,V[t])
                mu[t] <- mu[t-1]+d1[t]*nu[t]
                eps[t] <- gamma*eps[t-1] + u[t]
                V[t] <- V[t-1]*pow(omega,d2[t])
					# shift in variance mechanism
					}
predy263 <- mu[T] + eps[T]
predy262 <- mu[T-1] + eps[T-1]
predy261 <- mu[T-2] + eps[T-2]
mu259 <- mu[258]+d1[259]*nu[259]
mu260 <- mu[259]+d1[260]*nu[260]
V259 <- V[258]*pow(omega,d2[259])
V260 <- V[259]*pow(omega,d2[260])
# priors
tau.nu <- 10;  gamma ~ dnorm(0,0.1);
eta1 ~ dbeta(1,19); eta2 ~ dbeta(1,19)
omega ~ dexp(1)}"


# Running the model
library(rjags)
model5 <- jags.model(textConnection(model_string5),data = list(T=263,y=c(5.7,4.7,3.9,3.9,6.4,5.9,4.8,6.9,5.5,5.9,5.4,5.1,3.8,3.8,3.8,3.4,4.2,5.9,
                                                                         5.7,5,5.2,8.1,7.7,7.6,6.3,5.9,6.3,7.4,8.9,10.4,7.8,7.2,7.2,6.6,5.7,5.2,
                                                                         5.3,6.6,7.4,7.1,6.6,5.4,5.5,5.2,4.6,4.4,4.1,4.2,5.7,5.9,5.6,5.4,4.8,4.5,
                                                                         4.9,8.3,9.8,9,8.3,7.7,6.7,5.5,4.9,4.7,4.1,5.9,4.3,4.3,4.1,7.4,5.1,5.1,
                                                                         7.1,5.5,5.9,5.1,4.6,3.9,3.8,3.5,3.4,4.8,5.9,5.7,4.9,5.1,9,7.4,7,6,5.6,
                                                                         7.5,7.5,9.4,10.1,7.4,7.2,7.2,6.3,5.6,5.2,5.4,6.9,7.6,7.1,6.1,5.6,5.6,4.9,
                                                                         4.4,4.2,4,4.3,5.8,6.1,5.6,5.1,4.6,4.4,5.4,9.4,9.6,9,8.2,7.5,6.3,5.6,4.8,
                                                                         4.4,3.8,6,4.2,4.1,4.1,7.4,5.2,5.6,6.6,5.7,5.4,5,4.4,3.8,3.8,3.5,3.5,5.1,
                                                                         6.1,5.6,4.8,5.5,8.4,7.8,7,5.9,6,7.7,7.4,9.8,9.5,7.5,7.1,6.9,6,5.6,5.2,
                                                                         5.7,6.9,7.6,6.8,6,5.7,5.1,4.8,4.5,4.2,4.1,4.9,5.7,6.1,5.4,4.9,4.7,4.6,
                                                                         6.1,9.6,9.5,9,8.1,7.2,6.1,5.1,4.9,4.4,3.8,5.3,4.2,4.3,5.1,6.2,5.8,6.1,
                                                                         6.1,5.7,5.7,4.8,4.1,3.6,3.9,3.4,3.5,5.9,6,5.3,4.8,6.6,8.3,7.8,6.8,5.9,
                                                                         5.9,7.5,8.3,10.8,8.5,7.2,7,6.9,5.8,5.3,5.4,6.2,7,7.4,6.6,5.6,5.6,5.4,
                                                                         4.6,4.4,4.1,3.9,5.5,5.9,5.8,5.4,5,4.5,4.7,6.8,9.9,9.8,8.6,7.7,6.9,5.8,
                                                                         5.1,4.7,4.2,3.7,NA,NA,NA)), n.chains = 1, n.adapt= 10000)
                                                                                
                                                                
update(model5, 1000); # Burning for 1000 samples
mcmc_samples5 <- coda.samples(model5, variable.names=c("gamma","predy263","predy262","predy261","mu259","mu260","V259","V260"), n.iter=100000, thin = 100)



#### PLOT POSTERIOR STATISTICS
plot(mcmc_samples5)


#### ERGODIC MEAN
dev.new(width=16,height=12);par(mfrow=c(2,2))
M5<-as.matrix(mcmc_samples5)
plot(cumsum(M5[,1])/cumsum(rep(1,NROW(M5[,1]))),xlab = "Iterations",ylab = "mu",main = "gamma")


out <- capture.output(summary(mcmc_samples5))
cat("THstat5", out, file="stats5.txt", sep="\n", append=TRUE)

#### HPD CREDIBLE INTERVAL
library(HDInterval)
out <- capture.output(hdi(mcmc_samples5))
cat("HDI 5", out, file="stats5.txt", sep="\n", append=TRUE)




model_string4 <-"model{
  for (i in 1:n){
    y[i]~dpois(l1[i])
    J[i] <- step(i-k)
    l1[i]<- mu*J[i]+l*(1-J[i])
  }
  
  k~dcat(p[])
  for (j in 1:n){
    p[j]<-1/n
  }
  
  l<-U/g;U~dexp(1)
  g~dgamma(.5,1)
  
  mu<-V/d
  V~dexp(1)
  d~dgamma(.5,1)}"

# Running the model
model4 <- jags.model(textConnection(model_string4), data = list(n=172, y=c(1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0,
 0, 0, 1, 0, 0, 0,1, 0, 2, 3, 3, 1, 2, 2,0, 2, 2, 2, 3, 2, 3, 2, 2, 6, 5, 6, 3, 3,3, 5, 5, 4, 10, 6, 12, 8, 8, 5, 6, 5, 12, 12, 9, 9, 20, 15, 19, 12, 22, 25, 23, 16, 11, 16, 12, 13, 15, 6, 11, 12, 5, 17, 12, 11,
15, 17, 12,14, 9, 15, 6, 6, 2, 3, 4, 8, 7, 7,2, 6, 8, 9, 11, 4, 6, 2, 6, 6, 0, 1, 5, 4, 1, 2, 0, 0, 5, 3, 3, 1, 1, 2, 3, 0, 3,
2, 0, 3, 0,1, 1, 4, 0, 0, 0, 2, 1, 1, 1, 1, 3, 1, 1, 1, 0, 1,0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,0, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1)), n.chains = 1, n.adapt= 10000)

update(model4, 1000); # Burning for 1000 samples
mcmc_samples4 <- coda.samples(model4, variable.names=c("mu","l","k"), n.iter=100000, thin = 100)

summary(mcmc_samples4)


model_string5 <-"model  { 

#---------------Information to define y1 outside of the loop ----------------------------
V[1]<-0.1 # Here i just chose this value for V[1], Since no information is given on V0.
mu[1] ~ dnorm(0,1)
eps0~dnorm(0,1)
u1~dnorm(0,V[1])

eps[1]<- gamma*eps0 + u1
y1<-mu[1]+eps[1]


Ymu[1]<-mu[1]+gamma*eps0   # Because you are dealing with AR(1), define the y1, and mu[1] outside of the loop.

for (t in 2:T) {
y[t]~dnorm(Ymu[t], Ytau[t]) # Likelihood of y[t]

#Remember here that the distribution y[t] is actually a sum of two Normal distributions
#i.e mu[t]~Normal(mu[t-1], d1[t]* 1/tau_nu) and eps[t] ~Normal(gamma*eps[t-1], V[t]).
# For simplicity, i treated them as independent variables, therefore y[t] ~Normal(mu[t-1]+gamma*eps[t-1],d1[t]* 1/tau_nu +V[t] )
# If you are motivated enough and want to treat them as dependent, you need to add d1[t]* 1/tau_nu +V[t] + 2Cov(mu[t],eps[t])
# as the variance or the precision to Ytau below.

Ymu[t] <- mu[t-1]+gamma*eps[t-1] # THIS IS THE MEAN OF y[t]
Ytau[t]<- d1[t]*nu[t]+V[t] # THIS IS THE PRECISION OF y[t]

nu[t] ~ dnorm(0,tau.nu) # This is the distribution on NU
u[t] ~ dnorm(0,V[t]) # This is the distribution on NU

# What i am not sure is whether tan.nu and V[t] given in the question are the precision or the variance.
# Either way, be sure before you use it.

#--------Equations of the Relationship between mu[t],eps[t], and V[t]
eps[t] <- gamma*eps[t-1] + u[t]
V[t]<-V[t-1]*pow(omega, d2[t])
mu[t]<-mu[t-1] +d1[t]*nu[t]

d1[t]~dbern(eta1)  # This is the prior distribution given for delta1
d2[t] ~ dbern(eta2) # This is the prior distribution given for delta2

}


# priors
#tau.nu <- 0.1
tau.nu<-10
gamma ~ dnorm(0,0.1)
eta1 ~ dbeta(1,19)
eta2 ~ dbeta(1,19)
omega ~ dexp(1) 

}"
