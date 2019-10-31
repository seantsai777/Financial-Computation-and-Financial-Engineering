rm(list=ls())

#Basic requirement
#input
K <- 50
r <- 0.1
t <- 0.5
S0 <- c(100,120,130,150)
q <- c(0.05,0.02,0.03,0.04)
sigma <- c(0.4,0.2,0.3,0.5)
rho <- matrix(c(1.0,0.1,0.2,0.3,
                0,1.0,0.4,0.5,
                0,0,1.0,0.6,
                0,0,0,1.0),4,4,T)

n <- 4
times <- 20

#Output
#Monte Carlo simulation
C <- matrix(nrow = length(S0),ncol = length(S0))
for (i in 1:length(S0)) {
  for (j in 1:length(S0)) {
    if(i>j){rho[i,j] <- rho[j,i]}
    C[i,j] <- rho[i,j]*sigma[i]*sigma[j]
  }
}

A <- matrix(0,nrow = length(S0),ncol = length(S0))
A[1,1] <- sqrt(C[1,1])
A[1,2:length(S0)] <- C[1,2:length(S0)]/A[1,1]
for (g in 2:(length(S0)-1)) {
  A[g,g] <- sqrt(C[g,g]-sum(A[1:(g-1),g]^2))
  for (h in (g+1):length(S0)) {
    A[g,h] <- 1/A[g,g]*(C[g,h]-sum(A[1:(g-1),g]*A[1:(g-1),h]))
  }
}
A[length(S0),length(S0)] <- sqrt(C[length(S0),length(S0)]-sum(A[1:(length(S0)-1),length(S0)]^2))

Rb <- NULL
for (x in 1:times) {
  Z <- matrix(nrow = 10000,ncol = n)
  for (z in 1:n) {
    set.seed(10*x+z)
    Z[,z] <- rnorm(10000)
  }
  R <- Z %*% A

  ST <- matrix(nrow = 10000,ncol = n)
  for (y in 1:10000) {
   ST[y,] <- S0*exp((r-q-0.5*sigma^2)*t+sigma*sqrt(t)*R[y,]) 
  }
  
  payoff <- NULL
  for (a in 1:nrow(ST)) {
    m <- max(max(ST[a,])-K,0)
    payoff <- c(payoff,m)
  }
  v <- exp(-r*t)*mean(payoff)
  Rb <- c(Rb,v)
}
Rb_value <- mean(Rb)
value_var <- mean((Rb-Rb_value)^2)
value_CI <- c(Rb_value-2*sqrt(value_var),Rb_value+2*sqrt(value_var))

#Bonus 1
Rb_1 <- NULL
for (x in 1:times) {
  Z_1 <- matrix(nrow = 10000,ncol = n)
  for (z in 1:n) {
    set.seed(10*x+z)
    w <- rnorm(5000)
    Z_1[,z] <- c(w,-w)/sqrt(mean(c(w,-w)^2))
  }
  R_1 <- Z_1 %*% A
  
  ST_1 <- matrix(nrow = 10000,ncol = n)
  for (y in 1:10000) {
    ST_1[y,] <- S0*exp((r-q-0.5*sigma^2)*t+sigma*sqrt(t)*R_1[y,]) 
  }
  
  payoff_1 <- NULL
  for (a in 1:nrow(ST_1)) {
    m <- max(max(ST_1[a,])-K,0)
    payoff_1 <- c(payoff_1,m)
  }
  v <- exp(-r*t)*mean(payoff_1)
  Rb_1 <- c(Rb_1,v)
}
Rb_value_1 <- mean(Rb_1)
value_var_1 <- mean((Rb_1-Rb_value_1)^2)
value_CI_1 <- c(Rb_value_1-2*sqrt(value_var_1),Rb_value_1+2*sqrt(value_var_1))

#Bonus 2
Rb_2 <- NULL
for (x in 1:times) {
  Z_2 <- matrix(nrow = 10000,ncol = n)
  for (z in 1:n) {
    set.seed(10*x+z)
    w <- rnorm(5000)
    Z_2[,z] <- c(w,-w)/sqrt(mean(c(w,-w)^2))
  }
  Z_sigma <- colMeans((Z_2-matrix(colMeans(Z_2),nrow=10000,ncol=n,T))^2)
  ZC <- matrix(nrow = n,ncol = n)
  for (d in 1:n) {
    for (e in 1:n) {
      ZC[d,e] <- cor(Z_2[,d],Z_2[,e])*Z_sigma[d]*Z_sigma[e]
    }
  }
  ZA <- matrix(0,nrow = n,ncol = n)
  ZA[1,1] <- sqrt(ZC[1,1])
  ZA[1,2:n] <- ZC[1,2:n]/ZA[1,1]
  for (g in 2:(n-1)) {
    ZA[g,g] <- sqrt(ZC[g,g]-sum(ZA[1:(g-1),g]^2))
    for (h in (g+1):n) {
      ZA[g,h] <- 1/ZA[g,g]*(ZC[g,h]-sum(ZA[1:(g-1),g]*ZA[1:(g-1),h]))
    }
  }
  ZA[n,n] <- sqrt(ZC[n,n]-sum(ZA[1:(n-1),n]^2))
  
  R_2 <- Z_2 %*% solve(ZA) %*% A
  
  ST_2 <- matrix(nrow = 10000,ncol = n)
  for (y in 1:10000) {
    ST_2[y,] <- S0*exp((r-q-0.5*sigma^2)*t+sigma*sqrt(t)*R_2[y,]) 
  }
  
  payoff_2 <- NULL
  for (a in 1:nrow(ST_2)) {
    m <- max(max(ST_2[a,])-K,0)
    payoff_2 <- c(payoff_2,m)
  }
  v <- exp(-r*t)*mean(payoff_2)
  Rb_2 <- c(Rb_2,v)
}
Rb_value_2 <- mean(Rb_2)
value_var_2 <- mean((Rb_2-Rb_value_2)^2)
value_CI_2 <- c(Rb_value_2-2*sqrt(value_var_2),Rb_value_2+2*sqrt(value_var_2))



r3<-R_2 <- Z_2 %*% solve(ZA)
var(r3)
