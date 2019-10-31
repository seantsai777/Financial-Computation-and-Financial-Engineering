rm(list=ls())

#Basic requirement
#input
S0 <- 50
Smax <- 50
r <- 0.1
q <- 0
sigma <- 0.4
t <- 0.25
n <- 10000
times <- 20
steps <- 10

dt <- t/steps
u <- exp(sigma*sqrt(dt))
d <- 1/u
p <- (exp((r-q)*dt)-d)/(u-d)
mu <- exp(r*dt)
qq <- (mu*u-1)/(mu*(u-d))

Am <- 1

#Output
#CRR binomial tree model
tree <- matrix(nrow=steps+1, ncol=steps+1)
tree[1,1] <- S0
for (j in 2:(steps+1)) {
  tree[1,j] <- S0*u^(j-1)
  tree[j,j] <- S0*d^(j-1)
}
for (j in 3:(steps+1)) {
  for (i in 2:(j-1)) {
    tree[i,j] <- tree[i-1,j-2]
  }
}
a=NULL
b=NULL
A <- array(0,c(2,ceiling((steps+1)/2),steps+1,steps+1))
A[1,1,1,1] <- Smax
for (h in 1:steps) {
  for (g in 1:h) {
    for (f in 1:ceiling((steps+1)/2)) {
      if (A[1,f,g,h]>=tree[g,h+1] && sum(a==A[1,f,g,h])==0 && sum(b==A[1,f,g,h])==0) {
        b <- c(b,A[1,f,g,h])
      }else if (A[1,f,g,h]>0 && A[1,f,g,h]<tree[g,h+1] && sum(a==tree[g,h+1])==0 && sum(b==tree[g,h+1])==0){
        b <- c(b,tree[g,h+1])
      }
    }
    A[1,,g,h+1] <- c(a,b,seq(0,0,length=(ceiling((steps+1)/2)-length(a)-length(b))))
    a=NULL
    b=NULL
    for (e in 1:ceiling((steps+1)/2)) {
      if (A[1,e,g,h]>=tree[g+1,h+1] & sum(a==A[1,e,g,h])==0 && sum(b==A[1,e,g,h])==0) {
        a <- c(a,A[1,e,g,h]) 
      }else if (A[1,e,g,h]>0 && A[1,e,g,h]<tree[g+1,h+1] && sum(a==tree[g+1,h+1])==0 && sum(b==tree[g+1,h+1])==0){
        a <- c(a,tree[g+1,h+1])
      }
    }
  }
  A[1,,g+1,h+1] <- c(a,seq(0,0,length=(ceiling((steps+1)/2)-length(a))))
  a=NULL
}
  
for (g in 1:(steps+1)) {
  for (f in 1:ceiling((steps+1)/2)) {
    if (A[1,f,g,steps+1]>0){A[2,f,g,steps+1] <- max(A[1,f,g,steps+1]-tree[g,steps+1],0)}
  }
}

for (h in steps:1) {
  for (g in 1:h) {
    for (f in 1:ceiling((steps+1)/2)){
      if (A[1,f,g,h]>0 && sum(A[1,,g,h+1]==A[1,f,g,h])==1){
        x <- which(A[1,,g,h+1]==A[1,f,g,h])
        y <- which(A[1,,g+1,h+1]==A[1,f,g,h])
        if (Am==0){A[2,f,g,h] <- exp(-r*dt)*(p*A[2,x,g,h+1]+(1-p)*A[2,y,g+1,h+1])
        }else {A[2,f,g,h] <- max(exp(-r*dt)*(p*A[2,x,g,h+1]+(1-p)*A[2,y,g+1,h+1]),A[1,f,g,h]-tree[g,h])}
      }else if (A[1,f,g,h]>0 && sum(A[1,,g,h+1]==A[1,f,g,h])==0){
        x <- which(A[1,,g,h+1]==tree[g,h+1])
        y <- which(A[1,,g+1,h+1]==A[1,f,g,h])
        if (Am==0){A[2,f,g,h] <- exp(-r*dt)*(p*A[2,x,g,h+1]+(1-p)*A[2,y,g+1,h+1])
        }else {A[2,f,g,h] <- max(exp(-r*dt)*(p*A[2,x,g,h+1]+(1-p)*A[2,y,g+1,h+1]),A[1,f,g,h]-tree[g,h])}
      }
    }
  }
}

#Monte Carlo simulation 
MC_v=NULL
for (i in 1:times) {
  ST <- matrix(nrow = n, ncol = steps+1)
  ST[,1] <- Smax
  ST[,2] <- S0*exp((r-q-0.5*sigma^2)*dt+sigma*sqrt(dt)*rnorm(n))
  for (k in 3:(steps+1)) {
    for (m in 1:n) {
      ST[m,k] <- ST[m,k-1]*exp((r-q-0.5*sigma^2)*dt+sigma*sqrt(dt)*rnorm(1))
    }
  }
  MC_max=NULL
  for (m in 1:n) {
    z <- max(ST[m,])
    MC_max <- c(MC_max,z)
  }
  v <- exp(-r*t)*mean(pmax(MC_max-ST[,steps+1],0))
  MC_v <- c(MC_v,v)
}
MC_value <- mean(MC_v)
value_sd <- sqrt(mean((MC_v-MC_value)^2))
value_CI <- c(MC_value-2*value_sd,MC_value+2*value_sd)

#Bonus 1

#Bonus 2
C <- matrix(nrow = steps+1, ncol = steps+1)
for (i in 1:(steps+1)) {
  C[i,i:(steps+1)] <- u^(i-1)
}

CP <- matrix(nrow = steps+1, ncol = steps+1)
CP[,steps+1] <- pmax(C[,steps+1]-1,0)
for (j in steps:1) {
  CP[1,j] <- exp(-r*dt)*(qq*CP[1,j+1]+(1-qq)*CP[2,j+1])
  if (j>=2){
    for (i in 2:j) {
      CP[i,j] <- exp(-r*dt)*(qq*CP[i-1,j+1]+(1-qq)*CP[i+1,j+1])
    }
  }
}
CP[1,1]*S0

