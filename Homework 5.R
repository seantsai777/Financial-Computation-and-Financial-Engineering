rm(list=ls())

#Basic requirement
#input
St <- 50
Savet <- 50
r <- 0.1
q <- 0
sigma <- 0.4
K <- 50
t <- 2
tt <- 1
M <- 100
steps <- 100
n <- 10000
times <- 20

dt <- t/steps
tsteps <- tt/dt
u <- exp(sigma*sqrt(dt))
d <- 1/u
p <- (exp((r-q)*dt)-d)/(u-d)

Am <- 0

#Output
#CRR binomial tree model
A <- array(0,c(2,M+1,steps+1,steps+1))
for (h in 1:(steps+1)) {
  for (g in 1:h) {
    A[1,1,g,h] <- (Savet*tsteps+St*(1-u^(h-g+1))/(1-u)+St*u^(h-g)*d*(1-d^(g-1))/(1-d))/(tsteps+h)
    A[1,M+1,g,h] <- (Savet*tsteps+St*(1-d^g)/(1-d)+St*d^(g-1)*u*(1-u^(h-g))/(1-u))/(tsteps+h)
    for (f in 2:M) {
      A[1,f,g,h] <- (M-(f-1))/M*A[1,1,g,h]+(f-1)/M*A[1,M+1,g,h]
    }
  }
}

for (g in 1:(steps+1)) {
  A[2,,g,steps+1] <- pmax(A[1,,g,steps+1]-K,0)
}

a = 1
b = 1
for (h in steps:1) {
  for (g in 1:h) {
    for (f in 1:(M+1)){
      Au <- ((tsteps+h)*A[1,f,g,h]+St*u^(h+1-g)*d^(g-1))/(tsteps+h+1)
      while(signif(Au,10) < signif(A[1,a,g,h+1],10)){a = a+1}
      if (a > 1) {
        Wu <- (A[1,a-1,g,h+1]-Au)/(A[1,a-1,g,h+1]-A[1,a,g,h+1])
        Cu <- Wu*A[2,a,g,h+1]+(1-Wu)*A[2,a-1,g,h+1]
      } else {Cu <- A[2,1,g,h+1]}
      
      
      Ad <- ((tsteps+h)*A[1,f,g,h]+St*u^(h-g)*d^g)/(tsteps+h+1)
      while(signif(Ad,10) < signif(A[1,b,g+1,h+1],10)){b = b+1}
      if (b > 1) {
        Wd <- (A[1,b-1,g+1,h+1]-Ad)/(A[1,b-1,g+1,h+1]-A[1,b,g+1,h+1])
        Cd <- Wd*A[2,b,g+1,h+1]+(1-Wd)*A[2,b-1,g+1,h+1]
      } else {Cd <- A[2,1,g+1,h+1]}
      
      
      if (Am == 1) {A[2,f,g,h] <- max(exp(-r*dt)*(p*Cu+(1-p)*Cd),A[1,f,g,h]-K)
      } else {A[2,f,g,h] <- exp(-r*dt)*(p*Cu+(1-p)*Cd)}
      
      a = 1
      b = 1
    }
  }
}
A[2,1,1,1]

#Monte Carlo simulation 
MC_v=NULL
for (i in 1:times) {
  ST <- matrix(nrow = n, ncol = steps)
  ST[,1] <- St*exp((r-q-0.5*sigma^2)*dt+sigma*sqrt(dt)*rnorm(n))
  for (k in 2:steps) {
    for (m in 1:n) {
      ST[m,k] <- ST[m,k-1]*exp((r-q-0.5*sigma^2)*dt+sigma*sqrt(dt)*rnorm(1))
    }
  }
  MC_ave <- (rowSums(ST)+St+Savet*tsteps)/(steps+1+tsteps)
  
  v <- exp(-r*t)*mean(pmax(MC_ave-K,0))
  MC_v <- c(MC_v,v)
}
MC_value <- mean(MC_v)
value_sd <- sqrt(mean((MC_v-MC_value)^2))
value_CI <- c(MC_value-2*value_sd,MC_value+2*value_sd)
