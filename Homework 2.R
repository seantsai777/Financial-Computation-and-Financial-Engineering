rm(list=ls())

#Basic requirement
#input
S0 <- 100
r <- 0.05
q <- 0.01
sigma <- 0.5
t <- 0.4
K <- 80 
n <- 10000
times <- 20
steps <- 5

d1 <- (log(S0/K)+(r-q+0.5*sigma^2)*t)/(sigma*sqrt(t))
d2 <- d1-sigma*sqrt(t)
dt <- t/steps
u <- exp(sigma*sqrt(dt))
d <- 1/u
p <- (exp((r-q)*dt)-d)/(u-d)

#Output
#Black-Scholes formulas
BS_call <- S0*exp(-q*t)*pnorm(d1)-K*exp(-r*t)*pnorm(d2)
BS_put <- K*exp(-r*t)*pnorm(-d2)-S0*exp(-q*t)*pnorm(-d1)

#Monte Carlo simulation 
MC_v <- NULL
for (i in 1:times) {
  ST <- S0*exp((r-q-0.5*sigma^2)*t+sigma*sqrt(t)*rnorm(n))
  payoff <- NULL
  for (j in ST) {
      m <- c(max(j-K,0),max(K-j,0))
      payoff <- rbind(payoff,m)
  }
  v <- exp(-r*t)*colMeans(payoff)
  MC_v <- rbind(MC_v,v)
}
MC_value <- colMeans(MC_v)
setNames(MC_value,c("MC_call","MC_put"))

value_var <- colMeans((MC_v-matrix(MC_value,nrow=times,ncol=2,T))^2)
value_CI <- cbind(MC_value-2*sqrt(value_var),MC_value+2*sqrt(value_var))
rownames(value_CI) <- c("MC_call","MC_put")

#CRR binomial tree model
Am <- 0
Put <- 0
tree <- matrix(nrow=steps+1, ncol=steps+1)
for (h in 1:(steps+1)) {
  if (Put == 0)
  {tree[h,steps+1] <- max(S0*u^(steps+1-h)*d^(h-1)-K,0)}
  else {tree[h,steps+1] <- max(K-S0*u^(steps+1-h)*d^(h-1),0)}
}
for (j in steps:1) {
  for (i in 1:j) {
    tree[i,j] <- exp(-r*dt)*(p*tree[i,(j+1)]+(1-p)*tree[(i+1),(j+1)])
    if (Am == 1)
    {if (Put == 0) {tree[i,j] <- max(tree[i,j],S0*u^(j-i)*d^(i-1)-K)}
      else {tree[i,j] <- max(tree[i,j],K-S0*u^(j-i)*d^(i-1))}}
      
  }
}

#Bonus 1
tree_1col <- matrix(nrow=steps+1, ncol=1)
for (x in 1:(steps+1)) {
  if (Put == 0)
  {tree_1col[x] <- max(S0*u^(steps+1-x)*d^(x-1)-K,0)}
  else {tree_1col[x] <- max(K-S0*u^(steps+1-x)*d^(x-1),0)}
}
for (j in steps:1) {
  for (i in 1:j) {
    tree_1col[i] <- exp(-r*dt)*(p*tree_1col[i]+(1-p)*tree_1col[i+1])
    if (Am == 1)
    {if (Put == 0) {tree_1col[i] <- max(tree_1col[i],S0*u^(j-i)*d^(i-1)-K)}
      else {tree_1col[i] <- max(tree_1col[i],K-S0*u^(j-i)*d^(i-1))}}
  }
}

#Bonus 2
w <- 0
z <- 0
cnn <- 0
for (f in 1:(steps-1)) {
  if (f > steps/2){w <- steps-f} else {w <- f}
  for (g in steps:(steps+1-w)) {
    cnn <- cnn+log(g)-log(steps+1-g)
  }
  z <- z+exp(cnn+log(p)*(steps-f)+log(1-p)*f)*c(max(S0*u^(steps-f)*d^f-K,0),max(K-S0*u^(steps-f)*d^f,0))
  cnn <- 0
}
comb <- p^steps*c(max(S0*u^steps-K,0),max(K-S0*u^steps,0))+z+(1-p)^steps*c(max(S0*d^steps-K,0),max(K-S0*d^steps,0))
Comb_value <- exp(-r*t)*comb






