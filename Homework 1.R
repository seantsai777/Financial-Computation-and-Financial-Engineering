rm(list=ls())

#Basic requirement
#input
S0 <- 100
r <- 0.05
q <- 0.02
sigma <- 0.3
t <- 1
K1 <- 85
K2 <- 90
K3 <- 120
K4 <- 130
  
d1 <- (log(S0/K1)+(r-q-0.5*sigma^2)*t)/(sigma*sqrt(t))
d2 <- (log(S0/K2)+(r-q-0.5*sigma^2)*t)/(sigma*sqrt(t))
d3 <- (log(S0/K3)+(r-q-0.5*sigma^2)*t)/(sigma*sqrt(t))
d4 <- (log(S0/K4)+(r-q-0.5*sigma^2)*t)/(sigma*sqrt(t))

#Output
option_value <- exp(-r*t)*
                ((S0*exp((r-q)*t))*(pnorm(d1+sigma*sqrt(t))-pnorm(d2+sigma*sqrt(t)))-K1*(pnorm(d1)-pnorm(d2))
                +(K2-K1)*(pnorm(d2)-pnorm(d3))
                +((K2-K1)/(K3-K4))*(S0*exp((r-q)*t)*(pnorm(d3+sigma*sqrt(t))-pnorm(d4+sigma*sqrt(t)))-K4*(pnorm(d3)-pnorm(d4))))

#Bonus
n <- 10000
times <- 20
option_v <- NULL

for (k in 1:times) {
  ST <- S0*exp((r-q-0.5*sigma^2)*t+sigma*sqrt(t)*qnorm(runif(n)))
  payoff <- NULL

  for (i in c(1:length(ST))) {
    payoff1 <- max(ST[i]-K1,0)
    payoff2 <- K2-K1
    payoff3 <- max(((K2-K1)/(K3-K4))*(ST[i]-K4),0)
    if (ST[i] < K2) {p <- payoff1
    } else if (ST[i] >= K2 && ST[i] <= K3) {p <- payoff2
    } else {p <- payoff3}
    
  payoff <- c(payoff, p)
  }

  v <- exp(-r*t)*mean(payoff)
  option_v <- c(option_v, v)
}

value_monte <- mean(option_v)
value_var <- mean((option_v-value_monte)^2)
value_CI <- c(value_monte-2*sqrt(value_var),value_monte+2*sqrt(value_var))

