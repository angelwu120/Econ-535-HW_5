#2

k <- 10
n <- 200
#n <- 50  #for (b), rerun the experiments with n=50
N <- 200

## Generating beta matrix
beta0 <- 1
beta <- matrix(data = 0, nrow = 5, ncol = k)
for(i in 1:4) {
  beta[i+1,i] <- 1
}
beta[5,1] <- 0.5

## Generating X_i
mu <- matrix(data = 0, nrow = 1, ncol = k)

I_10 <- diag(k)
l <- matrix(data = 1, nrow = k, ncol = 1)

rho <- c(0,0.9)

sigma1 <- I_10 + (l %*% t(l) - I_10)*rho[1]
sigma2 <- I_10 + (l %*% t(l) - I_10)*rho[2]
X1 <- rmvnorm(n = n, mean = mu, sigma = sigma1)
X2 <- rmvnorm(n = n, mean = mu, sigma = sigma2)

# Initial setup ---- FLAG
coeff <- array(0:0, c(k+1,10,N))
tValue <- array(0:0, c(k+1,10,N))
pValue <- array(0:0, c(k+1,10,N))

for(j in 1:N) {
  ## Generating epsilon
  epsilon <- rnorm(n = n, mean = 0, sd = 1)
  
  ## Calculating y
  y <- array(0:0, c(n,10))
  #1~5 with rho=0
  for(k in 1:5){
    y[,k] <- beta0 + X1 %*% beta[k,] + epsilon
  } 
  #6~10 with rho=0.9
  for(k in 6:10){
    y[,k] <- beta0 + X2 %*% beta[k-5,] + epsilon
  } 
  
  ## Estimating beta
  fit1 <- lm(formula = y[,1] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + X1[,6] + X1[,7] + X1[,8] + X1[,9] + X1[,10])
  fit2 <- lm(formula = y[,2] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + X1[,6] + X1[,7] + X1[,8] + X1[,9] + X1[,10])
  fit3 <- lm(formula = y[,3] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + X1[,6] + X1[,7] + X1[,8] + X1[,9] + X1[,10])
  fit4 <- lm(formula = y[,4] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + X1[,6] + X1[,7] + X1[,8] + X1[,9] + X1[,10])
  fit5 <- lm(formula = y[,5] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + X1[,6] + X1[,7] + X1[,8] + X1[,9] + X1[,10])
  fit6 <- lm(formula = y[,6] ~ X2[,1] + X2[,2] + X2[,3] + X2[,4] + X2[,5] + X2[,6] + X2[,7] + X2[,8] + X2[,9] + X2[,10])
  fit7 <- lm(formula = y[,7] ~ X2[,1] + X2[,2] + X2[,3] + X2[,4] + X2[,5] + X2[,6] + X2[,7] + X2[,8] + X2[,9] + X2[,10])
  fit8 <- lm(formula = y[,8] ~ X2[,1] + X2[,2] + X2[,3] + X2[,4] + X2[,5] + X2[,6] + X2[,7] + X2[,8] + X2[,9] + X2[,10])
  fit9 <- lm(formula = y[,9] ~ X2[,1] + X2[,2] + X2[,3] + X2[,4] + X2[,5] + X2[,6] + X2[,7] + X2[,8] + X2[,9] + X2[,10])
  fit10 <- lm(formula = y[,10] ~ X2[,1] + X2[,2] + X2[,3] + X2[,4] + X2[,5] + X2[,6] + X2[,7] + X2[,8] + X2[,9] + X2[,10])
  
  c <- cbind(fit1$coefficients,
             fit2$coefficients,
             fit3$coefficients,
             fit4$coefficients,
             fit5$coefficients,
             fit6$coefficients,
             fit7$coefficients,
             fit8$coefficients,
             fit9$coefficients,
             fit10$coefficients)
  
  t <- cbind(as.data.frame(summary(fit1)$coefficients)$'t value',
             as.data.frame(summary(fit2)$coefficients)$'t value',
             as.data.frame(summary(fit3)$coefficients)$'t value',
             as.data.frame(summary(fit4)$coefficients)$'t value',
             as.data.frame(summary(fit5)$coefficients)$'t value',
             as.data.frame(summary(fit6)$coefficients)$'t value',
             as.data.frame(summary(fit7)$coefficients)$'t value',
             as.data.frame(summary(fit8)$coefficients)$'t value',
             as.data.frame(summary(fit9)$coefficients)$'t value',
             as.data.frame(summary(fit10)$coefficients)$'t value')
  
  p <- cbind(as.data.frame(summary(fit1)$coefficients)$'Pr(>|t|)',
             as.data.frame(summary(fit2)$coefficients)$'Pr(>|t|)',
             as.data.frame(summary(fit3)$coefficients)$'Pr(>|t|)',
             as.data.frame(summary(fit4)$coefficients)$'Pr(>|t|)',
             as.data.frame(summary(fit5)$coefficients)$'Pr(>|t|)',
             as.data.frame(summary(fit6)$coefficients)$'Pr(>|t|)',
             as.data.frame(summary(fit7)$coefficients)$'Pr(>|t|)',
             as.data.frame(summary(fit8)$coefficients)$'Pr(>|t|)',
             as.data.frame(summary(fit9)$coefficients)$'Pr(>|t|)',
             as.data.frame(summary(fit10)$coefficients)$'Pr(>|t|)')
  
  coeff[,,j] <- c
  tValue[,,j] <- t
  pValue[,,j] <- p
  
}

## Computing power1: t>1.96
t_1 <- array(0:0, c(4,10,N))
power1 <- array(0:0,c(4,10))
for(i in 2:5){
  for(j in 1:10){
    for(k in 1:N){
      if(tValue[i,j,k] > 1.96 ){
        t_1[i-1,j,k] <- 1
      }
    }
    power1[i-1,j] <- sum(t_1[i-1,j,])/N  # power(rho=0) > power(rho=0.9)
  }
}
power1

## Computing power2: t>2.81
t_2 <- array(0:0, c(4,10,N))
power2 <- array(0:0,c(4,10))
for(i in 2:5){
  for(j in 1:10){
    for(k in 1:N){
      if(tValue[i,j,k] > 2.81 ){
        t_2[i-1,j,k] <- 1
      }
    }
    power2[i-1,j] <- sum(t_2[i-1,j,])/N
  }
}
power2



