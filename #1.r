#1

alpha <- 1
beta <- matrix(data = c(0,1), nrow = 2, ncol = 1)

mu <- c(0,0)
I <- diag(x = 1,nrow = 2, ncol = 2)
n1 <- 50
n2 <- 100
N <- 300

## Generating X
install.packages("mvtnorm")
require("mvtnorm")
X.1 <- rmvnorm(n=n1, mean=mu, sigma = I)
X.2 <- rmvnorm(n=n2, mean=mu, sigma = I)

# Initial setup
c <- array(0:0, c(3,6,N))
t <- array(0:0, c(3,6,N))
p <- array(0:0, c(3,6,N))

for(i in 1:N) {
  ## Generating epsilon
  epsilon_1.1 <- rnorm(n=n1, mean = 0, sd = 1)
  epsilon_1.2 <- rnorm(n=n2, mean = 0, sd = 1)
  epsilon_2.1 <- rcauchy(n = n1, location = 0, scale = 1)
  epsilon_2.2 <- rcauchy(n = n2, location = 0, scale = 1)
  epsilon_3.1 <- rexp(n = n1, rate = 1) - 1
  epsilon_3.2 <- rexp(n = n2, rate = 1) - 1
  
  ## Calculating y
  y_1.1 <- alpha + X.1 %*% beta + epsilon_1.1
  y_1.2 <- alpha + X.2 %*% beta + epsilon_1.2
  y_2.1 <- alpha + X.1 %*% beta + epsilon_2.1
  y_2.2 <- alpha + X.2 %*% beta + epsilon_2.2
  y_3.1 <- alpha + X.1 %*% beta + epsilon_3.1
  y_3.2 <- alpha + X.2 %*% beta + epsilon_3.2

## Estimating beta
fit1.1 <- lm(formula = y_1.1 ~ X.1[,1] + X.1[,2])
fit1.2 <- lm(formula = y_1.2 ~ X.2[,1] + X.2[,2])
fit2.1 <- lm(formula = y_2.1 ~ X.1[,1] + X.1[,2])
fit2.2 <- lm(formula = y_2.2 ~ X.2[,1] + X.2[,2])
fit3.1 <- lm(formula = y_3.1 ~ X.1[,1] + X.1[,2])
fit3.2 <- lm(formula = y_3.2 ~ X.2[,1] + X.2[,2])

c[,,i] <- cbind(fit1.1$coefficients,
                fit1.2$coefficients,
                fit2.1$coefficients,
                fit2.2$coefficients,
                fit3.1$coefficients,
                fit3.2$coefficients)

t[,,i] <- cbind(as.data.frame(summary(fit1.1)$coefficients)$'t value',
                as.data.frame(summary(fit1.2)$coefficients)$'t value',
                as.data.frame(summary(fit2.1)$coefficients)$'t value',
                as.data.frame(summary(fit2.2)$coefficients)$'t value',
                as.data.frame(summary(fit3.1)$coefficients)$'t value',
                as.data.frame(summary(fit3.2)$coefficients)$'t value')

p[,,i] <- cbind(as.data.frame(summary(fit1.1)$coefficients)$'Pr(>|t|)',
                as.data.frame(summary(fit1.2)$coefficients)$'Pr(>|t|)',
                as.data.frame(summary(fit2.1)$coefficients)$'Pr(>|t|)',
                as.data.frame(summary(fit2.2)$coefficients)$'Pr(>|t|)',
                as.data.frame(summary(fit3.1)$coefficients)$'Pr(>|t|)',
                as.data.frame(summary(fit3.2)$coefficients)$'Pr(>|t|)')
}

## Plotting histograms
hist(t[3,1,]) #n=50,N(0,1)
hist(t[3,2,]) #n=100,N(0,1)
hist(t[3,3,]) #n=50,cauchy
hist(t[3,4,]) #n=100,cauchy
hist(t[3,5,]) #n=50,exp-1
hist(t[3,6,]) #n=100,exp-1

## Calculating quantiles
quant <- array(0:0,c(6,3))
for(i in 1:6) {
  quant[i,1] <- quantile(t[3,i,],0.1)
  quant[i,2] <- quantile(t[3,i,],0.05)
  quant[i,3] <- quantile(t[3,i,],0.01)
}
row.names(quant) <- c("n=50,N(0,1)",
                      "n=100,N(0,1)",
                      "n=50,cauchy",
                      "n=100,cauchy",
                      "n=50,exp-1",
                      "n=100,exp-1")
colnames(quant) <- c("10%","5%","1%")
quant
