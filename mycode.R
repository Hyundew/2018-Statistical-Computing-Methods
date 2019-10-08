rm(list=ls())

# original setting
set.seed(1)
n.sim <- 10^3
alpha <- 0.05 # confidence level

n <- 100
p <- 1

sigma <- 1
beta0 <- 3
beta1 <- 10
beta <- c(beta0, beta1)

x <- matrix(rnorm(n * p, 10, 10), n, p)
x <- cbind(rep(1, n), x^2)

y <- x %*% beta + sigma


# frequentist
# mmls
mmls <- function(x, y, max.iter = 100, eps = 1.0e-8, alpha = 0.05)
{
  beta <- c(0, 0)
  iter <- 1
  for (iter in 1:max.iter){
    a <- abs(x) / apply(abs(x), 1, sum)
    r <- c(y - x %*% beta) # y - t(x) %*% beta
    temp1 <- apply(x * r, 2, sum) # x(y-t(x) %*% beta)
    temp2 <- apply(x^2/a, 2, sum) 
    beta.new <- beta + temp1/temp2
    if (max(abs(beta.new - beta)) < eps) break
    beta <- beta.new
  }
  
  mse <- sum((y - x %*% beta.new)^2) / (n - ncol(x) - 1)

  # standard errors
  chol.obj <- chol(t(x) %*% x)
  inv.XX <- chol2inv(chol.obj)
  se <- sqrt(diag(inv.XX * mse))
  z.alpha <- qnorm(1 - alpha/2)

  lo <- beta.new[2] - z.alpha * se[2]
  up <- beta.new[2] + z.alpha * se[2]
  ci <- cbind(lo, up)

  obj <- c("estimate"= beta.new[2], "mse"=mse, "95% CI"=ci)
  return(obj)
}
result1 <- mmls(x, y, alpha=0.05)

# bootstrap
Householder <- function(x) 
{
  s <- sign(x[1]) * c(sqrt(t(x) %*% x)) # 부호는 상관없음
  nn <- length(x)
  e1 <- rep(0, nn); e1[1] <- 1
  u <- x + s * e1  
  d <- 2 / c(t(u) %*% u) # scalar
  U <- diag(nn) - outer(d*u, u)
  list(U = U, u = u/s)
}

housereg <- function(X, y) {
  p <- ncol(X)
  n <- nrow(X)
  Xy <- cbind(X, y)
  
  u <- as.list(1:p)   
  # list, dataframe 지양
  # u 저장
  
  # initial
  temp.hh <- Householder(X[,1])
  U <- temp.hh$U
  u[[1]] <- temp.hh$u
  
  
  for (i in 2:p) {
    x <- (U %*% X)[i:n,i]  # U, u를 구하고 나서 새로운 x update
    temp.hh <- Householder(x)
    u[[i]] <- temp.hh$u
    temp.U <- temp.hh$U
    temp.I <- diag(i-1)
    temp.01 <- matrix(0, i-1, n - i +1)
    temp.02 <- matrix(0, n - i + 1, i-1)
    U <- rbind(cbind(temp.I, temp.01), cbind(temp.02, temp.U)) %*% U
  }
  
  temp <- U %*% Xy
  R <- temp[1:p, 1:p]
  z1 <- temp[1:p, p+1]  # X & y cbind해서 y부분 계산 한번에 처리됨
  z2 <- temp[-(1:p), p+1]
  
  coef <- backsolve(R, z1)
  return(coef)
}

ci.boot <- function(x, y, alpha=0.05, B=2*10^3)
{
  
  p <- ncol(x)
  boot.beta <- matrix(0, B, p)
  
  # Bootstrapping - bias corrected
  for (b in 1:B)
  {
    id <- sample(n, replace = T)
    boot.x <- x[id,]
    boot.y <- y[id]
    boot.beta[b, ] <- boot.obj <- housereg(boot.x, boot.y)
  }
  
  hat.beta <- housereg(boot.x, boot.y)
  z.alpha <- qnorm(1-alpha/2)
  
  ci.boot <- matrix(0, p, 2)
  for (j in 1:p)
  {
    p0 <- mean(boot.beta[,j] <= hat.beta[j])
    z0 <- qnorm(p0)
    p.lo <- pnorm(2 * z0 - z.alpha)
    p.up <- pnorm(2 * z0 + z.alpha)
    ci.boot[j,] <- quantile(boot.beta[,j], prob = c(p.lo, p.up))
  }
  
  bb <- apply(boot.beta, 2, mean)
  mse <- sum((y- x %*% bb)^2) / (n-p-1)
   
  obj <- c('estimate'= bb[2], 'mse'= mse, '95% CI'=ci.boot[2, ])
  
  return(obj)
}
result2 <- ci.boot(x, y, alpha=0.05, B=10^3)


# bayesian
# gibbs
gibbs <- function(x, y, max.iter = 100, eps = 1.0e-8, alpha=0.05){
  # hyperparameters
  a <- 0.0001
  b <- 0.0001
  
  tau <- rep(1000, 2)
  mu <- rep(0, 2)
  
  # intial values:
  sigma2 <- 1
  beta <- rep(10, 2)
  
  n.samples <- 10000
  samples <- matrix(0, n.samples, 3)
  colnames(samples) <- c("beta0", "beta1", "sigma2")
  
  # Gibs_sampler
  for(i in 1:n.samples){
    # update sigma2:
    SSE <- sum((y - x %*% beta)^2)
    sigma2 <- 1/rgamma(1, n/2 + a, SSE/2 + b)
    
    # update beta1:
    v <- n/sigma2 + 1/tau[1]^2
    m <- sum(y - x[, 2] * beta[2])/sigma2 + mu[1]/tau[1]^2
    beta[1] <- rnorm(1, m/v, 1/sqrt(v))
    
    # update beta2:
    v <- sum(x^2)/sigma2 + 1/tau[2]^2
    m <- sum(x[, 2] * (y - beta[1]))/sigma2 + mu[2]/tau[2]^2
    beta[2] <- rnorm(1,m/v,1/sqrt(v))
    
    samples[i,] <- c(beta, sigma2)
  }
  
  post.sample <- samples[-(1:100), 1:2]
  est <- apply(post.sample, 2, mean)
  mse <- mse <- sum((y - x %*% est)^2) / (n - ncol(x) - 1)
  hat.beta1 <- post.sample[, 2]
  
  # credible intervals
  cppost <- cumsum(hat.beta1) / sum(hat.beta1)
  hi <- hat.beta1[which(cppost >= 0.025)[1]]
  lo <- hat.beta1[which(cppost >= 0.975)[1]-1]
  ci <- cbind(lo, hi)
 
  est <- cbind('estimate'=est[2], 'mse'=mse, '95% credible interval'= ci)
  return(est)
}
result3 <- gibbs(x, y, alpha=0.05)

list(result1, result2, result3)
