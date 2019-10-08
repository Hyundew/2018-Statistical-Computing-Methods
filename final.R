# set working directory
setwd("/Volumes/Hyunjoo/4-1/[대학원] 통계계산방법론/final")

# call functions
source('mycode.R')

# set up
n.sim <- 10^4
n <- 100
p <- 1
alpha <- 0.05 # confidence level
B <- 2*10^3

# parameter
true.beta0 <- 3
true.beta1 <- 10
beta <- c(true.beta0, true.beta1)
sigma <- 1


result1 <- result2 <- result3 <- matrix(0, n.sim, 4)
for (iter in 1:n.sim) {
  tic <- Sys.time()

  # Generating data set
  set.seed(iter)
  x <- matrix(rnorm(n * p, 10, 10), n, p)
  x <- cbind(rep(1, n), x^2)
  
  eps <- rgamma(n, 1, 10^(-1))
  sigma <- 1/eps
  
  y <- x %*% beta + eps
  
  ####################
  # comparison of CI #
  ####################

  result1[iter, ] <- tryCatch(mmls(x, y, alpha), error = function(e) matrix(0, p + 1 , 4))
  result2[iter, ] <- tryCatch(ci.boot(x, y, alpha, B), error = function(e) matrix(0, p + 1, 4))
  result3[iter, ] <- tryCatch(gibbs(x, y, alpha), error = function(e) matrix(0, p + 1, 4))

  # 시간 계산
  toc <- Sys.time()
  # 한번 돌린 이후 얼마나 걸리는지 확인
  # cat= print
  cat(iter, "th iteration takes ", round(toc - tic, 2), " secs. \n", sep = "")
}

# save results
write(result1, "result/result1.txt")
write(result2, "result/result2.txt")
write(result3, "result/result3.txt")


beta

# number of simulation
result1 <- matrix(scan("result/result1.txt"), nrow = n.sim)
result2 <- matrix(scan("result/result2.txt"), nrow = n.sim)
result <- matrix(scan("result/result3.txt"), nrow = n.sim)

##########################
# function for summaries #
##########################

get.summary <- function(result, beta){
  estimate <- mse <- length <- length.sd <- NULL
  estimate <- mean(result[, 1])
  mse <- mean(result[, 2])
  length   <- mean(result[, 4] - result[, 3])
  length.sd   <- sd(result[, 4] - result[, 3])
  obj <- list(estimate=estimate, mse=mse, length=length, length.sd= length.sd, n.sim = nrow(result))
  return(obj)
}

obj1 <- get.summary(result1, beta)
obj2 <- get.summary(result2, beta)
obj3 <- get.summary(result3, beta)

estimate <- rbind(obj1$estimate, obj2$estimate, obj3$estimate); estimate
MSE <- rbind(obj1$mse, obj2$mse, obj3$mse); MSE
CI <- rbind(obj1$length, obj2$length, obj3$length); CI
CI.sd <- rbind(obj1$length.sd, obj2$length.sd, obj3$length.sd); CI.sd



library(xtable)
xtable(estimate, 3)
xtable(length, 3)
