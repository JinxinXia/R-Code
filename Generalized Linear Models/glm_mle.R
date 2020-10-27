
## iterative method for computing maximum likelihood estimators of glm

# y
time=c(1.67,2.20,2.51,3.00,2.90,4.70,7.53,14.70,27.8,37.4,
       .8,1,1.37,2.25,2.95,3.70,6.07,6.65,7.05,7.37,
       .102, .18, .2, .24, .26, .32, .32, .42, .44, .88, 
       .073, .098, .117, .135, .175, .262, .270, .350, .386, .456)

# x
stress = c(rep(.87,10), rep(.99,10), rep(1.09,10), rep(1.18,10))


par(mfrow=c(1,2))
hist(time, breaks=20)
plot(stress, time)
dev.off()

glm1=glm(time ~ stress, family=Gamma(link="inverse"))
summary(glm1, dispersion=1)
summary(glm1, dispersion=1)$coef[,1] ## MLE
vcov(glm1) ## Variance

# y:time x:stress

U <- function(beta){
  score = c(0,0)
  lambda <- beta[1]+beta[2]*stress
  score[1] <- sum(1/lambda-time)
  score[2] <- sum((1/lambda-time)*stress)
  return(score)
}

J <- function(beta){
  obs_info <- rbind(c(0,0),c(0,0))
  lambda <- beta[1]+beta[2]*stress
  obs_info[1,1] <- sum(1/lambda^2)
  obs_info[1,2] <- sum(stress/lambda^2)
  obs_info[2,1] <- sum(stress/lambda^2)
  obs_info[2,2] <- sum(stress^2/lambda^2)
  return(obs_info)
}

beta_finding_fisher_scoring <- function(beta){
  steps <- 0
  cat(sprintf("iteration         beta_0           beta_1             Score_0                Score_1 \n"))
  while(U(beta) != 0){
    beta <- beta + solve(J(beta))%*%U(beta)
    steps <- steps + 1
    cat(sprintf("%6.2d          %6.4f          %6.4f            %6.4f                %6.4f \n", steps,beta[1],beta[2],U(beta)[1],U(beta)[2]))
  }
  return(beta)
}
beta_finding_fisher_scoring(c(1,0))