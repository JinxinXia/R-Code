library(ggplot2)


## load data
library(Stat2Data)
data(Putts1)
n=dim(Putts1)[1]
X=cbind(rep(1,n),Putts1$Length)
y=Putts1$Made

data1 <- cbind.data.frame(y,X[,2])
names(data1) <- c('Made','Length')


## Jitter adds a bit of random noise, to see points 
## that would otherwise fall on top of one another.
## Observe a tendency for putts to be made at short
## lengths, and missed at longer lengths.





## Make a scatterplot of Made vs. Length using the command >jitter() to add
## a bit of random scatter to the points to see the density of points. Add
## the empirical proportion made for each length to the plot. 

## Fit three GLMS: use the identity, logit, and probit link functions.
## response variable follows Bernoulli distribution


# 
#################################################################################################
# keep the original data untouched
Putts2 <- Putts1


# get Made proportion for each Length
uniq_Lenght <- unique(Putts2$Length)
len_uniq_Length <- length(uniq_Lenght)

uniq_Length_Made_proportion <- rep(0,len_uniq_Length)
uniq_Length_Made_success <- rep(0,len_uniq_Length)

for(i in 1:len_uniq_Length){
  uniq_Length_Made_success[i] <- sum(Putts2[which(Putts2$Length==uniq_Lenght[i]),]$Made)
  uniq_Length_Made_proportion[i] <- uniq_Length_Made_success[i]/length(Putts2[which(Putts2$Length==uniq_Lenght[i]),]$Made)
}

uniq_logit_prop <- log(uniq_Length_Made_proportion/(1-uniq_Length_Made_proportion))
emp_proportion <- as.data.frame(cbind(uniq_Lenght,uniq_Length_Made_proportion))
logit_emp_proportion <- as.data.frame(cbind(uniq_Lenght,uniq_logit_prop))


# fit glm models with logit, identity, probit
glm_logit <- glm(Made ~ Length, data = data1, family=binomial(link="logit"))
summary(glm_logit)

glm_identity <- glm(Made ~ Length, data = data1, family=binomial(link="identity"))
summary(glm_identity)

glm_probit <- glm(Made ~ Length, data = data1, family=binomial(link="probit"))
summary(glm_probit)

# construct new data frame for predicted values from 0 to 15
new_dat <- data.frame(Length = seq(0,15,0.1))
new_dat$logit_pred <- predict(glm_logit,newdata = new_dat,type='response')
new_dat$identity_pred <- predict(glm_identity,newdata = new_dat,type='response')
new_dat$probit_pred <- predict(glm_probit,newdata = new_dat,type='response')


# plot with only proportion
glmlines_Made_vs_Len_pro <- ggplot(data = Putts2,aes(x=jitter(Length),y=jitter(Made))) + geom_point() +
  geom_point(data=emp_proportion,aes(x=uniq_Lenght,y=uniq_Length_Made_proportion,color='gold')) +
  scale_color_discrete(name = "" ,labels = " ") + 
  labs(x="Length",y="Made",title = "Length Vs. Made and Empricial Proportion of Made") + 
  theme(plot.title = element_text(hjust = 0.5))
glmlines_Made_vs_Len_pro



# combine glm lines with scatter plots
glmlines_Made_vs_Len_glms <- ggplot(data = Putts2,aes(x=jitter(Length),y=jitter(Made))) + geom_point() +
  geom_line(data=new_dat,aes(x=Length,y=logit_pred,color='red'),size=1) + 
  geom_line(data=new_dat,aes(x=Length,y=identity_pred,color='blue'),size=1) + 
  geom_line(data=new_dat,aes(x=Length,y=probit_pred,color='grey'),size=1) + 
  labs(x="Length",y="Made",title = "Length Vs. Made with GLM Lines ") + 
  scale_color_discrete(name = "link functions",labels = c("identity", "probit","logit",'Made Success proportion')) +
  xlim(0,15) +
  theme(plot.title = element_text(hjust = 0.5))
glmlines_Made_vs_Len_glms


# plot the link function
#################################################################################################

glmlines_logit_Made_vs_Length <- ggplot() + geom_point(data = logit_emp_proportion,aes(x=uniq_Lenght,y=uniq_logit_prop,color='blue')) +
  geom_line(data=new_dat,aes(x=Length,y=logit_pred,color='red'),size=1) + 
  labs(x="Length",y="Logit Made",title = "Logit Made Vs. Length with GLM Lines")  +
  scale_color_discrete(name = "", labels = c("logit proportions", "fitted model")) +
  theme(plot.title = element_text(hjust = 0.5))
glmlines_logit_Made_vs_Length


# Fisher scoring J 
#################################################################################################


U1 <- function(beta){
  score = c(0,0)
  x <- Putts1$Length
  y <- Putts1$Made
  p <- exp(beta[1]+beta[2]*x)/(exp(beta[1]+beta[2]*x)+1)
  score[1] <- sum(y-p)
  score[2] <- sum((y-p)*x)
  return(score)
}

J1 <- function(beta){
  obs_info <- rbind(c(0,0),c(0,0))
  x <- Putts1$Length
  p <- exp(beta[1]+beta[2]*x)/(exp(beta[1]+beta[2]*x)+1)
  obs_info[1,1] <- sum(p*(1-p))
  obs_info[1,2] <- sum(p*(1-p)*x)
  obs_info[2,1] <- sum(p*(1-p)*x)
  obs_info[2,2] <- sum(p*(1-p)*x^2)
  return(obs_info)
}

num_store <- 15
beta_store <- cbind(rep(0,num_store),rep(0,num_store))


itera_beta <- ggplot() + geom_point(aes(x=1,y=0)) 

beta_finding_1 <- function(beta){
  steps <- 1
  cat(sprintf("iteration      beta_0        beta_1        Score_0           Score_1            U_res\n"))
  U_res <- 1
  while(U_res >1e-18 & steps<1+num_store){
    beta <- beta + MASS::ginv(J1(beta))%*%U1(beta)
    steps <- steps + 1
    U_res <- U1(beta)[1]^2 + U1(beta)[2]^2
    cat(sprintf("%6.2d         %6.4f        %6.4f        %6.4f          %6.4f         %6.4f \n", steps,beta[1],beta[2],U1(beta)[1],U1(beta)[2],U_res))
  }
  return(beta)
}


beta_store <- as.data.frame(cbind(c(3.5018,3.2371,3.2568,3.2568,3.2568,3.2568,3.2568,3.2568,3.2568,3.2568 ),
                                  c(-0.6425,-0.5616,-0.5661,-0.5661,-0.5661,-0.5661,-0.5661,-0.5661,-0.5661,-0.5661),
                                  seq(1,10,1)))

iter_beat <- ggplot(data = beta_store) + geom_point(aes(x=V3,y=V1,color='red')) 
iter_beat + geom_point(data=beta_store,aes(x=V3,y=V2,color='blue')) + 
  labs(x="iterations",y="Beta value",title = "Beta vs Iteration ")  +
  scale_color_discrete(name = "Beta", labels = c("beta1", "beta0")) +
  theme(plot.title = element_text(hjust = 0.5))

beta_finding_1(c(1,0))


# construct confidence interval
#################################################################################################

var.beta <- solve(J1(c(3.2568,-0.5661)))

#beta0
cbind(3.2568-1.96*sqrt(var.beta[1,1]),3.2568 + 1.96*sqrt(var.beta[1,1]))
#beta1
cbind(-0.5661-1.96*sqrt(var.beta[2,2]),-0.5661+1.96*sqrt(var.beta[2,2]))



## compute AIC for the logit linked model
#################################################################################################

log_likeli <- function(beta){
  x <- Putts1$Length
  y <- Putts1$Made
  p <- exp(beta[1]+beta[2]*x)/(1+exp(beta[1]+beta[2]*x))
  llike <-  sum(y*log(p) + (1-y)*log(1-p))
  return(llike)
}

logl_full <- log_likeli(c(glm_logit$coefficients[1],glm_logit$coefficients[2]))


# K =2 
AIC_logl_full = -2*logl_full + 2*2
AIC_logl_full
glm_logit$aic

# for onlu intercept, k = 1 
glm_logit_redu <- glm(Made ~ 1, data = data1, family=binomial(link="logit"))
logl_red <- log_likeli(c(glm_logit_redu$coefficients[1],0))
AIC_logl_red = -2*logl_full + 2*1
AIC_logl_red
glm_logit_redu$aic






##############################################################################################
