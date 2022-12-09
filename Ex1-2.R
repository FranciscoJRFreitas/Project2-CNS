## 1.2
survivalTimes = c(1552, 627, 884, 2183, 1354, 1354, 1014, 2420,  71, 3725,
	          2195, 2586, 1577, 1766, 1325, 1299, 159, 1825, 965, 695)
B = 100000

###(a)
mu0 = 1020
sd0 = sd(survivalTimes)
n = length(survivalTimes)

library(ggplot2)
df <- data.frame(survivalTimes)
p1 <- ggplot(df, aes(sample=survivalTimes)) +
    stat_qq( )+
    stat_qq_line()
    xlim(-2,2)
plot(p1)


shapiro.test(survivalTimes)$p.value

## Analysing the graphic and knowing that shapiro test p-value is greater than 0.05 we can affirm that the survival time estimates follow a normal Distribution.

p.value.exact = 1-pnorm((mean(survivalTimes)-mu0)/(sd(survivalTimes)/sqrt(n))); p.value.exact

set.seed(777)
t.obs = (mean(survivalTimes)-mu0)/(sd0/sqrt(n))
t.star = numeric(B)
z=survivalTimes-mean(survivalTimes)+mu0 
for(i in 1:B){
    z.star    = sample(z,n,replace=T)
    sd.z.star = sd(z.star)
    t.star[i] = (mean(z.star)-mu0)/(sd.z.star/sqrt(n))
}
# decision on H0 based on the p.value
p.value <- sum(t.star>t.obs)/B; p.value

## The exact p.value is greater than the p.value of the boostrap yet there are no significant changes. Both the exact p.value and the bootstrap p.value are less than 0.05 therefore we reject the null hypothesis.

###(b)
set.seed(777)
alpha = 0.1
t.star = numeric(B)
t= mean(survivalTimes)
for(i in 1:B){
    x.boot = sample(survivalTimes, length(survivalTimes), replace=T)
    t.star[i] = mean(x.boot)
}

### pivotal CI (review)
delta.star = t.star - t
d = quantile(delta.star, c(alpha/2,1-alpha/2))
ci.boot = t - c(d[2],d[1])
names(ci.boot) <- c("5%", "95%")
ci.boot

## draw histogram
library(ggplot2)
p1 <- ggplot(data.frame(bootstrap = t.star), aes(x = bootstrap)) +
      geom_histogram(aes(y = after_stat(density)))
plot(p1)

## percentile 90% CI
d = quantile(t.star, c(alpha/2,1-alpha/2)); d

##(c)
set.seed(777)
u <- c(alpha/2, 1-alpha/2) #desired quantiles
z <- qnorm(mean(t.star < t)) 
z.u <- qnorm(u)

## calculate accelaration using an estimator
t.star.jack = numeric(length(survivalTimes))
for(i in 1:length(survivalTimes)){
    t.star.jack[i] = mean(survivalTimes[-i])
}
mean.t.star.jack <- mean(t.star.jack)
a.hat <- sum((mean.t.star.jack - t.star.jack)^3)/(6*(sum((mean.t.star.jack - t.star.jack)^2))^(3/2))

#Adjusted quantiles
u_adjusted <- pnorm(z + (z+z.u)/(1-a.hat*(z+z.u))) 

#Accelerated Bootstrap CI
quantile(t.star, u_adjusted)

## (d)
set.seed(777)
library(boot)
boot.T <- function(data,indices){
    return(mean(data[indices,]))
}
boot.mean <- boot(data=as.data.frame(survivalTimes),statistic = boot.T,R=B)
boot.ci(boot.mean,type=c("basic","norm","perc", "bca"), conf=0.9)

## (e)
P.mean = mean(survivalTimes > 1100); P.mean

## (f)
### bootstrap
set.seed(777)
B= 10000
t.star.boot=numeric(B)
for(i in 1:B){
    survivalTimes.sample = sample(survivalTimes,length(survivalTimes),replace =T)
    t.star.boot[i] = mean(survivalTimes.sample > 1100)
}

t.star.boot.mean =mean(t.star.boot); t.star.boot.mean
t.star.boot.bias = t.star.boot.mean - P.mean; t.star.boot.bias
t.star.boot.var = var(t.star.boot); t.star.boot.var
t.star.boot.sd = sqrt(t.star.boot.var); t.star.boot.sd

### jackknife
set.seed(777)
t.star.jack=numeric(length(survivalTimes))
for(i in 1:length(survivalTimes)){
    t.star.jack[i] = mean(survivalTimes[-i] > 1100)
}

t.star.jack.mean = mean(t.star.jack); t.star.jack.mean
t.star.jack.bias = (t.star.jack.mean - P.mean) * (n-1); t.star.jack.bias
t.star.jack.var = mean((t.star.jack - t.star.jack.mean)^2) * (n-1); t.star.jack.var
t.star.jack.sd = sqrt(t.star.jack.var); t.star.jack.sd;

### (g)
set.seed(777)
B= 10000
t.star.boot=numeric(B)
t.star.jack=numeric(length(survivalTimes))
for(i in 1:B){
    z.star = sample(survivalTimes,length(survivalTimes),replace =T)
    t.star.boot[i] = mean(z.star > 1100)
}
for(i in 1:B){
    t.star.jack[i] = sd(t.star.boot[-i])
}

t.star.jack.var = mean((t.star.jack - mean(t.star.jack))^2) * (n-1); t.star.jack.var
t.star.jack.sd = sqrt(t.star.jack.var); t.star.jack.sd;
