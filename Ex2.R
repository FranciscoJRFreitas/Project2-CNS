#2.
  #2.1
    os = c(1.977866, 1.836622, 1.097168, 1.232889, 1.229526, 2.438342, 1.551389, 1.300618, 1.068584,
                1.183466, 2.179033, 1.535904, 1.323500, 1.458713, 1.013755, 3.602314, 1.087067, 1.014013,
                1.613929, 2.792161, 1.197081, 1.021430, 1.111531, 1.131036, 1.064926)
    
    n = length(os)
    
    #MLE
    MLE = function(x) {n/sum(log(x))}
    
    MLE.v = MLE(os)
    
    MLE.v
    
    #MME
    MME = function(x){mean(x)/(mean(x)-1)}
    
    MME.v = MME(os)
    
    MME.v
    
  #2.2
    In = n/(MLE.v^2)
    MLE.var = 1/In
    
    MLE.var
  #2.3
    #Likelihood function
    likh = function(theta){
      l = vector()
      for(t in theta)
        l = c(l, (t^n)*prod(1/(os^(t+1))))
      return(l)
    }
    
    #Log-Likelihood function
    loglikh = function(theta){
      ll = vector()
      for(t in theta)
        ll = c(ll, n*log(t)-((t + 1) * sum(log(os))))
      return(ll)
    }
    
    #Score function
    scr = function(theta) {
      s = vector()
      for(t in theta)
        s = c(s, (n/t) - sum(log(os)))
      return(s)
    }
    
    #Graphical Display
    
    library(ggplot2)
    base <-
      ggplot() +
      xlim(0, 4)
    
    base + xlab(expression(theta)) + ylab("likelihood") + geom_function(fun = likh)
    
    base + xlab(expression(theta)) + ylab("log-likelihood") + geom_function(fun = loglikh)
    
    base + xlab(expression(theta)) + ylab("score") + geom_function(fun = scr)
    
    
    #If we use interval estimation (2,4):
    
    #Likelihood function maximum
    optimize(likh,c(2,4),maximum=T)$maximum
    
    #Log-Likelihood function maximum
    optimize(loglikh,c(2,4),maximum=T)$maximum
    
    #Score function root
    uniroot(scr,c(2,4))$root
    
  #2.4
    #(install.packages("maxLik"))
    library(maxLik)
    maxLik(loglikh,start=2.81418)
    # Maximum Likelihood estimation
    # Newton-Raphson maximisation, 1 iterations
    # Return code 1: gradient close to zero (gradtol)
    # Log-Likelihood: -8.016816 (1 free parameter(s))
    # Estimate(s): 2.814179
    
    #The ML estimate of alpha by R function maxLik() was 2.814179.
    
  #2.5
    #Bissection algorithm to estimate alpha
    
    bisection <- function(a,b,eps){
      # [a,b]: interval where s verifies Bolzano’s theorem
      # eps : is the stopping rule
      alpha.it = vector(); alpha.it[1] = (a+b)/2
        k = 1; diff = 1
        diffs= vector()
        while(diff>eps){
          #If multiplication is negative, reduce b
          if(scr(alpha.it[k])*scr(a)<0){
            b = alpha.it[k]
            alpha.it[k+1] = (a+b)/2
          }
          #If multiplication is positive, raise a
          else{if(scr(alpha.it[k])*scr(a)>0){
            a = alpha.it[k]
            alpha.it[k+1] = (a+b)/2
            #If multiplication result is 0, we apply the stopping rule
          }else{alpha.it[k+1]=alpha.it[k]}
          }
          #Update the iteration difference for the stopping rule
          diff = abs(alpha.it[k+1]-alpha.it[k])
          diffs[k]=diff
          k = k+1
        }
        result = as.matrix(alpha.it)
        colnames(result)<-"iterations"
        rownames(result)<-1:length(alpha.it)
        retList<-list("res"=result,"error"=diffs)
        retList
      }
      #Bisection result for each iteration
      l1=bisection(2,4,1e-06)
      m1=t(l1$res)
      m1
      
      #Error result for each iteration
      l1$error
      
      #(Newton-Raphson aux function)
      prime <- function(alpha){
        out = numeric(length(alpha))
        if(length(alpha)==1){out =- n/alpha^2}
        if(length(alpha)!=1){
          for(i in 1:length(alpha)){out[i] = - n/alpha[i]^2}
        }
        return(out)
      }
      
    #Newton-Raphson algorithm to estimate alpha
     NR <- function(x,alpha0,eps){
       # x: observed sample
       # alpha0: first value to iterate over the vector
       # eps: the stopping rule
       alpha.it = vector()
       alpha.it[1] = alpha0
       k = 1
       diff = 1
       diffs=vector()
       #iterations
       while(diff>eps){
         alpha.it[k+1] = alpha.it[k]-scr(alpha.it[k])/prime(alpha.it[k])
         diff = abs(alpha.it[k+1]-alpha.it[k])
         diffs[k]=diff
         k = k+1
       }
       result = as.matrix(alpha.it)
       colnames(result)<-"iterations"
       rownames(result)<-1:length(alpha.it)
       result
       retList<-list("res"=result,"error"=diffs)
       retList
     }
     #NR result for each iteration
     l2=NR(os,MLE.var,1e-06)
     m2=t(l2$res)
     m2
      
     #Error result for each iteration
     l2$error
      
     
     #Secant algorithm to estimate alpha
     secant <- function(x, alpha0, alpha1, eps) {
       # x: observed sample
       # alpha0: interval left limit
       # alpha1: interval right limit
       # eps: stopping rule value
       alpha.it = vector()
       alpha.it[1] = alpha0
       32
       alpha.it[2] = alpha1
       k = 2
       diff = 1
       diffs = vector()
       diffs[1] = diff
       # iterations
       while (diff > eps) {
         alpha.it[k + 1] = alpha.it[k] - scr(alpha.it[k]) * 
           ((alpha.it[k] - alpha.it[k-1])/(scr(alpha.it[k]) - scr(alpha.it[k - 1])))
         diff = abs(alpha.it[k + 1] - alpha.it[k])
         diffs[k] = diff
         k = k + 1
       }
       result = as.matrix(alpha.it)
       colnames(result) <- "iterations"
       rownames(result) <- 1:length(alpha.it)
       retList <- list(res = result, error = diffs)
       retList
     }
     #Secant result for each iteration
     l3 = secant(os, 2, 4, 1e-06)
     m3 = t(l3$res)
     m3
     
     #Error result for each iteration
     l3$error
      
    
  #2.6
     #I(alpha)
     Ialpha<-function(alpha){
       n/(alpha^2)
     }
     #Fisher Scoring algorithm
     fisherScoring <- function(x,alpha0,eps){
       # x: observed sample
       # alpha0: first value to iterate over the vector
       # eps: the stopping rule
       alpha.it = vector()
       alpha.it[1] = alpha0
       k = 1
       diff = 1
       diffs=vector()
       diffs[1]=diff
       # iterations
       while(diff>eps){
         alpha.it[k+1] = alpha.it[k]+(1/Ialpha(alpha.it[k]))*scr(alpha.it[k])
         diff = abs(alpha.it[k+1]-alpha.it[k])
         k = k+1
         diffs[k]=diff
       }
       result = as.matrix(alpha.it)
       colnames(result)<-"iterations"
       rownames(result)<-1:length(alpha.it)
       retList<-list("res"=result,"error"=diffs)
       retList
     }
     #Fisher Scoring result for each iteration
     l4=fisherScoring(os,MLE.var,0.000001)
     m4=t(l4$res)
     m4
     
     #Error result for each iteration
     l4$error
