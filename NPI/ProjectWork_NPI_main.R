## THIS FILE CONTAINS CODE FOR GENERATING RANDOM SAMPLES AND USEFUL FUNCTIONS ONLY ,
  ## PLOTS & COMPARISONS WILL BE IN SEPARATE FILE.

## Can obtain sampling distribution of various nonparametric LMP test statistic, can perform randomized test of exact size alpha.


#### libs----------------------------------------
library(SuppDists)
library(nimble)
library(latticeExtra)
library(car)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggiraph)
library(pracma)
library(MASS)
library(fMultivar)




#### Defining various Test Statitics-------------------------------

TH<-function(sample1 , sample2){ ## Terry-Hoeffding / Fisher-Yates test
  rank_comb = rank(c(sample1, sample2), ties.method = "random")
        # theoritically in a sample from continuous distribution ,there is no tie 
          # if there is tie in sample ,randomly allocate integer rank
  n1=length(sample1)
  n2=length(sample2)
  N_sum = n1+n2
  Rx = rank_comb[1:n1]
  
  statistic = sum((normOrder(N_sum))[Rx])
  return(statistic)
}

VdW<-function(sample1 , sample2){  ## Van der Waerden stat
  rank_comb = rank(c(sample1, sample2),ties.method = "random")
  n1=length(sample1)
  n2=length(sample2)
  N_sum = n1+n2
  Rx = rank_comb[1:n1]
  
  statistic = sum(qnorm(Rx/(N_sum+1)))
  return(statistic)
}

Med<-function(sample1 , sample2){  ## Median test
  rank_comb = rank(c(sample1, sample2), ties.method = "random")
  n1=length(sample1)
  n2=length(sample2)
  N_sum = n1+n2
  Rx = rank_comb[1:n1]
  
  statistic = 0.5*(sum(sign(Rx - 0.5*(N_sum+1))) + n1)
  return(statistic)
}

Wilcox<-function(sample1 , sample2){  ## Wilcoxon Rank-Sum test
  rank_comb = rank(c(sample1, sample2), ties.method = "random")
  n1=length(sample1)
  n2=length(sample2)
  N_sum = n1+n2
  Rx = rank_comb[1:n1]
  
  statistic = sum(Rx)/(N_sum+1)
  return(statistic)
}


t_two<-function(sample1 , sample2){   ## two sample t statistic
  n1 = length(sample1)
  n2 = length(sample2)
  N_sum = n1+n2
  
  pooled_v = ((n1-1)*var(sample1) + (n2-1)*var(sample2))/(N_sum -2)
  statistic = (mean(sample1) - mean(sample2))/sqrt(pooled_v*sum(1/c(n1,n2)))
  
}


C_LMPR<-function(sample1 , sample2){  ## LMPR stat for cauchy distribution 
  rank_comb = rank(c(sample1, sample2),ties.method = "random" )
  n1=length(sample1)
  n2=length(sample2)
  N_sum = n1+n2
  Rx = rank_comb[1:n1]
  
  statistic = sapply(Rx, 
                 function(r) {
                   integral(function(x) (2*x/(x^2 + 1))*(N_sum*choose(N_sum-1 ,r-1)*pow(pcauchy(x),r-1)*dcauchy(x)*pow(pcauchy(x,lower.tail = F),N_sum-r)),
                            -10000,10000) 
                   #(-Inf,Inf) is replaced by very large interval , so probability of that interval almost 1
           
          }) %>% sum
  
  return(statistic)
}


logE_LMPR<-function(sample1 , sample2){ ## LMPR for log(Exp) ,location family
  
  rank_comb = rank(c(sample1, sample2),ties.method = "random" )
  n1=length(sample1)
  n2=length(sample2)
  N_sum = n1+n2
  Rx = rank_comb[1:n1]
  
  g<-function(x) exp(x-exp(x)) # PDF for log(Exp)
  G<-function(x) 1-exp(-exp(x)) # CDF for log(Exp)
  
  
  statistic = sapply(Rx, 
                     function(r) {
                       integral(function(x) exp(x)*(N_sum*choose(N_sum-1 ,r-1)*pow(G(x),r-1)*g(x)*pow(1-G(x),N_sum-r)),
                                -100,100)
                       
                     }) %>% sum
  
  return(statistic)
}

paired_t<-function(sample1 , sample2){ ## paired sample t test
  d = sample1 - sample2
  n = length(d)
  N_sum = 2*n
  
  statistic = sqrt(n)*mean(d)/sqrt(var(d))
  
  return(statistic)
}


#### Generating Sampling Distribution of Rank Statistics from required distribution-------------------------------------------

generate<-function(n, rand_samp1 = rnorm, m=n, rand_samp2 = rand_samp1, Names = c(TH,VdW,Med), I=1000){ 
  ## n : size of sample1 , m : size of sample2 (by default m=n)
  ## rand_samp1 or rand_samp2 (by default they are same) : expression for random sample from given distribution ; 
  #e.g.-  function(size) rnorm(size , 2 , 1) etc.
  ## Names : Name of the Rank statistic required out of Terry-Hoeffding test, Van der Waerden test, Median test and Wilcoxon Rank-Sum test etc. defined above
  # by default Names = c(TH,VdW,Med) , mention any specific subset or other statistic if needed
  ## I : No. of iterations (by default 1000 ; if takes too much time , reduce it ; if permits increase it)
  
  
  X = matrix(rand_samp1(n*I),nrow = I) # taking random sample of size n for I times
  # is same as taking n*I random samples
  # for easier thinking we can put byrow = T , but has no effect in pattern
  Y = matrix(rand_samp2(m*I),nrow = I)
  SAMP = cbind(X,Y)
  Rank_Stats = sapply(Names, function(Rank_stat) {
    SAMP %>% 
      apply(MARGIN = 1, 
            FUN = function(row_i) Rank_stat(row_i[1:n],row_i[(n+1):(n+m)]))
  
    }  # in each row of SAMP : first n are sample1 ,remaining m are sample2 ; apply Rank_stat on them
  
  )
  return(data.frame(Rank_Stats)) 
  # each returned column corresponds to I obs. of a Rank_stat
}



generate2<-function(n, rand_biv_samp = function(s) mvrnorm(s,c(0,0),matrix(c(1,0.5,0.5,1),nrow = 2)), Names = c(TH,VdW,Med), I=1000){ 
  ## Bivariate extension of generate()
  
  SAMP = sapply(1:I, function(i) rand_biv_samp(n) %>% as.numeric() ) %>% t()
        # rand_biv_samp(n) is a n*2 matrix , each column has n obs
          # we are making in 2n length vector through as.numeric , first n component component corresponds X, last n corresponds to Y
            # SAMP contains I replications of such paired sample
  Rank_Stats = sapply(Names, function(Rank_stat) {
    SAMP %>% 
      apply(MARGIN = 1, 
            FUN = function(row_i) Rank_stat(row_i[1:n],row_i[(n+1):(2*n)]))
    
  }  
  
  )
  return(data.frame(Rank_Stats)) 
  
}





#### exact distribution for small n,m  -----------------------------------
Exact<-function(n = 4,m = 5,nam = c(TH)){
  N = m+n
  X_ranks = combn(N,n) # for large  n, m this step can take too much time to give all possible rank combinations of X
  exact_distn = 
  sapply(1:ncol(X_ranks),function(j) {
    fun = nam[[1]]
    fun(X_ranks[,j],(1:N)[-c(X_ranks[,j])])
                      ## if we treat X_ranks[,j] as X obs. and remaining integers from 1:N ,i.e.(1:N)[-c(X_ranks[,j])] as Y obs.
                            # based on that combined sample , ranks of X obs. are X_ranks[,j] itself 
    
  })
  return((table(exact_distn))/ncol(X_ranks))
}









#### Cut Off values for rejection region ------------------------------------------------
  # [NOTE: based on our setup H_0 will be rejected for small observed value of statistic] 

Cutoff_RankStat<-function(n,m = n,alpha = 0.05,Names = c(TH,VdW,Med,Wilcox),Itr = 10000){
 
  
  ##[HINT: simulate large number of values (by default 10000) of the statistic under H_0 to mimic its sampling distribution
    # take the lower alpha-th quantile as the cutoff for rejection region ]
  
  cutoff = (generate(n,function(s) rnorm(s),m,Names = Names,I = Itr)) %>% 
    # (as distribution free under H_0 ,for simplicity start with N(0,1))
                apply(MARGIN = 2, 
                      FUN = function(x){
                        eCDF = as.numeric(cumsum(table(x))/Itr)
                        id = sum(eCDF <= alpha)
                         ## to make our tests conservative ,
                          ## we are taking (j)-th obs. as cutoff when eCDF(j)<= alpha<eCDF(j+1)   
                            ## this function will fail only when Itr is so small that for 1st obs itself eCDF > alpha
                        mass = as.numeric(names(table(x))) # all mass points
                        
                        return(c(
                          mass[id], #cutoff
                          mass[id+1], # immediate next point of cutoff
                          (alpha - eCDF[id])/(eCDF[id+1] - eCDF[id]) # value of randomization parameter based on difference in level & size
                          ))
                        
                      })
  
    
    
  return(cutoff)
  ## output is a 3*length(Names) matrix
}


Cutoff_RankStat2<-function(n,rbiv_dist,alpha = 0.05,Names = c(TH,VdW,Med,Wilcox),Itr = 10000){ 
  ## Bivariate extension of cutoff()
  
  cutoff = (generate2(n,function(s) rbiv_dist(s),Names = Names,I = Itr)) %>% 
    apply(MARGIN = 2, 
          FUN = function(x){
            eCDF = as.numeric(cumsum(table(x))/Itr)
            id = sum(eCDF <= alpha)
            mass = as.numeric(names(table(x))) 
            
            return(c(
              mass[id], 
              mass[id+1], 
              (alpha - eCDF[id])/(eCDF[id+1] - eCDF[id]) 
            ))
            
          })
   return(cutoff)
  
}









#### POWER CURVES ---------------------------------------------


PowerComparison <- function(n,m = n,rdist = rnorm,Del_seq = (pow(1.1,1:10)-1.1)/(pow(1.1,10)),stat_names = c(1:4,6) ,alpha = 0.05,repli = 5000 , sensi = 10000, randomiz = F , zooming = F ,plotting = T ,return_data = !plotting){
  ## n , m : sample sizes
  ## rdist : expression needed to generate random data from required distribution (with proper location & scale ,if any)
    #  e.g. - rnorm (by default) , function(s) rdexp(s,-2,2) etc.
  ## Del_seq : sequence of Delta values to evaluate power curve ( by default 10 points between 0 & 1 , with GP increament)
  ## stat_names : name of the statistics to be compared ,
    #put TH:1 ,VdW:2 ,Med:3 ,Wilcox:4 ,C_LMPR:5 ,t_two:6 ,logE_LMPR:7
  ## alpha : level of significance ( 5% by default )
  ## repli : No. of replication (by default 5000) to estimate power , i.e estimated_power = (No. of times the statistic is in rejection region)/repli
  ## sensi : No. of iterations to calculate cut-off values (by default 10000) 
  ## randomiz : whether or not to perform exact size alpha test by randomized testing rule
      # by default F ,put T for very small sample sizes
  ## plotting : whether power curve will be drawn (T by default)
  ## zooming : whether the plot is fixed or can be zoomed (F by default)
  ## return_data : whether the numerical value of powers will be returned (!plotting by default)
  
  
  ### naming ........................
  if (identical(rdist,rnorm)) { a="Normal" }
  else if (identical(rdist,rcauchy)) { a="Cauchy" }
  else if (identical(rdist,rlogis)) { a="Logistic" }
  else if (identical(rdist,rdexp)) { a="Laplace" }
  else { a=" " }
  all = c(TH,VdW,Med,Wilcox,C_LMPR,t_two,logE_LMPR)
  all_names = c("delta","FisherY","VDW","Median","WRS","Cauchy","t_test","logE_LMPR")
  
  th_seq = Del_seq
  
  ### cutoffs .......................
  Cutoff = Cutoff_RankStat(n,m,alpha,Names = (all)[stat_names] , Itr = sensi)
      # a 3*length(stat_names) matrix ,for each col 3rows are respectively cutoff, next value , level-size
  
  
  cutoff = Cutoff[1,]
  nextval = Cutoff[2,]
  lam = Cutoff[3,]
  
  ### calculating power .........................................................
  
  power_data = matrix(th_seq) %>% apply(MARGIN = 1,FUN = function(del) generate(n,function(s) rdist(s),m,function(s) (del + rdist(s)), Names = (all)[stat_names], I = repli))
  # generating (data from dist.with location = del) is same as 
  # generating (del + (data from dist.with location = 0))

  if(!randomiz){power_values = matrix(1:length(th_seq)) %>% 
                 apply(MARGIN = 1,
                       FUN = function(i) { # for one particular value of theta 
                         (t(power_data[[i]]) <= cutoff) %>% apply(MARGIN = 1,FUN = mean)
                         ## counting how many times observed statistic is in rejection region
                            ## & estimating probability as relative frequency
                       }            
                 )
  }
  else{
    power_values = matrix(1:length(th_seq)) %>% 
      apply(MARGIN = 1,
            FUN = function(i) { 
            power_data[[i]] %>% 
                apply(MARGIN = 1 ,
                      function(x) {
                          ((x<=cutoff) +
                          (x==nextval)*(sapply(lam, function(p) rbinom(1,1,p)))) 
                      ## each position is now 1 ,if the obs. is in rejection region or a bernoulli(lam) experiment gives Success              
                      }) %>%  apply(MARGIN = 1,mean)
                  
                         
              }            
           )
    
    
  }
  ## power_values is a matrix of length(stat_names)*length(th_seq) now
  ## each col corresponds to a theta value and contains power of those tests at that theta
 

  POWER_DATA = data.frame(th_seq,t(power_values))
  colnames(POWER_DATA) <- all_names[c(1,(stat_names+1))]
  
  if(plotting) {
    ### plotting ..............................
    plot_data <- melt(POWER_DATA,id="delta")
    
    POWER_CURVE <- ggplot(plot_data,aes(x=delta,y=value,color=variable)) + geom_line() + 
      geom_line(aes(y = alpha), color = "black" , lty = 3) +
      scale_color_manual(values=(c("red","darkblue","darkgreen","darkorange","darkmagenta","purple"))[1:length(stat_names)]) +
      labs(title=paste(a,":n=",n,",m=",m ), y = "Power") + theme_bw() +
      theme(plot.title = element_text(colour="steelblue",size=15),axis.title.y = element_text(size=8)) + annotate("text", x= th_seq[length(th_seq)]/2, y=alpha, label= paste("level = ",alpha))
    
    if(zooming) {
      POWER_CURVE = girafe_options(girafe(ggobj = POWER_CURVE),
                                   opts_zoom( max = 5) )
    }
   
  }
  
  if(plotting*return_data) return(list(POWER_CURVE,POWER_DATA))
  else if(plotting) return(POWER_CURVE)
  else return(POWER_DATA)
 
}




PowerComparison2 <- function(n,rbiv_dist = function(s,del) mvrnorm(s,c(0,del),matrix(c(1,0.5,0.5,1),nrow = 2)),Del_seq = (pow(1.1,1:10)-1.1)/(pow(1.1,10)),stat_names = c(1:4,6) ,alpha = 0.05,repli = 5000 , sensi = 10000, randomiz = F , zooming = F ,plotting = T ,return_data = !plotting){
  ## Bivariate extension of PowerComparison()
  
  all = c(TH,VdW,Med,Wilcox,C_LMPR,t_two,logE_LMPR,paired_t)
  all_names = c("delta","FisherY","VDW","Median","WRS","Cauchy","t_test","logE_LMPR","paired_t")
  
  th_seq = Del_seq
  
  Cutoff = Cutoff_RankStat2(n,function(s) rbiv_dist(s,0),alpha,Names = (all)[stat_names] , Itr = sensi)
  cutoff = Cutoff[1,]
  nextval = Cutoff[2,]
  lam = Cutoff[3,]
  
  
  power_data = matrix(th_seq) %>% apply(MARGIN = 1,FUN = function(del) generate2(n,function(s) rbiv_dist(s,del), Names = (all)[stat_names], I = repli))
  
  
  if(!randomiz){power_values = matrix(1:length(th_seq)) %>% 
    apply(MARGIN = 1,
          FUN = function(i) { 
            (t(power_data[[i]]) <= cutoff) %>% apply(MARGIN = 1,FUN = mean)
            
          }            
    )
  }
  else{
    power_values = matrix(1:length(th_seq)) %>% 
      apply(MARGIN = 1,
            FUN = function(i) { 
              power_data[[i]] %>% 
                apply(MARGIN = 1 ,
                      function(x) {
                        ((x<=cutoff) +
                           (x==nextval)*(sapply(lam, function(p) rbinom(1,1,p)))) 
                      }) %>%  apply(MARGIN = 1,mean)
              
              
            }            
      )
    
    
  }
  
  POWER_DATA = data.frame(th_seq,t(power_values))
  colnames(POWER_DATA) <- all_names[c(1,(stat_names+1))]
  
  if(plotting) {
    plot_data <- melt(POWER_DATA,id="delta")
    
    POWER_CURVE <- ggplot(plot_data,aes(x=delta,y=value,color=variable)) + geom_line() + 
      geom_line(aes(y = alpha), color = "black" , lty = 3) +
      scale_color_manual(values=(c("red","darkblue","darkgreen","darkorange","darkmagenta","purple"))[1:length(stat_names)]) +
      labs(title=paste(":n=",n ), y = "Power") + theme_bw() +
      theme(plot.title = element_text(colour="steelblue",size=15),axis.title.y = element_text(size=8)) + annotate("text", x= th_seq[length(th_seq)]/2, y=alpha, label= paste("level = ",alpha))
    
    if(zooming) {
      POWER_CURVE = girafe_options(girafe(ggobj = POWER_CURVE),
                                   opts_zoom( max = 5) )
    }
    
  }
  
  if(plotting*return_data) return(list(POWER_CURVE,POWER_DATA))
  else if(plotting) return(POWER_CURVE)
  else return(POWER_DATA)
  
}

