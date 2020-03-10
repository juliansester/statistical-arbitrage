# Ar_1 Function for simulations
AR_1<-function(t=0.1,N=10000,mu=100,k=0.1,sigma=20,plot=T,plot_threshold=2*sigma,S_0=100){
  norm_rv<-rnorm(N,mean=0,sd=t*sigma)
  #  norm_rv<-rep(1,N)
  Y<-numeric(N)
  Y[1]<-S_0
  for(i in 2:N){
    Y[i]<-Y[i-1]+k*(mu-Y[i-1])*t+norm_rv[i]
  }
  # Y<-arima.sim(list( ar=(1-k*t),innov=(norm_rv+rep(k*mu*t,N))),n.start = 1,start.innov = mu,n=N)[1:N]
  if(plot){
    plot(1:N,Y,type='l',col="black",lwd=0.5,xlab="Time",ylab="Spread")
    abline(h=mu+plot_threshold,lty=2)
    abline(h=mu,lty=2)
    abline(h=mu-plot_threshold,lty=2)
  }
  return(Y)
  
}





N<-5 # different Parameters
n<-100 # time steps 
n_discr<-50
time_diff<-1
kk<-seq(0.005,0.1,length=N)
# Setting Parameters
mu<-seq(0,0,length=N)
sigma<-seq(1,1,length=N)
S_0<-0
## Values
x<-matrix(0,N,n_discr)
y<-matrix(0,N,n_discr)
standard_dev<-0
for(j in 0:(n/2)){
  standard_dev<-standard_dev+((1-kk*time_diff)^j)*sigma*time_diff
}
standard_dev<-sqrt(standard_dev)
## randomize
kk<-sample(kk,N)
mu<-sample(mu,N)
sigma<-sample(sigma,N)

# probabilities
prob_set<-array(0,dim=c(N,n_discr,n_discr))
for(i in 1:N){
  # Assigning Values
  x[i,]<-seq(S_0-50,S_0+50,length=n_discr)
  y[i,]<-seq(S_0-50,S_0+50,length=n_discr)
  for(j in 1:n_discr){
    for(k in 1:n_discr){
      # Assigning Probabilities
      prob_set[i,j,k]<-dnorm(x[i,j],mean=mu[i]*(1-(1-kk[i]*time_diff)^(n/2))+((1-kk[i]*time_diff)^(n/2))*S_0,sd=standard_dev[i])*dnorm(y[i,k],mean=mu[i]*(1-(1-kk[i]*time_diff)^(n/2))+((1-kk[i]*time_diff)^(n/2))*x[i,j],sd=standard_dev[i])
    }
  }
  prob_set[i,,]<-prob_set[i,,]/sum(prob_set[i,,])
}

### Removing too low probabilities
epsilon<-0.0001
for(i in 1:N){
  for(j in 1:n_discr){
    for(k in 1:n_discr){
      if((prob_set[i,j,k])<epsilon){
        prob_set[i,j,k]<-0
      }
    }
  }
  prob_set[i,,]<-prob_set[i,,]/sum(prob_set[i,,])
}

persp(x[1,],y[1,],prob_set[1,,],xlab="S_50",ylab="S_100",zlab="Prob.",theta=30)
AR_1(t=time_diff,N=n+1,mu=mu[N],k=kk[N],sigma=sigma[N],plot=T,plot_threshold=2*sigma[N],S_0=S_0)
# Function Definition

func1<-function(x,y){0}

# Computation of the Bounds
limit<-10

# Robust Bounds
check_statarb(x[1,],y[1,],prob_set[1,,],S_0)
super_rep_fixp(x[1,],y[1,],prob_set[1,,],func=function(x,y){0},stat_arb=T,lower=F,S_0,limit)
lower_bound_nosa<-super_rep(x,y,prob_set,func=func1,stat_arb=F,lower=T,S_0,limit)
lower_bound<-super_rep_robust(x,y,prob_set,func=func1,stat_arb=T,lower=T,S_0,limit)
upper_bound<-super_rep_robust(x,y,prob_set,func=func1,stat_arb=T,lower=F,S_0,limit)
upper_bound_nosa<-super_rep(x,y,prob_set,func=func1,stat_arb=F,lower=F,S_0,limit)
# for(loop_ii in 1:10){
#   kk<-sample(kk,N)
#   mu<-sample(mu,N)
#   sigma<-sample(sigma,N)
#   prob_set<-array(0,dim=c(N,n_discr,n_discr))
#   for(i in 1:N){
#     # Assigning Values
#     x[i,]<-seq(S_0-50,S_0+50,length=n_discr)
#     y[i,]<-seq(S_0-50,S_0+50,length=n_discr)
#     for(j in 1:n_discr){
#       for(k in 1:n_discr){
#         # Assigning Probabilities
#         prob_set[i,j,k]<-dnorm(x[i,j],mean=mu[i]*(1-(1-kk[i]*time_diff)^(n/2))+((1-kk[i]*time_diff)^(n/2))*S_0,sd=standard_dev[i])*dnorm(y[i,k],mean=mu[i]*(1-(1-kk[i]*time_diff)^(n/2))+((1-kk[i]*time_diff)^(n/2))*x[i,j],sd=standard_dev[i])
#       }
#     }
#     prob_set[i,,]<-prob_set[i,,]/sum(prob_set[i,,])
#   }
#   
#   ### Removing too low probabilities
#   epsilon<-0.0001
#   for(i in 1:N){
#     for(j in 1:n_discr){
#       for(k in 1:n_discr){
#         if((prob_set[i,j,k])<epsilon){
#           prob_set[i,j,k]<-0
#         }
#       }
#     }
#     prob_set[i,,]<-prob_set[i,,]/sum(prob_set[i,,])
#   }
#   lower_bound_nosa_1<-super_rep(x,y,prob_set,func=func1,stat_arb=F,lower=T,S_0,limit)
#   if(lower_bound_nosa_1$d<lower_bound_nosa$d){lower_bound_nosa<-lower_bound_nosa_1}
#   lower_bound_1<-super_rep(x,y,prob_set,func=func1,stat_arb=T,lower=T,S_0,limit)
#   if(lower_bound_1$d<lower_bound$d){lower_bound<-lower_bound_1}
#   upper_bound_1<-super_rep(x,y,prob_set,func=func1,stat_arb=T,lower=F,S_0,limit)
#   if(upper_bound_1$d>upper_bound$d){upper_bound<-upper_bound_1}
#   upper_bound_nosa_1<-super_rep(x,y,prob_set,func=func1,stat_arb=F,lower=F,S_0,limit)
#   if(upper_bound_nosa_1$d>upper_bound_nosa$d){upper_bound_nosa<-upper_bound_nosa_1}
# }


# Simulations
## Statistical Arbitrage
par(mfrow=c(1,2))
par(mar=c(2.5, 4.5,1.5, 1.5))

Delta_0_upper<-upper_bound$Delta_0
kk<-sample(kk,N)
mu<-sample(mu,N)
sigma<-sample(sigma,N)
Delta_1_upper<-splinefun(upper_bound$x,upper_bound$Delta_1,method="natural")
#plot(upper_bound$x,Delta_1_upper(upper_bound$x),type="l")
nr_gains<-0
nr_losses<-0
d_upper<-upper_bound$d
gain<-0
Nr_sim<-1000000
best<-0
worst<-0
ff<-numeric(Nr_sim)
for(i in 1:Nr_sim){
  ii<-sample(1:N,1) # choosing random parameter
  jj<-sample(1:N,1) 
  ll<-sample(1:N,1)
  S<-AR_1(t=time_diff,N=n+1,mu=mu[ii],k=kk[jj],sigma=sigma[ll],plot=F,plot_threshold=2*sigma,S_0=S_0)
  f<-(Delta_0_upper*(S[n/2]-S_0)+Delta_1_upper(S[n/2])*(S[n]-S[n/2]))
  ff[i]<-f
  gain<-gain+f
  if(f<0){nr_losses<-nr_losses+1}
  if(f>0){nr_gains<-nr_gains+1}
  if(f>best){best<-f}
  if(f<worst){worst<-f}
}
gain/Nr_sim
best
worst
nr_losses/Nr_sim
nr_gains/Nr_sim
hist(ff,xlab="Average Gain",main="Statistical Arbitrage with K=10",breaks=50)

########## Delta Strategies ###############
################# PAIRS TRADING ###########

Delta_strat<-function(N,eps=20,delta=0.001,Y,mu=0,intermediate_NR=N-1){
  ## Define Trading Points
  trading_times<-round(seq(N/(intermediate_NR+1),N*intermediate_NR/(intermediate_NR+1),length=intermediate_NR)+1,0)
  Delta<-numeric(N)
  for(i in 2:N){
    if((Y[i]>=(mu+eps))*(Delta[i-1]==0)*(i %in% trading_times)){
      Delta[i]<-(Delta [i-1]-1)
    }
    else if((Y[i]>=(mu-delta))*(Delta[i-1]==1)*(i %in% trading_times)){
      Delta[i]<-(Delta [i-1]-1)
    }
    else if((Y[i]<=(mu-eps))*(Delta[i-1]==0)*(i %in% trading_times)){
      Delta[i]<-(Delta [i-1]+1)
    }
    else if((Y[i]<=(mu+delta))*(Delta[i-1]==(-1))*(i %in% trading_times)){
      Delta[i]<-(Delta [i-1]+1)
    }
    else{
      Delta[i]<-(Delta [i-1])
    }
  }
  return(Delta)
}
Gain<-function(Delta,N,S){
  G<-numeric(N)
  for(i in 2:N){
    G[i]<-Delta[i-1]*(S[i]-S[i-1])
  }
  G<-cumsum(G)
  # print(G[N])
  return(G)
}
# par(mfrow=c(2,1))
# par(mar=c(2.5, 4.5,1.5, 1.5))
# set.seed(103) # Good Scenario
# #set.seed(434) # Bad Scenario
# Y<-AR_1(time_diff,N=n+1,mu[1],k=kk[1],sigma[1],plot=T,plot_threshold = 2*sigma,S_0=S_0)
# #D1<-Delta_strat(N,eps=2*sigma,delta=0.001,Y)
# D1<-Delta_strat(N=n+1,eps=2*sigma[1],delta=0.001,Y,mu[1],intermediate_NR=99)
# 
# # Trading-Points einzeichnen
# points((1:(n+1))[c(F,diff(D1)!=0)],Y[c(F,diff(D1)!=0)],pch=19,cex=1.5,col="black")
# 
# # Anzahl der Long/shortPositionen bestimmen
# long_entry<-(1:(n+1))[c(F,(diff(D1)==1)&(D1[-1]!=0))]
# long_exit<-(1:(n+1))[c(F,(diff(D1)==-1)&(D1[-(n+1)]==1))]
# short_entry<-(1:(n+1))[c(F,(diff(D1)==-1)&(D1[-1]!=0))]
# short_exit<-(1:(n+1))[c(F,(diff(D1)==1)&(D1[-(n+1)]==-1))]
# nlong<-length(long_entry)
# nshort<-length(short_entry)
# if(nlong>length(long_exit)){long_exit<-c(long_exit,(n+1))}
# if(nshort>length(short_exit)){short_exit<-c(short_exit,(n+1))}
# # Wenn investiert, dann kolorieren
# for(i in 1:nlong){
#   lines((1:(n+1))[seq(long_entry[i],long_exit[i],by=1)],Y[seq(long_entry[i],long_exit[i],by=1)],lty=1,lwd=0.5,col="red")
# }
# for(i in 1:nshort){
#   lines((1:(n+1))[seq(short_entry[i],short_exit[i],by=1)],Y[seq(short_entry[i],short_exit[i],by=1)],lty=1,lwd=0.5,col="blue")
# }
# #Gain plotten
# plot(1:(n+1),Gain(D1,(n+1),Y),type="l",xlab="Time",ylab="Gain")
# for(i in 1:nlong){
#   lines((1:(n+1))[seq(long_entry[i],long_exit[i],by=1)],Gain(D1,(n+1),Y)[seq(long_entry[i],long_exit[i],by=1)],lty=1,lwd=0.5,col="red")
# }
# for(i in 1:nshort){
#   lines((1:(n+1))[seq(short_entry[i],short_exit[i],by=1)],Gain(D1,(n+1),Y)[seq(short_entry[i],short_exit[i],by=1)],lty=1,lwd=0.5,col="blue")
# }

##############Simulations
#epsilon=2*sigma

nr_gains<-0
nr_losses<-0
gain<-0
Nr_sim<-1000000
best<-0
worst<-0
fff<-numeric(Nr_sim)
for(sim in 1:Nr_sim){
  ii<-sample(1:N,1) # choosing random parameter
  jj<-sample(1:N,1) 
  ll<-sample(1:N,1)
  S<-AR_1(t=time_diff,N=n+1,mu=mu[ii],k=kk[jj],sigma=sigma[ll],plot=F,plot_threshold=2*sigma,S_0=S_0)
  D1<-Delta_strat(N=n+1,eps=2*sigma[ii],delta=0.001,S,mu[ii])
  f<-Gain(D1,(n+1),S)[n+1]
  fff[sim]<-f
  gain<-gain+f
  if(f<0){nr_losses<-nr_losses+1}
  if(f>=0){nr_gains<-nr_gains+1}
  if(f>best){best<-f}
  if(f<worst){worst<-f}
}
round(gain/Nr_sim,2)
round(best,2)
round(worst,2)
round(100*nr_losses/Nr_sim,2)
round(100*nr_gains/Nr_sim,2)
hist(fff,xlab="Average Gain",main="Pairs Trading with Epsilon = 2sigma",breaks=50)

################ More trading points #########################


more_trading<-function(Nr_sim=1000,nr_trading_points=1,bad_scenarios=F){

f<-0
nr_gains<-0
nr_losses<-0
gain<-0
best<-0
worst<-0
standard_dev<-0
for(j in 0:(n/(nr_trading_points+1))){
  standard_dev<-standard_dev+((1-kk*time_diff)^j)*sigma*time_diff
}
standard_dev<-sqrt(standard_dev)

if(bad_scenarios){
  ff<-NULL
  for(s in bad_seeds){
    set.seed(s)
    S<-AR_1(time_diff,N=n+1,mu[1],k=kk[1],sigma[1],plot=F,plot_threshold = 2*sigma,S_0=S_0)
  for(h in 1:nr_trading_points){
    ## Values
    x<-matrix(0,N,n_discr)
    y<-matrix(0,N,n_discr)
    
    # probabilities
    prob_set<-array(0,dim=c(N,n_discr,n_discr))
    for(i in 1:N){
      # Assigning Values
      x[i,]<-seq(S_0-50/nr_trading_points,S_0+50/nr_trading_points,length=n_discr)
      y[i,]<-seq(S_0-50/nr_trading_points,S_0+50/nr_trading_points,length=n_discr)
      for(j in 1:n_discr){
        for(k in 1:n_discr){
          # Assigning Probabilities
          prob_set[i,j,k]<-dnorm(x[i,j],mean=mu[i]*(1-(1-kk[i]*time_diff)^((n/(nr_trading_points+1))))+((1-kk[i]*time_diff)^((n/(nr_trading_points+1))))*S_0,sd=standard_dev[i])*dnorm(y[i,k],mean=mu[i]*(1-(1-kk[i]*time_diff)^(n/(nr_trading_points+1)))+((1-kk[i]*time_diff)^(n/(nr_trading_points+1)))*x[i,j],sd=standard_dev[i])
        }
      }
      prob_set[i,,]<-prob_set[i,,]/sum(prob_set[i,,])
    }
    
    ### Removing too low probabilities
    # epsilon<-0.0001
    # for(i in 1:N){
    #   for(j in 1:n_discr){
    #     for(k in 1:n_discr){
    #       if((prob_set[i,j,k])<epsilon){
    #         prob_set[i,j,k]<-0
    #       }
    #     }
    #   }
    #   prob_set[i,,]<-prob_set[i,,]/sum(prob_set[i,,])
    # }
    
    # Function Definition
    
    
    
    # Robust Bounds
    upper_bound<-super_rep_robust(x,y,prob_set,func=function(x,y){0},stat_arb=T,lower=F,S_0,limit=10)
    Delta_0_upper<-upper_bound$Delta_0
    Delta_1_upper<-splinefun(upper_bound$x,upper_bound$Delta_1,method="natural")
    
    
    
    f<-f+(Delta_0_upper*(S[(h*n)/(nr_trading_points+1)]-S_0)+Delta_1_upper(S[((h)*n)/(nr_trading_points+1)])*(S[((h+1)*n)/(nr_trading_points+1)]-S[(h*n)/(nr_trading_points+1)]))
    S_0<-S[(h*n)/(nr_trading_points+1)]

  }
  ff<-c(ff,f)
  gain<-gain+f
  if(f<0){nr_losses<-nr_losses+1}
  if(f>0){nr_gains<-nr_gains+1}
  if(f>best){best<-f}
  if(f<worst){worst<-f}
  f<-0
  S_0<-0
  
  }
  hist(ff,xlab="Average Gain",main=paste("Stat. Arb with ",  nr_trading_points, " Trades", sep=" "),breaks=50)
  return(list(
    AVERAGE_GAIN=round(gain/length(bad_seeds),2),
    BEST=round(best,2),
    WORST=round(worst,2),
    LOSSES=round(100*nr_losses/length(bad_seeds),2),
    GAINS=round(100*nr_gains/length(bad_seeds),2)
  )
  )
}
else{
  ff<-numeric(Nr_sim)

for(sim in 1:Nr_sim){
  ii<-sample(1:N,1) # choosing random parameter
  S<-AR_1(t=time_diff,N=n,mu=mu[ii],k=kk[ii],sigma=sigma[ii],plot=F,plot_threshold=2*sigma,S_0=S_0)
  for(h in 1:nr_trading_points){
    ## Values
    x<-matrix(0,N,n_discr)
    y<-matrix(0,N,n_discr)
    
    # probabilities
    prob_set<-array(0,dim=c(N,n_discr,n_discr))
    for(i in 1:N){
      # Assigning Values
      x[i,]<-seq(S_0-50/nr_trading_points,S_0+50/nr_trading_points,length=n_discr)
      y[i,]<-seq(S_0-50/nr_trading_points,S_0+50/nr_trading_points,length=n_discr)
      for(j in 1:n_discr){
        for(k in 1:n_discr){
          # Assigning Probabilities
          prob_set[i,j,k]<-dnorm(x[i,j],mean=mu[i]*(1-(1-kk[i]*time_diff)^((n/(nr_trading_points+1))))+((1-kk[i]*time_diff)^((n/(nr_trading_points+1))))*S_0,sd=standard_dev[i])*dnorm(y[i,k],mean=mu[i]*(1-(1-kk[i]*time_diff)^(n/(nr_trading_points+1)))+((1-kk[i]*time_diff)^(n/(nr_trading_points+1)))*x[i,j],sd=standard_dev[i])
        }
      }
      prob_set[i,,]<-prob_set[i,,]/sum(prob_set[i,,])
    }
    
    ### Removing too low probabilities
    # epsilon<-0.0001
    # for(i in 1:N){
    #   for(j in 1:n_discr){
    #     for(k in 1:n_discr){
    #       if((prob_set[i,j,k])<epsilon){
    #         prob_set[i,j,k]<-0
    #       }
    #     }
    #   }
    #   prob_set[i,,]<-prob_set[i,,]/sum(prob_set[i,,])
    # }
    
    # Function Definition
    
    
    
    # Robust Bounds
    upper_bound<-super_rep_robust(x,y,prob_set,func=function(x,y){0},stat_arb=T,lower=F,S_0,limit=10)
    Delta_0_upper<-upper_bound$Delta_0
    Delta_1_upper<-splinefun(upper_bound$x,upper_bound$Delta_1,method="natural")
    
    
    
    f<-f+(Delta_0_upper*(S[(h*n)/(nr_trading_points+1)]-S_0)+Delta_1_upper(S[((h)*n)/(nr_trading_points+1)])*(S[((h+1)*n)/(nr_trading_points+1)]-S[(h*n)/(nr_trading_points+1)]))
    S_0<-S[(h*n)/(nr_trading_points+1)]
    ff[sim]<-f
  }
  gain<-gain+f
  if(f<0){nr_losses<-nr_losses+1}
  if(f>0){nr_gains<-nr_gains+1}
  if(f>best){best<-f}
  if(f<worst){worst<-f}
  f<-0
  S_0<-0
}
}
hist(ff,xlab="Average Gain",main=paste("Stat. Arb with ",  nr_trading_points, " Trades", sep=" "),breaks=50)
return(list(
AVERAGE_GAIN=round(gain/Nr_sim,2),
BEST=round(best,2),
WORST=round(worst,2),
LOSSES=round(100*nr_losses/Nr_sim,2),
GAINS=round(100*nr_gains/Nr_sim,2)
)
)

}
par(mfrow=c(2,2))
par(mar=c(2.5, 4.5,1.5, 1.5))

more_trading(10000,1,bad_scenarios = F)
more_trading(10000,4)
more_trading(10000,9)
more_trading(10000,19)
#more_trading(10000,99)


########################################################
########### Investigating Bad Scenarios ################
########################################################

par(mfrow=c(1,2))
########## Picking the worst Scenarios
#Nr of bad scenarios
bad_seeds<-NULL
for(s in 1:10000000){
set.seed(s)
Y<-AR_1(time_diff,N=n+1,mu[1],k=kk[1],sigma[1],plot=F,plot_threshold = 2*sigma,S_0=S_0)
if(abs(Y[n+1])>12){bad_seeds<-c(bad_seeds,s)}
}
length(bad_seeds)

bad_seeds_copy<-bad_seeds
bad_seeds<-bad_seeds_copy[1:1000]
### Investigating Success on these scenarios

### Pairs Trading ####
nr_gains<-0
nr_losses<-0
gain<-0
best<-0
worst<-0
ff<-numeric(length(bad_seeds))
for(s in bad_seeds){
  set.seed(s)
  S<-AR_1(time_diff,N=n+1,mu[1],k=kk[1],sigma[1],plot=F,plot_threshold = 2*sigma,S_0=S_0)
  D1<-Delta_strat(N=n+1,eps=2*sigma[ii],delta=0.001,S,mu[ii],intermediate_NR=3)
  f<-Gain(D1,(n+1),S)[n+1]
  ff[s]<-f
  gain<-gain+f
  if(f<0){nr_losses<-nr_losses+1}
  if(f>=0){nr_gains<-nr_gains+1}
  if(f>best){best<-f}
  if(f<worst){worst<-f}
}

round(gain/length(bad_seeds),2)
round(best,2)
round(worst,2)
round(100*nr_losses/length(bad_seeds),2)
round(100*nr_gains/length(bad_seeds),2)
hist(ff,xlab="Average Gain",main=paste("Pairs Trade, Bad Scenarios"),breaks=50)
############### Stat Arbitrage ###########
more_trading(1000,3,bad_scenarios=T)


