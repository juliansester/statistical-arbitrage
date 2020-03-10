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





N<-10 # different Parameters
n<-100 # time steps 
n_discr<-50
time_diff<-1
kk<-rep(0.001,N)
# Setting Parameters
mu<-seq(-1,1,length=N)
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
AR_1(t=time_diff,N=n+1,mu=mu[1],k=kk[1],sigma=sigma[1],plot=T,plot_threshold=2*sigma[1],S_0=S_0)
# Function Definition

func1<-function(x,y){max(y-S_0,0)}

# Computation of the Bounds
limit<-Inf

# Robust Bounds
check_statarb(x[1,],y[1,],prob_set[1,,],S_0)
super_rep_fixp(x[1,],y[1,],prob_set[1,,],func=function(x,y){0},stat_arb=T,lower=F,S_0,limit)
lower_bound_nosa<-super_rep(x,y,prob_set,func=func1,stat_arb=F,lower=T,S_0,limit)
lower_bound<-super_rep_robust(x,y,prob_set,func=func1,stat_arb=T,lower=T,S_0,limit)
upper_bound<-super_rep_robust(x,y,prob_set,func=func1,stat_arb=T,lower=F,S_0,limit)
upper_bound_nosa<-super_rep(x,y,prob_set,func=func1,stat_arb=F,lower=F,S_0,limit)

# Simulations
## Statistical Arbitrage
par(mfrow=c(2,2))
par(mar=c(2.5, 4.5,1.5, 1.5))
Delta_0_upper<-upper_bound$Delta_0
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
  S<-AR_1(t=time_diff,N=n+1,mu=mu[ii],k=kk[ii],sigma=sigma[ii],plot=F,plot_threshold=2*sigma,S_0=S_0)
  f<-(6-func1(S[n/2],S[n])+Delta_0_upper*(S[n/2]-S_0)+Delta_1_upper(S[n/2])*(S[n]-S[n/2]))
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
hist(ff,xlab="Average Gain",main="Statistical Arbitrage",breaks=50)

########## Delta Strategies
Delta_strat<-function(N,eps=20,delta=0.001,Y,mu=0){
  Delta<-numeric(N)
  for(i in 2:N){
    if((Y[i]>=(mu+eps))*(Delta[i-1]==0)){
      Delta[i]<-(Delta [i-1]-1)
    }
    else if((Y[i]>=(mu-delta))*(Delta[i-1]==1)){
      Delta[i]<-(Delta [i-1]-1)
    }
    else if((Y[i]<=(mu-eps))*(Delta[i-1]==0)){
      Delta[i]<-(Delta [i-1]+1)
    }
    else if((Y[i]<=(mu+delta))*(Delta[i-1]==(-1))){
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
# Y<-AR_1(time_diff,N=n+1,mu[1],k=kk[1],sigma[1],plot=F,plot_threshold = 2*sigma,S_0=S_0)
# #D1<-Delta_strat(N,eps=2*sigma,delta=0.001,Y)
# D1<-Delta_strat(N=n+1,eps=2*sigma[1],delta=0.001,Y,mu[1])
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
#1)epsilon=sigma

nr_gains<-0
nr_losses<-0
gain<-0
Nr_sim<-1000000
best<-0
worst<-0
ff<-numeric(Nr_sim)
for(i in 1:Nr_sim){
  ii<-sample(1:N,1) # choosing random parameter
  S<-AR_1(t=time_diff,N=n+1,mu=mu[ii],k=kk[ii],sigma=sigma[ii],plot=F,plot_threshold=2*sigma[ii],S_0=S_0)
  D1<-Delta_strat(N=n+1,eps=1*sigma[ii],delta=0.001,S,mu[ii])
  f<-Gain(D1,(n+1),S)[n+1]
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
hist(ff,xlab="Average Gain",main="Pairs trading with Epsilon =sigma",breaks=50)

#2)epsilon=2*sigma

nr_gains<-0
nr_losses<-0
gain<-0
Nr_sim<-1000000
best<-0
worst<-0
ff<-numeric(Nr_sim)
for(i in 1:Nr_sim){
  ii<-sample(1:N,1) # choosing random parameter
  S<-AR_1(t=time_diff,N=n+1,mu=mu[ii],k=kk[ii],sigma=sigma[ii],plot=F,plot_threshold=2*sigma[ii],S_0=S_0)
  D1<-Delta_strat(N=n+1,eps=2*sigma[ii],delta=0.001,S,mu[ii])
  f<-Gain(D1,(n+1),S)[n+1]
  gain<-gain+f
  ff[i]<-f
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
hist(ff,xlab="Average Gain",main="Pairs trading with Epsilon =2*sigma",breaks=50)

#3) epsilon=3*sigma

nr_gains<-0
nr_losses<-0
gain<-0
Nr_sim<-1000000
best<-0
worst<-0
ff<-numeric(Nr_sim)
for(i in 1:Nr_sim){
  ii<-sample(1:N,1) # choosing random parameter
  S<-AR_1(t=time_diff,N=n+1,mu=mu[ii],k=kk[ii],sigma=sigma[ii],plot=F,plot_threshold=2*sigma[ii],S_0=S_0)
  D1<-Delta_strat(N=n+1,eps=3*sigma[ii],delta=0.001,S,mu[ii])
  f<-Gain(D1,(n+1),S)[n+1]
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
hist(ff,xlab="Average Gain",main="Pairs trading with Epsilon =3*sigma",breaks=50)



