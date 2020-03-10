# Setting some fix values
n<-10 # Nr.of time steps
N<-100 # Nr. of different P's

S_0<-100

# Defining the Model

### Choice
u<-rep(1.1,N)
d<-rep(1/u,N)
p<-seq(0.4,0.6,length=N)

## Values
x<-matrix(0,N,n/2+1)
y<-matrix(0,N,n+1)


# probabilities
prob_set<-array(0,dim=c(N,(n/2)+1,n+1))
for(i in 1:N){
  # Assigning Values
  for(l in 0:(n/2)){
    x[i,l+1]<-S_0*u[i]^l*d[i]^(n/2-l)
  }
  for(l in 0:(n)){
    y[i,l+1]<-S_0*(u[i]^l)*d[i]^(n-l)
  }
  for(j in 1:((n/2)+1)){
    for(k in 0:(n/2)){
      # Assigning Probabilities
      prob_set[i,j,k+j]<-dbinom(j-1,n/2,p[i])*dbinom(k,n/2,p[i])
    }
  }
}

check_statarb(x[10,],y[10,],prob_set[10,,],S_0)

# Function Definition
func1<-function(x,y){max(0.5*(x+y)-S_0,0)}

# Computation of the Bounds
limit<-Inf
# Robust Bounds
lower_bound_nosa<-super_rep(x,y,prob_set,func=func1,stat_arb=F,lower=T,S_0,limit)
lower_bound<-super_rep_robust(x,y,prob_set,func=func1,stat_arb=T,lower=T,S_0,limit)
upper_bound<-super_rep_robust(x,y,prob_set,func=func1,stat_arb=T,lower=F,S_0,limit)
upper_bound_nosa<-super_rep(x,y,prob_set,func=func1,stat_arb=F,lower=F,S_0,limit)

print(upper_bound)
print(lower_bound)


Delta_0_lower<-lower_bound$Delta_0
Delta_1_lower<-splinefun(lower_bound$x,lower_bound$Delta_1,method="natural")
plot(lower_bound$x,Delta_1_lower(lower_bound$x),type="l")

#Simulation

#
# 


## Statistical Arbitrage
nr_gains<-0
nr_losses<-0
d_lower<-lower_bound$d
gain<-0
Nr_sim<-1000000
best<-0
worst<-0
for(i in 1:Nr_sim){
  random_p<-rbinom(n,1,runif(1,0.4,0.6))
  S<-S_0*cumprod(u[1]*(random_p==1)+d[1]*(random_p==0))
  f<-(func1(S[n/2],S[n])-Delta_0_lower*(S[n/2]-S_0)-Delta_1_lower(S[n/2])*(S[n]-S[n/2])-9)
  gain<-gain+f
  if(f<0){nr_losses<-nr_losses+1}
  if(f>0){nr_gains<-nr_gains+1}
  if(f>best){best<-f}
  if(f<worst){worst<-f}
}
gain/Nr_sim
best
worst
nr_losses/1000000
nr_gains/1000000

#### Too high price
Delta_0_upper<-upper_bound$Delta_0
Delta_1_upper<-splinefun(upper_bound$x,upper_bound$Delta_1,method="natural")
plot(upper_bound$x,Delta_1_lower(upper_bound$x),type="l")

#Simulation

#
# 


## Statistical Arbitrage
nr_gains<-0
nr_losses<-0
d_upper<-upper_bound$d
gain<-0
Nr_sim<-1000000
best<-0
worst<-0
for(i in 1:Nr_sim){
  random_p<-rbinom(n,1,runif(1,0.4,0.6))
  S<-S_0*cumprod(u[1]*(random_p==1)+d[1]*(random_p==0))
  f<-(-func1(S[n/2],S[n])+Delta_0_upper*(S[n/2]-S_0)+Delta_1_upper(S[n/2])*(S[n]-S[n/2])+20)
  gain<-gain+f
  if(f<0){nr_losses<-nr_losses+1}
  if(f>0){nr_gains<-nr_gains+1}
  if(f>best){best<-f}
  if(f<worst){worst<-f}
}
gain/Nr_sim
best
worst
nr_losses/1000000
nr_gains/1000000


# Pure Arbitrage
Delta_0_upper_nosa<-upper_bound_nosa$Delta_0
Delta_1_upper_nosa<-splinefun(upper_bound_nosa$x,upper_bound_nosa$Delta_1,method="natural")
plot(upper_bound_nosa$x,Delta_1_upper_nosa(upper_bound_nosa$x),type="l")

d_upper_nosa<-upper_bound_nosa$d
gain<-0
Nr_sim<-100000
best<-0
worst<-0
for(i in 1:Nr_sim){
  random_p<-rbinom(n,1,runif(1,0.4,0.6))
  S<-S_0*cumprod(u[1]*(random_p==1)+d[1]*(random_p==0))
  f<-(60-func1(S[n/2],S[n])+Delta_0_upper_nosa*(S[n/2]-S_0)+Delta_1_upper_nosa(S[n/2])*(S[n]-S[n/2]))
  gain<-gain+f
  if(f>best){best<-f}
  if(f<worst){worst<-f}
}
gain/Nr_sim
best
worst
# 
# # ##### Plotting
# random_p<-rbinom(n,1,p[1])
# S<-S_0*cumprod(u[1]*(random_p==1)+d[1]*(random_p==0))
# plot(1:n,S,type="l",ylim=c(50,150))
# points(rep(n,n+1),y[1,])
# points(rep(n/2,n/2+1),x[1,])
# for(j in 1:10){
#   random_p<-rbinom(n,1,p[1])
#   S<-S_0*cumprod(u[1]*(random_p==1)+d[1]*(random_p==0))
#   lines(1:n,S)
# }
# random_p<-rbinom(n,1,p[N])
# S<-S_0*cumprod(u[N]*(random_p==1)+d[N]*(random_p==0))
# plot(1:n,S,type="l")
# 
# 

