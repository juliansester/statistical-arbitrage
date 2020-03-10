# Loading necessary packages
library('linprog')
library('gurobi')


############################################################################
########################## Defining the Function for fixed P ###############
############################################################################
stat_arb_bounds_fixp<-function(x,y,prob,func=function(x,y){abs(x-y)},stat_arb=T,min=F,S_0=1,epsilon=0){

n1<-length(x) # Length of the vector with the values of the first marginal
n2<-length(y) # Length of the vector with the values of the 2nd marginal

r<-numeric(n1+3)# R.H.S vector of the Matrix
A<-matrix(0,ncol=(n1*n2),nrow=(n1+3)) #L.H.S. Vector / Setting the size

# Martingale Conditions
for(i in 1:n1){
  a<-array(0,dim=c(n1,n2)) # Dummy
  for(j in 1:n2){
    a[i,j]<-(y[j]-x[i])
  }
  A[i,]<-as.vector(a)
}

# Martingale constraint at 0
a<-array(0,dim=c(n1,n2)) # Dummy
for(i in 1:n1){
  for(j in 1:n2){
    a[i,j]<-x[i]
  }
}
A[n1+1,]<-as.vector(a)
r[n1+1]<-S_0

# Measure constraint
A[n1+2,]<-rep(1,n1*n2)
r[n1+2]<-1


# not equivalent
A[n1+3,]<-as.vector(rep(0,n1*n2)+(as.vector(prob)==0))
r[n1+3]<-0


nr_nonzero<-0

if(stat_arb){
#measurability constraint
for(i in 1:n1){
  for(j in 1:n2){
    for(k in 1:n1){
      if(k != i){
      a<-array(0,dim=c(n1,n2))
      a[i,j]<-prob[k,j]
      a[k,j]<-(-prob[i,j])
      A<-rbind(A,as.vector(a))
      r<-c(r,0)
      }
    }
  }
}
  #ensuring non-zero entries

  for(i in 1:(n1*n2)){
    if(as.vector(prob)[i]>0){
      a<-rep(0,n1*n2)
      a[i]<-1
      A<-rbind(A,a)
      nr_nonzero<-nr_nonzero+1
      r<-c(r,epsilon)
    }
  }
}





#Cost Function
costs<-array(0,dim=c(n1,n2))
for(i in 1:n1){
  for(j in 1:n2){
    costs[i,j]<-func(x[i],y[j]) 
  }
}

costs<-as.vector(costs)
if(1-min){
  opt<-gurobi( list(A=A,obj=costs,lb=rep(0,n1*n2),ub=rep(1,n1*n2),modelsense="max",rhs=r,sense=c(rep("=",dim(A)[[1]]-nr_nonzero),rep(">",nr_nonzero))), params=  list( OutputFlag=0))
}
if(min){
  opt<-gurobi( list(A=A,obj=costs,lb=rep(0,n1*n2),ub=rep(1,n1*n2),modelsense="min",rhs=r,sense=c(rep("=",dim(A)[[1]]-nr_nonzero),rep(">",nr_nonzero))), params=  list( OutputFlag=0))
}

q<-array(opt$x,dim=c(n1,n2))
price<-opt$objval


return(list(Q=q,Price=price))
                      
}


##################################################################################
############## Generalization to a range of probability measures #################
##################################################################################

stat_arb_bounds<-function(x,y,prob_set,func=function(a,b){abs(a-b)},stat_arb=T,min=F,S_0=1){
  # Initial value
  s<-stat_arb_bounds_fixp(x[1,],y[1,],prob_set[1,,],func,stat_arb,min,S_0)
  values<-s[[2]]
  measures<-s[[1]]
  p<-prob_set[1,,]
  xx<-x[1,]
  yy<-y[1,]
  if(min){
  for(i in 2:dim(prob_set)[1]){
    s<-stat_arb_bounds_fixp(x[i,],y[i,],prob_set[i,,],func,stat_arb,min,S_0)
    if(s[[2]]<values){
      values<-s[[2]]
      measures<-s[[1]]
      p<-prob_set[i,,]
      xx<-x[i,]
      yy<-y[i,]
    }
  }
  }
  if(min==F){
    for(i in 2:dim(prob_set)[1]){
      s<-stat_arb_bounds_fixp(x[i,],y[i,],prob_set[i,,],func,stat_arb,min,S_0)
      if(s[[2]]>values){
        values<-s[[2]]
        measures<-s[[1]]
        p<-prob_set[i,,]
        xx<-x[i,]
        yy<-y[i,]
      }
    }
  }
  return(list(P=p,Q=measures,Optimal_Value=values,x=xx,y=yy))
}

##################################################################################
### Computation of Super-Replication Strategies
##################################################################################

super_rep_fixp<-function(x,y,prob,func=function(x,y){abs(x-y)},stat_arb=T,lower=T,S_0=1,limit=Inf){
  library('linprog')
  n1<-length(x) # Length of the vector with the values of the first marginal
  n2<-length(y) # Length of the vector with the values of the 2nd marginal
  #costs
  costs<-array(0,dim=c(n1,n2))
  for(i in 1:n1){
    for(j in 1:n2){
      costs[i,j]<-func(x[i],y[j]) 
    }
  }
  
  if(stat_arb){
  #Right-Hand-Side Vector
  r<-rep(0,n2)
  #lhs
  A<-matrix(0,n2,n1+2)
  #Conditions
  for(j in 1:n2){
    if(sum(prob[,j])>0)
      {
    for(i in 1:n1){
      #Delta_1(x_i)
      A[j,i]<-((prob[i,j]/sum(prob[,j]))*(y[j]-x[i]))
    }
    #Delta_0
    A[j,n1+1]<-(sum(x*prob[,j]/sum(prob[,j]))-S_0)
    #d
    A[j,n1+2]<-1
    # right-hand side: E_P[c~|~S_2]
    r[j]<-(sum(costs[,j]*prob[,j]/sum(prob[,j])))
    }
  }
  # Additional Constraint E_P[f]>E_P[c]

  if(lower){
  opt<-gurobi( list(A=A,obj=c(rep(0,n1+1),1),lb=rep(-limit,n1+2),ub=rep(limit,n1+2),modelsense="max",rhs=r,sense='<'),
               params=  list( DualReductions=0, OutputFlag=0))
  }
  else{
    opt<-gurobi( list(A=A,obj=c(rep(0,n1+1),1),lb=rep(-limit,n1+2),ub=rep(limit,n1+2),modelsense="min",rhs=r,sense='>'),
                 params=  list(DualReductions=0, OutputFlag=0))
  }
  list(d=opt$objval,Delta_0=opt$x[n1+1],Delta_1=opt$x[1:n1])
  }
  else{
    #Right-Hand-Side Vector
    r<-rep(0,n1*n2)
    #lhs
    A<-matrix(0,n1*n2,n1+2)
    #Conditions
    for(i in 1:n1){
      for(j in 1:n2){
        #Delta_1(x_i)
        A[i+n1*(j-1),i]<-(y[j]-x[i])

      #Delta_0
      A[i+n1*(j-1),n1+1]<-(x[i]-S_0)
      #d
      A[i+n1*(j-1),n1+2]<-1
      r[i+n1*(j-1)]<-costs[i,j]
      }
    }
    

    
    if(lower){
      opt<-gurobi( list(A=A,obj=c(rep(0,n1+1),1),lb=rep(-limit,n1+2),ub=rep(limit,n1+2),modelsense="max",rhs=r,sense='<'),
                   params=  list( OutputFlag=0))
    }
    else{
      opt<-gurobi( list(A=A,obj=c(rep(0,n1+1),1),lb=rep(-limit,n1+2),ub=rep(limit,n1+2),modelsense="min",rhs=r,sense='>'),
                   params=  list( OutputFlag=0))
    }
    list(d=opt$objval,Delta_0=opt$x[n1+1],Delta_1=opt$x[1:n1])
  }
}



##################################################################################
### Computation of Robust Super-Replication Strategies
##################################################################################

super_rep_robust<-function(x,y,prob_set,func=function(x,y){abs(x-y)},stat_arb=T,lower=F,S_0=1,limit=Inf){
  library('linprog')
  n1<-length(x[1,]) # Length of the vector with the values of the first marginal: All x have the same length
  n2<-length(y[1,]) # Length of the vector with the values of the 2nd marginal: All y have the same length
  N<-dim(prob_set)[1] # Number of considered probabilities
  #costs
  costs<-array(0,dim=c(N,n1,n2))
  for (k in 1:N){
    for(i in 1:n1){
      for(j in 1:n2){
        costs[k,i,j]<-func(x[k,i],y[k,j]) 
      }
    }
  }
    #Right-Hand-Side Vector
    r<-rep(0,n2*N)
    #lhs
    A<-matrix(0,n2*N,n1+2)
    #Conditions
    for(k in 1:N){
      for(j in 1:n2){
        if(sum(prob_set[k,,j])>0)
        {
          for(i in 1:n1){
            A[(k-1)*n2+j,i]<-((prob_set[k,i,j]/sum(prob_set[k,,j]))*(y[k,j]-x[k,i]))
          }
          #Delta_0
          A[(k-1)*n2+j,n1+1]<-(sum(x[k,]*prob_set[k,,j]/sum(prob_set[k,,j]))-S_0)
          #d
          A[(k-1)*n2+j,n1+2]<-1
          # right-hand side: E_P[c~|~S_2]
          r[(k-1)*n2+j]<-(sum(costs[k,,j]*prob_set[k,,j]/sum(prob_set[k,,j])))

        }
      }
    }
    # Additional Constraint E_P[f]>E_P[c]
    
    if(lower){
      opt<-gurobi( list(A=A,obj=c(rep(0,n1+1),1),lb=rep(-limit,n1+2),ub=rep(limit,n1+2),modelsense="max",rhs=r,sense='<'),
                   params=  list( DualReductions=0, OutputFlag=0))
    }
    else{
      opt<-gurobi( list(A=A,obj=c(rep(0,n1+1),1),lb=rep(-limit,n1+2),ub=rep(limit,n1+2),modelsense="min",rhs=r,sense='>'),
                   params=  list(DualReductions=0, OutputFlag=0))
    }
    return(list(d=opt$objval,Delta_0=opt$x[n1+1],Delta_1=opt$x[1:n1],x=x[1,],y=y[1,]))
}


##################################################################################
# Function for non-robust super-replication strategies
##################################################################################

super_rep<-function(x,y,prob_set,func=function(x,y){abs(x-y)},stat_arb=T,lower=F,S_0=1,limit=Inf){
  s<-super_rep_fixp(x[1,],y[1,],prob_set[1,,],func,stat_arb,lower,S_0,limit)
  d<-s$d
  Delta_0<-s$Delta_0
  Delta_1<-s$Delta_1
  xx<-x[1,]
  yy<-y[1,]
  p<-prob_set[1,,]
  if(lower){
    for(i in 2:dim(prob_set)[1]){
      s<-super_rep_fixp(x[i,],y[i,],prob_set[i,,],func,stat_arb,lower,S_0,limit)
      if(s$d<d){
        d<-s$d
        Delta_0<-s$Delta_0
        Delta_1<-s$Delta_1
        p<-prob_set[i,,]
        xx<-x[i,]
        yy<-y[i,]
      }
    }
  }
  if(lower==F){
    for(i in 2:dim(prob_set)[1]){
      s<-super_rep_fixp(x[i,],y[i,],prob_set[i,,],func,stat_arb,lower,S_0,limit)
      if(s$d>d){
        d<-s$d
        Delta_0<-s$Delta_0
        Delta_1<-s$Delta_1
        p<-prob_set[i,,]
        xx<-x[i,]
        yy<-y[i,]
      }
    }
  }
  return(list(d=d,P=p,Delta_0=Delta_0,Delta_1=Delta_1,x=xx,y=yy))
}


################################################
### Function checking for Statistical Arbitrage
################################################
check_statarb<-function(x,y,prob,S_0=1,epsilon=0){
    if(super_rep_fixp(x,y,prob,func=function(x,y){epsilon},stat_arb=T,lower=F,S_0,limit=10e+10)$d<0){print("There is Statistical Arbitrage")}
  else{print("No Statistical Arbitrage!")}
  }