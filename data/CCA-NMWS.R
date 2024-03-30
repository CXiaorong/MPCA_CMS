#CCA-NMWS
#1.First select the function(1)-function(14) function and run it.
#2.load data
#3.Initialization parameters.
#4.Iterative loop


#function(1): Calculate the SI in the network
cor_sum4 <- function(Like,x){ 
  k <- length(x)
  result <- 0
  if(k!=1){
    for(i in 1:(k-1)){
      for(j in (i+1):k){
        temp <- as.numeric(Like[x[i],x[j]])
        result <- result+abs(temp)
      }
    }
  }
  s1=0
  for(i in 1:k){
    s1=s1+sum(Like[x[i],])
  }
  ss=result/(k*(k-1)/2)+(s1-2*result)/(k*maxl)
  return(ss)
}


#function(2): fitness function--CMS model
fitness <- function(A,E,nets_weight_data,Like,P,x){
  temp <- matrix(0,1,n)
  temp[x] <- 1
  y <- temp
  index <- which(y==1) 
  m1=nrow(A)
  a <- as.matrix(A[,index])
  a_indexsum <- rowSums(a)
  a_indexsum2=a_indexsum[which(a_indexsum>1)] 
  a_indexsum3=a_indexsum[which(a_indexsum>0)]
  
  a_colsum2= colSums(a)
  a_colsum <- colSums(a)/m1
  l1=length(a_indexsum3)
  l2=length(a_indexsum2)
  
  f_s=matrix(0,1,k+4)
  f_s[1,1:k]=index
    if(l1==0){ 
      f_s[1,k+2]=0 
  }else{
    f_s[1,k+2]=sum(1/(a_indexsum3))/(m1)
  }
  second_like=cor_sum4(second_like_data,index)
  f_s[1,k+1]=length(a_indexsum3)/m1
  f_s[1,k+3]=second_like
  beta=log10(k)
  f_s[1,k+4]=f_s[1,k+1]+f_s[1,k+2]+f_s[1,k+3]
  return(f_s)
}

fitness4 <- function(x){
  temp <- matrix(0,1,n)
  temp[x] <- 1
  y <- temp
  index <- which(y==1) 
  m1=nrow(SNVdata)
  a <- as.matrix(SNVdata[,index])
  a_indexsum <- rowSums(a)
  a_indexsum2=a_indexsum[which(a_indexsum>1)] 
  a_indexsum3=a_indexsum[which(a_indexsum>0)] 
  
  a_colsum2= colSums(a)
  a_colsum <- colSums(a)/m1 
  l1=length(a_indexsum3)
  l2=length(a_indexsum2)
  
  f_s=matrix(0,1,k+4)
  f_s[1,1:k]=index
    if(l1==0){   
      f_s[1,k+2]=0
  }else{
    f_s[1,k+2]=sum(1/(a_indexsum3))/(m1)
  }
  second_like=cor_sum4(second_like_data,index)
  f_s[1,k+1]=length(a_indexsum3)/m1
  f_s[1,k+3]=second_like
  beta=log10(k)
  f_s[1,k+4]=f_s[1,k+1]+f_s[1,k+2]+f_s[1,k+3]
  return(f_s)
}


#function(3): Compute the selection probability function.
select_order_fitness2 <- function(fit_vector){
  n <- length(fit_vector) 
  I <- order(fit_vector,decreasing = 'T')
  p <- matrix(0,1,n)
  s=sum(fit_vector)
  fun<-function(i){
    return(fit_vector[i]/s)
  }
  p <- sapply(1:n, fun)
  p_cumsum <- cumsum(p)#p=1,2,3,4;   p_cumsum=1,3,6,10
  random_data <- runif(1)
  temp <- which(p_cumsum>=random_data)
  index1 <- temp[1] 
  return(index1)
}

#function(4): crossover function1 
crossover2 <- function(parent1,parent2,n){ 
  temp22=parent1
  temp <- matrix(0,1,n) 
  temp[parent1] <- 1
  parent1 <- temp
  temp <- matrix(0,1,n)
  temp[parent2] <- 1
  parent2 <- temp
  newpop <- matrix(0,1,n)
  index <- which((parent1+parent2)==2)  
  newpop[index] <- 1
  parent1[index] <- 0
  parent2[index] <- 0
  temp <- which(parent1+parent2==1) 
  index <- sample(1:sum(parent1+parent2==1),sum(parent1+parent2))
  newpop[temp[index[1:(k-sum(newpop))]]]=1
  newpop <- which(newpop==1)
  return(newpop)
}

#function(5): Mutation function
mutation_SA <- function(SNVdata,GEdata,nets_weight_data,second_like_data,MinP,popTemp,n,N){
  popTemp2=matrix(0,2,k+4)
  popTemp2[1,]=popTemp
  pop_i=popTemp[1:k]
  
  houxuan<-c(1:n) 
  houxuan <- houxuan[!houxuan%in%pop_i] 
  houxuan=sample(houxuan,n-k)
  indx=sample(k,N)
  pop_j=pop_i
  for (i in 1:N) {
    pop_j[indx[i]]=houxuan[i]
  }
  popTemp2[2,1:(k+4)]=fitness(SNVdata,GEdata,nets_weight_data,second_like_data,MinP,pop_j)
  if (popTemp2[2,(k+4)]>=popTemp2[1,(k+4)])
    pop_i <- pop_j
  return(pop_i)
}

#function(7): select function
select2 <- function(pop3,n){
  I <- order(pop3[,k+4],decreasing = T) 
  pop3 <- pop3[I,]
  p <- vector(length = popsize)
  pop_next <- matrix(0, popsize,k+4)
  pop_next[1,] <- pop3[1,] 
  p=as.numeric(pop3[,k+4])/sum(as.numeric(pop3[,k+4]))
  ff <- (popsize%/%10)*8
  for(i in 2:ff){ 
    random_data=runif(1)
    p_cumsum <- cumsum(p)
    temp <- which(p_cumsum>=random_data)
    index1 <- temp[1]
    pop_next[i,] <- pop3[index1,]
  }
  ff2=(popsize-ff)*5
  mut_pop=matrix(0,ff2,k)
  for(i in 1:ff2){
    mut_pop[i,]=sample(1:n,k)
  }
  res <- parApply(cl,mut_pop,1,fitness4)
  pop_temp=t(res)
  I=order(pop_temp[,k+4],decreasing = 'T')
  pop_temp=pop_temp[I,]
  pop_next[(ff+1):popsize,]=pop_temp[1:(ff2/5),]
  return(pop_next)
}

select3 <- function(pop3){
  I <- order(pop3[,k+4],decreasing = T) 
  pop3 <- pop3[I,]
  p <- vector(length = popsize)
  pop_next <- matrix(0, popsize,k+4)
  pop_next[1,] <- pop3[1,] 
  p=as.numeric(pop3[,k+4])/sum(as.numeric(pop3[,k+4]))
  for(i in 2:popsize){ 
    random_data=runif(1)
    p_cumsum <- cumsum(p)
    temp <- which(p_cumsum>=random_data)
    index1 <- temp[1]
    pop_next[i,] <- pop3[index1,]
  }
  return(pop_next)
}

##function(8): GA2 function
GA2 <- function(pop){  
  I <- order(pop[,k+4],decreasing = T) 
  pop <- pop[I,]
  p <- vector(length = popsize)
  pop_next <- matrix(0, popsize,k+4)
  pop_next[1,] <- pop[1,] 
  fit_vector <-as.numeric(pop[,k+4]) 
  n=ncol(SNVdata)
  f1 <- popsize%/%2
  for(i in 2:f1){ 
    index <- select_order_fitness2(fit_vector) 
    index1 <- index[1]
    index2 <- index[2]
    new_ind <- crossover2(pop[index1,(1:k)],pop[index2,1:k],n)
     pop_next[i,1:k]=new_ind 
  }
  for(i in f1:popsize){ 
    fit_vector <-as.numeric(pop_next[,k+4])
    n <- length(fit_vector) 
    I <- order(fit_vector) 
    p <- matrix(0,1,n)
    s=sum(fit_vector)
    for(j in 1:n){
      p[I[j]] <- fit_vector[j]/s 
    }
    p_cumsum <- cumsum(p)
    random_data <- runif(1)
    temp <- which(p_cumsum>=random_data)
    index1 <- temp[1] 
    if(index1==1){
      index1=index1+1
    }
    n4 <- ncol(SNVdata)
    new_ind=mutation_SA(SNVdata,GEdata,nets_weight_data,second_like_data,MinP,pop[index1,],n4,1)
    pop_next[i,1:k]=new_ind
    # }
  }
  res <- parApply(cl,pop_next[,1:k],1,fitness4)
  pop_temp=t(res)
  return(pop_temp)
}

#function(9): GA6 function
GA6 <- function(pop){  
  r1=0.5 #交叉概率
  r2=0.6 #突变概率
  r3=0.6 
  I <- order(pop[,k+4],decreasing = T) 
  pop <- pop[I,]

  pop_next <- matrix(0, popsize,k+4)
  pop_next1 <- matrix(0, popsize,k+4)
  pop_next2 <- matrix(0, popsize,k+4)
  pop_next3 <- matrix(0, popsize,k+4)

  pop_next1[1,] <- pop[1,] 
  fit_vector <-as.numeric(pop[,k+4])
  n=ncol(SNVdata)

  for(i in 1:popsize){
    index <- select_order_fitness2(fit_vector)
    index1 <- index[1]
    index2 <- index[2]
    pop_next1[i,1:k] <- crossover2(pop[index1,(1:k)],pop[index2,1:k],n)
  }

  n4 <- ncol(SNVdata)
  for(i in 1:popsize){
    r5=runif(1)
    if(r5>r2){
      Kgene=pop[i,1:k] 
      beixuan<-c(1:n) 
      beixuan <- beixuan[!beixuan%in%Kgene] 
      new_ind=pop[i,1:k]
      if(k<3){
        index=sample(1:k,1)
        new_ind[index]=sample(beixuan,1)
      }else{
        index=sort(sample(1:(k-1),2)+1) 
        c1=index[1]
        c2=index[2]
        new_ind[c1:c2]=sample(beixuan,(c2-c1+1))
      }
      pop_next2[i,1:k]=new_ind
      # }
    }else{
      pop_next2[i,1:k]=pop[i,1:k]
    }
  }

  for(i in 1:popsize){
    r6=runif(1)
    if(r6>r3){
      index=sample(1:k,1) 
      c1=pop[,index]
      c2=base::unique(c1) 
      if(length(c2)<=k){
        c2=table(pop[,1:k])
        c2=as.numeric(names(c2))
      }
      new_ind=sample(c2,k)
      pop_next3[i,1:k]=new_ind 
    }else{
      pop_next3[i,1:k]=pop[i,1:k]
    }
  }

  #Three matrix connections
  pop_next=rbind(pop_next1,pop_next2,pop_next3)
  res <- parApply(cl,pop_next[,1:k],1,fitness4)
  pop_next=t(res)
  I <- order(pop_next[,k+4],decreasing = T)
  pop_next <- pop_next[I,]
  pop_next4=base::unique(pop_next) 
  n3=nrow(pop_next4)
  if(n3<popsize){
    n4=popsize-n3+1
    pop_next[n4:popsize,]=pop_next4
  }else{
    pop_next[1:popsize,]=pop_next4[1:popsize,]
  }
  return(pop_next[1:popsize,]) 
}



#function(10): Gene evolution of cooperative pool
GA9 <- function(pop){  
  I <- order(pop[,k+4],decreasing = T) 
  pop <- pop[I,]
  p <- vector(length = popsize)
  pop_next <- matrix(-1000, popsize,k+4)
  fit_vector <-as.numeric(pop[,k+4]) 
  n=ncol(SNVdata)
  f1 <- nrow(pop)
  pop_next[1,] <- pop[1,] #Retain the optimal individual
  
  cro_pop=matrix(0,popsize,k)
  for(i in 2:popsize){ 
    index <- select_order_fitness2(fit_vector) 
    index1 <- index[1]
    index2 <- index[2]
    cro_pop[i,] <- crossover2(pop[index1,(1:k)],pop[index2,1:k],n)
  }
  res <- parApply(cl,cro_pop[2:popsize,],1,fitness4)
  pop_next[2:popsize,]=t(res)
  
  fit_vector <-as.numeric(pop_next[,k+4])
  n2 <- length(fit_vector) 
  I <- order(fit_vector) #In order from small to large, mutation is carried out on those with poor adaptation value
  p <- matrix(0,1,n2)
  s=sum(fit_vector)
  for(j in 1:n2){
    p[I[j]] <- fit_vector[j]/s #Probability of being selected
  }
  p_cumsum <- cumsum(p)#p=1,2,3,4;   p_cumsum=1,3,6,10
  u=1
  for (j in 1:u) {
    random_data <- runif(1)#There are two functions that generate random numbers, and they are runif(),rnorm()
    temp <- which(p_cumsum>=random_data)
    index1 <- temp[1] #Take the first number that is closest to greater than the random number
    Kgene=pop_next[index1,1:k]
    beixuan<-c(1:n) 
    beixuan <- beixuan[!beixuan%in%Kgene] 
    beixuan=sample(beixuan,n-k)
    
    popTemp2=matrix(-1000,2,k+4)
    popTemp2[1,]=pop_next[index1,]
    temp=sample(1:k,1)
    # ff=floor((n-k)/2)
    ff=floor(n-k)
    # ff=n-k
    Kgene2=matrix(0,ff,k)
    for(i in 1:ff){
      Kgene[temp] <- beixuan[i]
      Kgene2[i,]=Kgene
    }
    
    res <- parApply(cl,Kgene2,1,fitness4)
    pop2=t(res)
    I <- order(pop2[,k+4],decreasing = T)
    pop2 <- pop2[I,]
    if(pop2[1,k+4]>pop_next[index1,k+4]){
      pop_next[index1,] <- pop2[1,]
    }
  }
  return(pop_next) #Select the pre-popsize individual
}

#function(11): Initial population
initPop2<-function(SNVdata,GEdata,nets_weight_data,second_like_data,MinP,k,n,popsize){
  x <- 1:n
  temp=matrix(0,popsize,k)
  for(i in 1:(popsize)){
    temp[i,] <- sample(x,k) #The sequence of n genes is shuffled so that only the first 1: k genes are taken at a time
    arr<<- append(arr,temp[i,])
  }
  pop=parApply(cl,temp,1,fitness4)
  pop=t(pop)
  return(pop) 
}

#function(12): significance test function
significance <- function(A,E,nets_weight_data,Like,P,subset_M){
  m <- nrow(A)
  w <- matrix(0,1,1000)
  n <- length(subset_M)
  for(j in 1:1000){
    A_temp <- A
    A_temp[,subset_M] <- 0
    for(i in 1:n){
      temp <- sum(A[,subset_M[i]])
      index <- round(runif(temp,min = 1,max = m))
      A_temp[index,subset_M[i]] <- 1
    }
    W1=fitness(A_temp,E,nets_weight_data,Like,P,subset_M)
    W[j]=W1[1,k+4]
    
  }
  w2=fitness(A,E,nets_weight_data,Like,P,subset_M)
  p <- sum(w>=w2[1,k+4])/1000
  return(p)
}


#start
#1. read data
#GBM

SNV_data<-read.csv('E:/data/GBM/90x440/SNVdata_440.csv')
GE_data<-read.csv('E:/data/GBM/90x440/GEdata_440.csv')
second_like_data<-read.csv('E:/sourceFiles-MPCA-CMS/data/GBM/90x440/second_like_data_GBM_440.csv')
rownames(second_like_data)=second_like_data[,1]
second_like_data=second_like_data[,-1]
second_like_data=as.matrix(second_like_data)

#OVCA
SNV_data<-read.csv('E:/data/OV/313x2547/OV_SNV_2547.csv')
GE_data<-read.csv('E:/data/OV/313x2547/OV_GE_2547.csv')
second_like_data<-read.csv('E:/data/OV/313x2547/second_like_data_hint_2547.csv')
rownames(second_like_data)=second_like_data[,1]
second_like_data=second_like_data[,-1]
second_like_data=as.matrix(second_like_data)

#KIRC
SNV_data<-read.csv('E:/data/KIRC/332x5804/SNV_332x5804.csv')
GE_data<-read.csv('E:/data/KIRC/332x5804/GE_332x5804.csv')
second_like_data<-read.csv('E:/data/KIRC/332x5804/second_like_data_332x5804.csv')
rownames(second_like_data)=second_like_data[,1]
second_like_data=second_like_data[,-1]
second_like_data=as.matrix(second_like_data)

#2. processing data
rownames(SNV_data)<-SNV_data[,1]
SNV_data<-SNV_data[,-1]
rownames(GE_data)<-GE_data[,1]
GE_data<-GE_data[,-1]
snv_colsum<-colSums(SNV_data) 
SNVdata<-SNV_data[,snv_colsum>1]
GEdata<-GE_data[,which(colnames(GE_data)%in%colnames(SNVdata))]

#2.Initialization Parameters
n <- ncol(SNVdata)
m <- nrow(SNVdata)
geneName <- colnames(SNVdata)
iteration <- 1000    #iterations
# install.packages('doParallel') #If the multi-core package does not exist, install it first
# library(doParallel) #You need to import it after installation
cl <- makeCluster(getOption("cl.cores", 4))      
registerDoParallel(cl)       #Process registration
clusterExport(cl, "n") 
clusterExport(cl, "SNVdata") 
clusterExport(cl, "GEdata")
clusterExport(cl, "cor_sum4") 
clusterExport(cl, "second_like_data") 

#Solve for the maximum degree
b=second_like_data
c=b
l=nrow(b)
maxl=0
ll=0
for (i in 1:l) {
  cc=sum(b[i,])
  if(maxl<cc){
    maxl=cc
    ll=i
  }
}
clusterExport(cl, "maxl")

v=1
result_parameter2=matrix(0,10,30)

for(k in 2:4){
  k <-3 #set K
  clusterExport(cl, "k") 
  print(k)
  popsize <- floor(log2(n^k))
  #3.Initial population
  t1 <- Sys.time()   #start counting
  arr=c()
  result_pop=initPop2(SNVdata,GEdata,nets_weight_data,second_like_data,MinP,k,n,2*popsize)
  P1=result_pop[1:popsize,]
  P2=result_pop[(popsize+1):(2*popsize),]
  x <- 1:popsize
  
  j=1
  R=1
  t=0
  I=order(P1[,k+4],decreasing = 'T')
  Pbest=P1[I,]
  source_pop=matrix(-1000,popsize,k+4)
  while(j<=iteration & R<=10) { 
    #Determine whether convergence
    best_ind=Pbest[1,1:k]
    temp_pop <- as.matrix(P1[!duplicated(P1[,1:k]),,drop=T])
    temp_pop2 <- as.matrix(P2[!duplicated(P2[,1:k]),,drop=T])
    R1=nrow(temp_pop)
    R2=nrow(temp_pop2)
    
    np1=ncol(temp_pop)
    np2=ncol(temp_pop2)
    
    if(R1*np1==(k+4)||R2*np2==(k+4)){
      #print('finish').
      break
    }
    condition=(j%%5==0) || (R1<=popsize/2 || R2<=popsize/2)
    if(condition){
      I=order(P1[,k+4],decreasing = 'T')
      P1=P1[I,]
      I=order(P2[,k+4],decreasing = 'T')
      P2=P2[I,]
      source_pop[1:(popsize/2),]=P1[1:(popsize/2),]
      source_pop[(popsize/2+1):popsize,]=P2[1:(popsize/2),]
      source_pop=GA9(source_pop)
      
      P1=source_pop[1:popsize,]
      I=order(source_pop[,k+4],decreasing = 'T')
      source_pop=source_pop[I,]
      Pbest=source_pop[1:popsize,]
      ll=length(which(best_ind==Pbest[1,1:k])==TRUE) #Determine whether it is equal
      if(ll==k){
        R=R+1
      }else{
        R=0
      }
    }
    P1=GA2(P1)
    P2=GA6(P2)

    P1 <- select2(P1,n)
    P2 <- select2(P2,n)
    #print(j)
    j=j+1
  }
  t2 <- Sys.time() 
  diff <- difftime(t2,t1,units="secs")
  
  pop=matrix(0,2*popsize,k+4)
  pop[1:popsize,]=P1
  pop[(popsize+1):(2*popsize),]=P2
  I <- order(pop[,(k+4)],decreasing = T) 
  pop <- pop[I,]
  maxpop <- pop
  maxpop[,1:k] <- apply(pop[,1:k,drop=T],2,function(x) {geneName[x]})

  #output the optimal solution.
  print(maxpop[1,1:k])
  print(maxpop[1,(k+4)])
 
  v=v+1
}
stopCluster(cl)
