################################## loading all the possible required package
library(pixmap)
library(RGCCA)
library(cluster)
library(kernlab)
library(CCA)
library(mclust)
library(class)
library("R.matlab")
library(R1magic)
library(psych)
library(fpc)
library(caret) ##### confusionMatrix
library(RSNNS)
library(MASS)  ###lda
library(e1071)
library(ggplot2)

######################################################################
#################OBJ gap as stopping criteria  #######################
######################################################################
  
######################################################################
######################################################################
######################   R1 norm realization  ########################
######################################################################
######################################################################
####################
BFPCA.obj <- function(Data, p, q, gamma0, epsilon, init ){  ##Data is the sample list, epsilon is the threldhold 
  ## p*q is the dims of feature, gamma0 not necessary
  timestart<-Sys.time()
  N <- length(Data)  ### sample size
  m <- dim(Data[[1]])[1]  ## number of rows
  n <- dim(Data[[1]])[2]  ## number of columns
  ######initialize U,V
  if ( init == "way1" )  {# way 1
    U <- matrix(0, p, m)
    V <- matrix(0, n, q)
    for(i in 1:p)
      U[i,i] <- 1
    
    for(i in 1:q)
      V[i,i] <- 1  }
  #########initialization with slicing learning i.e. HOSVD
  ## way 2
  else { X.col <- rep(0, n)
  for (j in 1:N)
    X.col <- rbind(X.col, Data[[j]])
  X.col <- X.col[-1,]
  
  G.col <- t(X.col) %*% X.col
  ev.col <- eigen(G.col)
  V <- Re(ev.col$vectors[,1:q])
  
  X.row <- rep(0, m)
  for (j in 1:N)
    X.row <- cbind(X.row, Data[[j]])
  X.row <- X.row[,-1]
  
  G.row <- X.row %*% t(X.row)
  ev.row <- eigen(G.row)
  U <- Re(t(ev.row$vectors[,1:p])) }
  ##################### other initialization
  U.past <- matrix(0, p, m)
  V.past <- matrix(0, n, q)
  
  t <- 1
  d <- rep(0,N)
  s <- 0
  
  s.past <- 0
  for (i in 1:N)
    s.past <- s.past + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")
  
  obj.track <- NULL
  UVGap.track <- NULL
  while (round(abs(s - s.past),20 )> epsilon ) { ## abs(s.past - s)
    ######output summary
    
    print(paste("number of iteration:",t))
    print(paste("past objective value:",s.past) )
    # print(paste("objective gap:",s.past - s  ))
    # cat("\n")
    # 
    s.past <- 0
    for (i in 1:N)
      s.past <- s.past + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")
    
    obj.track[t] <- s.past
    #####get d
    U.past <- U
    V.past <- V
    for(i in 1:N){
      d[i] <- 1/(norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + gamma0) }
    
    ##### solve optimization problem1
    H <- matrix(0, n, n)
    for(i in 1:N){
      H <- H + t(Data[[i]]) %*% t(U) %*% U %*% Data[[i]] * d[i] }
    ev.H <- eigen(H)
    
    # nor.matrix <- matrix(0,q,q)
    # for (i in 1:q){
    #   nor.matrix[i,i] <- sign(Re( ev.H$vectors[1, i]))
    # }
    #   
    # V <- as.matrix(Re( ev.H$vectors[,1:q]) %*% nor.matrix)  ###update U
    V <- as.matrix(Re( ev.H$vectors[,1:q]))
    ##### solve optimization problem2
    G <- matrix(0, m, m)
    for(i in 1:N){
      G <- G + Data[[i]] %*% V %*% t(V) %*% t(Data[[i]]) * d[i] }
    ev.G <- eigen(G)
    
    #  nor.matrix2 <- matrix(0,p,p)
    # for (i in 1:p){
    #   nor.matrix2[i,i] <- sign(Re( ev.G$vectors[1, i]))
    # }
    # U <- as.matrix(Re( t(ev.G$vectors[,1:p] %*%  nor.matrix2)))  ####update V
    U <- as.matrix(Re( t(ev.G$vectors[,1:p])))
    #####update object value
    s <- 0
    for (i in 1:N)
      s <- s + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")
    
    UVGap.track[t] <- norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
      norm(V %*% t(V) - V.past %*% t(V.past), type = "F")
    
    
    print(paste("now objective value:", s))
    print(paste("objective gap:",round(s.past - s, 20  )))
    print(paste("U,V gap",UVGap.track[t]))
    cat("\n")
    t <- t + 1
    if(t > 10^4) break
  }
  
  timeend<-Sys.time()
  runningtime<-timeend-timestart
    
  result <- list(U = U, V = V, t = t - 1, Obj = obj.track, UV.Gap = UVGap.track,Time = runningtime)
  return(result)
}
######################################################################
######################################################################
######################   R1 norm realization  ########################
######################     regularization     ########################
######################################################################
######################################################################
BFPCA.obj.Re <- function(Data, p, q, gamma0, epsilon, init, init.reg, lambda){  ##Data is the sample list, epsilon is the threldhold 
  timestart<-Sys.time()
  ## p*q is the dims of feature, gamma0 not necessary
  N <- length(Data)  ### sample size
  m <- dim(Data[[1]])[1]  ## number of rows
  n <- dim(Data[[1]])[2]  ## number of columns
  ######initialize U,V
  if ( init == "way1" )  {# way 1
    U <- matrix(0, p, m)
    V <- matrix(0, n, q)
    for(i in 1:p)
      U[i,i] <- 1
    
    for(i in 1:q)
      V[i,i] <- 1  }
  #########initialization with slicing learning i.e. HOSVD
  ## way 2
  else { X.col <- rep(0, n)
  for (j in 1:N)
    X.col <- rbind(X.col, Data[[j]])
  X.col <- X.col[-1,]
  
  G.col <- t(X.col) %*% X.col
  ev.col <- eigen(G.col)
  V <- Re(ev.col$vectors[,1:q])
  
  X.row <- rep(0, m)
  for (j in 1:N)
    X.row <- cbind(X.row, Data[[j]])
  X.row <- X.row[,-1]
  
  G.row <- X.row %*% t(X.row)
  ev.row <- eigen(G.row)
  U <- Re(t(ev.row$vectors[,1:p])) }
  ##################### other initialization
  U.past <- matrix(0, p, m)
  V.past <- matrix(0, n, q)
  
  t <- 1
  d <- rep(0,N)
  s <- 0
  
  s.past <- 0
  for (i in 1:N)
    s.past <- s.past + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")
  
  obj.track <- NULL
  UVGap.track <- NULL
  while (round(abs(s - s.past),20 )> epsilon ) { ## abs(s.past - s)
    ######output summary
    
    print(paste("number of iteration:",t))
    print(paste("past objective value:",s.past) )
    # print(paste("objective gap:",s.past - s  ))
    # cat("\n")
    # 
    if (init.reg == "adaptive"){
      lambda[2] <- 1/norm(t(U) %*% U - t(U.past) %*% U.past, type = "F")
      lambda[1] <- 1/norm(V %*% t(V) - V.past %*% t(V.past), type = "F")
    }
    
    s.past <- 0
    for (i in 1:N)
      s.past <- s.past + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")
    
    obj.track[t] <- s.past
    #####get d
    U.past <- U
    V.past <- V
    for(i in 1:N){
      d[i] <- 1/(norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + gamma0) }
    
    ##### solve optimization problem1
    H <- - lambda[1] * (diag(n) - V %*% t(V))
    for(i in 1:N){
      H <- H + t(Data[[i]]) %*% t(U) %*% U %*% Data[[i]] * d[i] }
    ev.H <- eigen(H)
    
    # nor.matrix <- matrix(0,q,q)
    # for (i in 1:q){
    #   nor.matrix[i,i] <- sign(Re( ev.H$vectors[1, i]))
    # }
    #   
    # V <- as.matrix(Re( ev.H$vectors[,1:q]) %*% nor.matrix)  ###update U
    V <- as.matrix(Re( ev.H$vectors[,1:q]))
    ##### solve optimization problem2
    G <- - lambda[2] * (diag(m) - t(U) %*% U)
    for(i in 1:N){
      G <- G + Data[[i]] %*% V %*% t(V) %*% t(Data[[i]]) * d[i] }
    ev.G <- eigen(G)
    
    #  nor.matrix2 <- matrix(0,p,p)
    # for (i in 1:p){
    #   nor.matrix2[i,i] <- sign(Re( ev.G$vectors[1, i]))
    # }
    # U <- as.matrix(Re( t(ev.G$vectors[,1:p] %*%  nor.matrix2)))  ####update V
    U <- as.matrix(Re( t(ev.G$vectors[,1:p])))
    #####update object value
    s <- 0
    for (i in 1:N)
      s <- s + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")
    
    UVGap.track[t] <- norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
      norm(V %*% t(V) - V.past %*% t(V.past), type = "F")
    
    
    print(paste("now objective value:", s))
    print(paste("objective gap:",round(s.past - s, 20  )))
    print(paste("U,V gap",UVGap.track[t]))
    cat("\n")
    t <- t + 1
    if(t > 10^4) break
  }
  timeend<-Sys.time()
  runningtime<-timeend-timestart
  
  result <- list(U = U, V = V, t = t - 1, Obj = obj.track, UV.Gap = UVGap.track, Time = runningtime)
  return(result)
}


######################################################################
######################################################################
##############   R1 norm + F2 norm  realization  #####################
######################################################################
######################################################################
####################
## Data is the sample list, epsilon is the threldhold 
## p*q is the dims of feature, gamma0 not necessary, alpha is the balance term 

CombR.obj <- function(Data, p, q, gamma0, epsilon, init, alpha){  
  timestart<-Sys.time()
  N <- length(Data)  ### sample size
  m <- dim(Data[[1]])[1]  ## number of rows
  n <- dim(Data[[1]])[2]  ## number of columns
  ######initialize U,V
  if ( init == "way1" )  {# way 1
    U <- matrix(0, p, m)
    V <- matrix(0, n, q)
    for(i in 1:p)
      U[i,i] <- 1
    
    for(i in 1:q)
      V[i,i] <- 1  }
  #########initialization with slicing learning i.e. HOSVD
  ## way 2
  else { X.col <- rep(0, n)
  for (j in 1:N)
    X.col <- rbind(X.col, Data[[j]])
  X.col <- X.col[-1,]
  
  G.col <- t(X.col) %*% X.col
  ev.col <- eigen(G.col)
  V <- Re(ev.col$vectors[,1:q])
  
  X.row <- rep(0, m)
  for (j in 1:N)
    X.row <- cbind(X.row, Data[[j]])
  X.row <- X.row[,-1]
  
  G.row <- X.row %*% t(X.row)
  ev.row <- eigen(G.row)
  U <- Re(t(ev.row$vectors[,1:p])) }
  ##################### other initialization
  U.past <- matrix(0, p, m)
  V.past <- matrix(0, n, q)
  
  t <- 1
  d <- rep(0,N)
  s <- 0
  
  s.past <- 0
  for (i in 1:N)
    s.past <- s.past + alpha * norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + 
    norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
  obj.track <- NULL
  UVGap.track <- NULL
  while (round(abs(s - s.past),20 )> epsilon ) {
    ######output summary
    
    print(paste("number of iteration:",t))
    print(paste("past objective value:",s.past) )
    # print(paste("objective gap:",s.past - s  ))
    # cat("\n")
    # 
    s.past <- 0
    for (i in 1:N)
      s.past <- s.past + alpha * norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + 
      norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
    
    obj.track[t] <- s.past
    #####get d
    U.past <- U
    V.past <- V
    for(i in 1:N){
      d[i] <- 1/(norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + gamma0) }
    
    ##### solve optimization problem1
    H <- matrix(0, n, n)
    for(i in 1:N){
      H <- H + t(Data[[i]]) %*% t(U) %*% U %*% Data[[i]] * (alpha * d[i] + 2) }
    ev.H <- eigen(H)
    
    # nor.matrix <- matrix(0,q,q)
    # for (i in 1:q){
    #   nor.matrix[i,i] <- sign(Re( ev.H$vectors[1, i]))
    # }
    #   
    # V <- as.matrix(Re( ev.H$vectors[,1:q]) %*% nor.matrix)  ###update U
    V <- as.matrix(Re( ev.H$vectors[,1:q]))
    ##### solve optimization problem2
    G <- matrix(0, m, m)
    for(i in 1:N){
      G <- G + Data[[i]] %*% V %*% t(V) %*% t(Data[[i]]) * (alpha * d[i] + 2)  }
    ev.G <- eigen(G)
    
    #  nor.matrix2 <- matrix(0,p,p)
    # for (i in 1:p){
    #   nor.matrix2[i,i] <- sign(Re( ev.G$vectors[1, i]))
    # }
    # U <- as.matrix(Re( t(ev.G$vectors[,1:p] %*%  nor.matrix2)))  ####update V
    U <- as.matrix(Re( t(ev.G$vectors[,1:p])))
    #####update object value
    s <- 0
    for (i in 1:N)
      s <- s + alpha * norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + 
      norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
    
    
    UVGap.track[t] <- norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
      norm(V %*% t(V) - V.past %*% t(V.past), type = "F")
    print(paste("now objective value:", s))
    print(paste("objective gap:",round(s.past - s, 20  )))
    print(paste("U,V gap", UVGap.track[t]))
    cat("\n")
    t <- t + 1
    if(t > 10^4) break
  }
  
  timeend<-Sys.time()
  runningtime<-timeend-timestart
  
  result <- list(U = U, V = V, t = t - 1, Obj = obj.track, UV.Gap = UVGap.track,Time = runningtime)
  return(result)
}

######################################################################
##############   R1 norm + F2 norm  realization  #####################
######################     regularization     ########################
######################################################################

CombR.obj.Re <- function(Data, p, q, gamma0, epsilon, init, alpha, init.reg, lambda){  
  timestart<-Sys.time()
  N <- length(Data)  ### sample size
  m <- dim(Data[[1]])[1]  ## number of rows
  n <- dim(Data[[1]])[2]  ## number of columns
  ######initialize U,V
  if ( init == "way1" )  {# way 1
    U <- matrix(0, p, m)
    V <- matrix(0, n, q)
    for(i in 1:p)
      U[i,i] <- 1
    
    for(i in 1:q)
      V[i,i] <- 1  }
  #########initialization with slicing learning i.e. HOSVD
  ## way 2
  else { X.col <- rep(0, n)
  for (j in 1:N)
    X.col <- rbind(X.col, Data[[j]])
  X.col <- X.col[-1,]
  
  G.col <- t(X.col) %*% X.col
  ev.col <- eigen(G.col)
  V <- Re(ev.col$vectors[,1:q])
  
  X.row <- rep(0, m)
  for (j in 1:N)
    X.row <- cbind(X.row, Data[[j]])
  X.row <- X.row[,-1]
  
  G.row <- X.row %*% t(X.row)
  ev.row <- eigen(G.row)
  U <- Re(t(ev.row$vectors[,1:p])) }
  ##################### other initialization
  U.past <- matrix(0, p, m)
  V.past <- matrix(0, n, q)
  
  t <- 1
  d <- rep(0,N)
  s <- 0
  
  s.past <- 0
  for (i in 1:N)
    s.past <- s.past + alpha * norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + 
    norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
  obj.track <- NULL
  UVGap.track <- NULL
  
  while (round(abs(s - s.past),20 )> epsilon ) {
    ######output summary
    print(paste("number of iteration:",t))
    print(paste("past objective value:",s.past) )
    # print(paste("objective gap:",s.past - s  ))
    # cat("\n")
    # 
    if (init.reg == "adaptive"){
      lambda[2] <- 1/norm(t(U) %*% U - t(U.past) %*% U.past, type = "F")
      lambda[1] <- 1/norm(V %*% t(V) - V.past %*% t(V.past), type = "F")
    }
    
    s.past <- 0
    for (i in 1:N)
      s.past <- s.past + alpha * norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + 
      norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
    
    obj.track[t] <- s.past
    #####get d
    U.past <- U
    V.past <- V
    for(i in 1:N){
      d[i] <- 1/(norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + gamma0) }
    
    ##### solve optimization problem1
    H <- - lambda[1] * (diag(n) - V %*% t(V))
    for(i in 1:N){
      H <- H + t(Data[[i]]) %*% t(U) %*% U %*% Data[[i]] * (alpha * d[i] + 2) }
    ev.H <- eigen(H)
    
    # nor.matrix <- matrix(0,q,q)
    # for (i in 1:q){
    #   nor.matrix[i,i] <- sign(Re( ev.H$vectors[1, i]))
    # }
    #   
    # V <- as.matrix(Re( ev.H$vectors[,1:q]) %*% nor.matrix)  ###update U
    V <- as.matrix(Re( ev.H$vectors[,1:q]))
    ##### solve optimization problem2
    G <- - lambda[2] * (diag(m) - t(U) %*% U)
    for(i in 1:N){
      G <- G + Data[[i]] %*% V %*% t(V) %*% t(Data[[i]]) * (alpha * d[i] + 2)  }
    ev.G <- eigen(G)
    
    #  nor.matrix2 <- matrix(0,p,p)
    # for (i in 1:p){
    #   nor.matrix2[i,i] <- sign(Re( ev.G$vectors[1, i]))
    # }
    # U <- as.matrix(Re( t(ev.G$vectors[,1:p] %*%  nor.matrix2)))  ####update V
    U <- as.matrix(Re( t(ev.G$vectors[,1:p])))
    #####update object value
    s <- 0
    for (i in 1:N)
      s <- s + alpha * norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + 
      norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
    
    
    UVGap.track[t] <- norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
      norm(V %*% t(V) - V.past %*% t(V.past), type = "F")
    print(paste("now objective value:", s))
    print(paste("objective gap:",round(s.past - s, 20  )))
    print(paste("U,V gap", UVGap.track[t]))
    cat("\n")
    t <- t + 1
    if(t > 10^4) break
  }
  timeend<-Sys.time()
  runningtime<-timeend-timestart
  
  result <- list(U = U, V = V, t = t - 1, Obj = obj.track, UV.Gap = UVGap.track, Time = runningtime)
  return(result)
}

######################################################################
######################################################################
##############   Tucker decomposition   ##############################
######################################################################
######################################################################
####################
## Data is the sample list, epsilon is the threldhold 
## p*q is the dims of feature, init is the way of initialization

TuckerALS.obj <- function(Data, p, q, epsilon, init ){  
  timestart<-Sys.time()
  N <- length(Data)  ### sample size
  m <- dim(Data[[1]])[1]  ## number of rows
  n <- dim(Data[[1]])[2]  ## number of columns
  ######initialize U,V
  if ( init == "way1" )  {# way 1
    U <- matrix(0, p, m)
    V <- matrix(0, n, q)
    for(i in 1:p)
      U[i,i] <- 1
    
    for(i in 1:q)
      V[i,i] <- 1  }
  #########initialization with slicing learning i.e. HOSVD
  ## way 2
  else { X.col <- rep(0, n)
  for (j in 1:N)
    X.col <- rbind(X.col, Data[[j]])
  X.col <- X.col[-1,]
  
  G.col <- t(X.col) %*% X.col
  ev.col <- eigen(G.col)
  V <- Re(ev.col$vectors[,1:q])
  
  X.row <- rep(0, m)
  for (j in 1:N)
    X.row <- cbind(X.row, Data[[j]])
  X.row <- X.row[,-1]
  
  G.row <- X.row %*% t(X.row)
  ev.row <- eigen(G.row)
  U <- Re(t(ev.row$vectors[,1:p])) }
  
  ##################### other initialization
  U.past <- matrix(0, p, m)
  V.past <- matrix(0, n, q)
  
  t <- 1
  # d <- rep(0,N)
  s <- 0
  
  s.past <- 0
  for (i in 1:N)
    s.past <- s.past + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
  
  obj.track <- NULL
  UVGap.track <- NULL
  while (round(abs(s - s.past),20 )> epsilon ) {
    ######output summary
    
    print(paste("number of iteration:",t))
    print(paste("past objective value:",s.past) )
    # print(paste("objective gap:",s.past - s  ))
    # cat("\n")
    # 
    s.past <- 0
    for (i in 1:N)
      s.past <- s.past + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
    
    obj.track[t] <- s.past
    U.past <- U
    V.past <- V
    # #####get d
    # for(i in 1:N){
    #   U.past <- U
    #   V.past <- V
    #   d[i] <- 1/(norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + gamma0) }
    
    ##### solve optimization problem1
    H <- matrix(0, n, n)
    for(i in 1:N){
      H <- H + t(Data[[i]]) %*% t(U) %*% U %*% Data[[i]]}
    ev.H <- eigen(H)
    
    # nor.matrix <- matrix(0,q,q)
    # for (i in 1:q){
    #   nor.matrix[i,i] <- sign(Re( ev.H$vectors[1, i]))
    # }
    #   
    # V <- as.matrix(Re( ev.H$vectors[,1:q]) %*% nor.matrix)  ###update U
    V <- as.matrix(Re( ev.H$vectors[,1:q]))
    ##### solve optimization problem2
    G <- matrix(0, m, m)
    for(i in 1:N){
      G <- G + Data[[i]] %*% V %*% t(V) %*% t(Data[[i]])}
    ev.G <- eigen(G)
    
    #  nor.matrix2 <- matrix(0,p,p)
    # for (i in 1:p){
    #   nor.matrix2[i,i] <- sign(Re( ev.G$vectors[1, i]))
    # }
    # U <- as.matrix(Re( t(ev.G$vectors[,1:p] %*%  nor.matrix2)))  ####update V
    U <- as.matrix(Re( t(ev.G$vectors[,1:p])))
    #####update object value
    s <- 0
    for (i in 1:N)
      s <- s + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
    
    
    UVGap.track[t] <- norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
      norm(V %*% t(V) - V.past %*% t(V.past), type = "F")
    print(paste("now objective value:", s))
    print(paste("objective gap:",round(s.past - s, 20  )))
    print(paste("U,V gap",UVGap.track[t]) )
    cat("\n")
    
    t <- t + 1
    if(t > 2*10^3) break
  }
  
  timeend<-Sys.time()
  runningtime<-timeend-timestart
  
  result <- list(U = U, V = V, t = t - 1, Obj = obj.track, UV.Gap = UVGap.track, Time = runningtime)
  return(result)
}



######################################################################
######################################################################
#####    Tucker decomposition with regularization  ###################
######################################################################
######################################################################
####################
##Data is the sample list, epsilon is the threldhold 
## p*q is the dims of feature, lambda is the initial regularization weight

TuckerALS_Re.obj <- function(Data, p, q, epsilon, init, lambda, U_initial, V_initial ){  ##Data is the sample list, epsilon is the threldhold 
  ## p*q is the dims of feature, gamma0 not necessary
  timestart<-Sys.time()
  N <- length(Data)  ### sample size
  m <- dim(Data[[1]])[1]  ## number of rows
  n <- dim(Data[[1]])[2]  ## number of columns
  ######initialize U,V
  if ( init == "way1" )  {# way 1
    U <- matrix(0, p, m)
    V <- matrix(0, n, q)
    for(i in 1:p)
      U[i,i] <- 1
    
    for(i in 1:q)
      V[i,i] <- 1  }
  #########initialization with slicing learning i.e. HOSVD
  ## way 2
  else if (init == "way2") { X.col <- rep(0, n)
  for (j in 1:N)
    X.col <- rbind(X.col, Data[[j]])
  X.col <- X.col[-1,]
  
  G.col <- t(X.col) %*% X.col
  ev.col <- eigen(G.col)
  V <- Re(ev.col$vectors[,1:q])
  
  X.row <- rep(0, m)
  for (j in 1:N)
    X.row <- cbind(X.row, Data[[j]])
  X.row <- X.row[,-1]
  
  G.row <- X.row %*% t(X.row)
  ev.row <- eigen(G.row)
  U <- Re(t(ev.row$vectors[,1:p])) }
  
  else {U <- U_initial
  V <- V_initial}
  
  ##################### other initialization
  U.past <- matrix(0, p, m)
  V.past <- matrix(0, n, q)
  
  t <- 1
  # d <- rep(0,N)
  s <- 0
  
  s.past <- 0
  for (i in 1:N)
    s.past <- s.past + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
  
  obj.track <- NULL
  UVGap.track <- NULL
  while (round(abs(s - s.past),20 )> epsilon ) {
    ######output summary
    
    print(paste("number of iteration:",t))
    print(paste("past objective value:",s.past) )
    # print(paste("objective gap:",s.past - s  ))
    # cat("\n")
    # 
    s.past <- 0
    for (i in 1:N)
      s.past <- s.past + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
    
    obj.track[t] <- s.past
    U.past <- U
    V.past <- V
    # #####get d
    # for(i in 1:N){
    #   U.past <- U
    #   V.past <- V
    #   d[i] <- 1/(norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + gamma0) }
    
    ##### solve optimization problem1
    H <- - lambda[1] * (diag(n) - V %*% t(V))
    for(i in 1:N){
      H <- H + t(Data[[i]]) %*% t(U) %*% U %*% Data[[i]]}
    ev.H <- eigen(H)
    V <- as.matrix(Re( ev.H$vectors[,1:q]))
    # L1 <- length(which(ev.H$values>=0))
    # if ( L1 >= q ) {V <- as.matrix(Re( ev.H$vectors[,1:q]))}
    # else  V <- as.matrix(cbind(Re( ev.H$vectors[,1:L1]), Re( ev.H$vectors[, (n - L1 + 1):n]) )) 
    # nor.matrix <- matrix(0,q,q)
    # for (i in 1:q){
    #   nor.matrix[i,i] <- sign(Re( ev.H$vectors[1, i]))
    # }
    #   
    # V <- as.matrix(Re( ev.H$vectors[,1:q]) %*% nor.matrix)  ###update U
    
    ##### solve optimization problem2
    G <- - lambda[2]  * (diag(m) - t(U) %*% U)
    for(i in 1:N){
      G <- G + Data[[i]] %*% V %*% t(V) %*% t(Data[[i]])}
    ev.G <- eigen(G)
    
    # L2 <- length(which(ev.G$values>=0))
    # if ( L2 >= p ) {U <- as.matrix(Re( t(ev.G$vectors[,1:p])))}
    # else  U <- as.matrix(rbind(Re( t(ev.G$vectors[,1:L2])), Re( t(ev.G$vectors[, m + 1 - L1:m])) )) 
    
    #  nor.matrix2 <- matrix(0,p,p)
    # for (i in 1:p){
    #   nor.matrix2[i,i] <- sign(Re( ev.G$vectors[1, i]))
    # }
    # U <- as.matrix(Re( t(ev.G$vectors[,1:p] %*%  nor.matrix2)))  ####update V
    U <- as.matrix(Re( t(ev.G$vectors[,1:p])))
    #####update object value
    s <- 0
    for (i in 1:N)
      s <- s + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
    
    UVGap.track[t] <- norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
      norm(V %*% t(V) - V.past %*% t(V.past), type = "F")
    
    print(paste("now objective value:", s))
    print(paste("objective gap:",round(s.past - s, 20  )))
    print(paste("U,V gap",UVGap.track[t]))
    cat("\n")
    t <- t + 1
    if(t > 10^4) break
  }
  
  timeend<-Sys.time()
  runningtime<-timeend-timestart
  
  result <- list(U = U, V = V, t = t - 1, Obj = obj.track, UV.Gap = UVGap.track, Time = runningtime)
  return(result)
}


######################################################################
######################################################################
#####    Tucker decomposition with regularization  ###################
#####         wtih lambda=1/|U_tU_t^T-XXXXX|       ###################
######################################################################
######################################################################
##Data is the sample list, epsilon is the threldhold 
## p*q is the dims of feature, lambda is the initial regularization weight

TuckerALS_Re2.obj <- function(Data, p, q, epsilon, init, lambda, U_initial, V_initial ){  ##Data is the sample list, epsilon is the threldhold 
  ## p*q is the dims of feature, gamma0 not necessary
  timestart<-Sys.time()
  N <- length(Data)  ### sample size
  m <- dim(Data[[1]])[1]  ## number of rows
  n <- dim(Data[[1]])[2]  ## number of columns
  ######initialize U,V
  if ( init == "way1" )  {# way 1
    U <- matrix(0, p, m)
    V <- matrix(0, n, q)
    for(i in 1:p)
      U[i,i] <- 1
    
    for(i in 1:q)
      V[i,i] <- 1  }
  #########initialization with slicing learning i.e. HOSVD
  ## way 2
  else if (init == "way2") { X.col <- rep(0, n)
  for (j in 1:N)
    X.col <- rbind(X.col, Data[[j]])
  X.col <- X.col[-1,]
  
  G.col <- t(X.col) %*% X.col
  ev.col <- eigen(G.col)
  V <- Re(ev.col$vectors[,1:q])
  
  X.row <- rep(0, m)
  for (j in 1:N)
    X.row <- cbind(X.row, Data[[j]])
  X.row <- X.row[,-1]
  
  G.row <- X.row %*% t(X.row)
  ev.row <- eigen(G.row)
  U <- Re(t(ev.row$vectors[,1:p])) }
  
  else {U <- U_initial
  V <- V_initial}
  
  ##################### other initialization
  U.past <- matrix(0, p, m)
  V.past <- matrix(0, n, q)
  
  t <- 1
  # d <- rep(0,N)
  s <- 0
  
  s.past <- 0
  for (i in 1:N)
    s.past <- s.past + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
  
  obj.track <- NULL
  UVGap.track <- NULL
  while (round(abs(s - s.past),20 )> epsilon ) {
    ######output summary
    
    print(paste("number of iteration:",t))
    print(paste("past objective value:",s.past) )
    # print(paste("objective gap:",s.past - s  ))
    # cat("\n")
    # 
    lambda[2] <- 1/norm(t(U) %*% U - t(U.past) %*% U.past, type = "F")
    lambda[1] <- 1/norm(V %*% t(V) - V.past %*% t(V.past), type = "F")
    
    s.past <- 0
    for (i in 1:N)
      s.past <- s.past + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
    
    obj.track[t] <- s.past
    U.past <- U
    V.past <- V
    # #####get d
    # for(i in 1:N){
    #   U.past <- U
    #   V.past <- V
    #   d[i] <- 1/(norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F") + gamma0) }
    
    ##### solve optimization problem1
    H <- - lambda[1] * (diag(n) - V %*% t(V))
    for(i in 1:N){
      H <- H + t(Data[[i]]) %*% t(U) %*% U %*% Data[[i]]}
    ev.H <- eigen(H)
    V <- as.matrix(Re( ev.H$vectors[,1:q]))
    # L1 <- length(which(ev.H$values>=0))
    # if ( L1 >= q ) {V <- as.matrix(Re( ev.H$vectors[,1:q]))}
    # else  V <- as.matrix(cbind(Re( ev.H$vectors[,1:L1]), Re( ev.H$vectors[, (n - L1 + 1):n]) )) 
    # nor.matrix <- matrix(0,q,q)
    # for (i in 1:q){
    #   nor.matrix[i,i] <- sign(Re( ev.H$vectors[1, i]))
    # }
    #   
    # V <- as.matrix(Re( ev.H$vectors[,1:q]) %*% nor.matrix)  ###update U
    
    ##### solve optimization problem2
    G <- - lambda[2] * (diag(m) - t(U) %*% U)
    for(i in 1:N){
      G <- G + Data[[i]] %*% V %*% t(V) %*% t(Data[[i]])}
    ev.G <- eigen(G)
    
    # L2 <- length(which(ev.G$values>=0))
    # if ( L2 >= p ) {U <- as.matrix(Re( t(ev.G$vectors[,1:p])))}
    # else  U <- as.matrix(rbind(Re( t(ev.G$vectors[,1:L2])), Re( t(ev.G$vectors[, m + 1 - L1:m])) )) 
    
    #  nor.matrix2 <- matrix(0,p,p)
    # for (i in 1:p){
    #   nor.matrix2[i,i] <- sign(Re( ev.G$vectors[1, i]))
    # }
    # U <- as.matrix(Re( t(ev.G$vectors[,1:p] %*%  nor.matrix2)))  ####update V
    U <- as.matrix(Re( t(ev.G$vectors[,1:p])))
    #####update object value
    s <- 0
    for (i in 1:N)
      s <- s + norm(Data[[i]] - t(U) %*% U %*% Data[[i]] %*% V %*% t(V), type = "F")^2
    
    
    UVGap.track[t] <- norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
      norm(V %*% t(V) - V.past %*% t(V.past), type = "F")
    
    print(paste("now objective value:", s))
    print(paste("objective gap:",round(s.past - s, 20  )))
    print(paste("U,V gap",UVGap.track[t]))
    cat("\n")
    
    t <- t + 1
    if(t > 10^4) break
  }
  
  timeend<-Sys.time()
  runningtime<-timeend-timestart
  
  result <- list(U = U, V = V, t = t - 1, Obj = obj.track, UV.Gap = UVGap.track, Time = runningtime)
  return(result)
}



##########################
####################################################################
################# Numerical study 1: Objective value w/ or w/o noise
####################################################################
#####################################################
########### w/ noise is good
#####################################################
library(rmutil)
path <- ("C:/Users/boshen.HOKIES/Desktop/image")
pathname <- file.path(path, "labimage.mat")
Dat <- readMat(pathname)
data1 <- Dat$dat.raw
Data.raw <- NULL
for (i in 1:655) {
  Data.raw[[i]] <- matrix(data1[i,], 80, 80)/255 
}

set.seed(20)
###### train.index <- sort(sample(1:655,300))
train.index <- 1:300
Train.raw <- Data.raw[train.index ]
Test.raw <- Data.raw[c(1:655)[-train.index ]]

###### add pure noise in the training step
set.seed(10)    
#####Case1: 300 samples + 25 noise samples has ok result with c(10,10) 
#####Case2: 300 samples + 20 noise samples has good result with c(10,10), ok result with c(20,20)
#####Case3: 300 samples + 10 noise samples has ok result c(5,5), not good result with c(10,10)
#####Case4: 300 samples + 19 noise samples has ok result with c(10,10)
#####Case5: 280 samples + 20 noise samples has ok result with c(10,10) similar with case 3:c(5,5)
   for (j in 301:320)       #for (j in 281:300)  
  Train.raw[[j]]  <- matrix(runif(80*80), 80, 80)


set.seed(205)
######HOOI
Tucker <- TuckerALS(Train.raw, 30, 30, 10^-6, "way1")
######lambda_n is constant c(10,10) is ag good choice for 300 + 20 noise
Tucker.Re1 <- TuckerALS_Re(Train.raw,30,30, 10^-6, "way1", c(10,10), Tucker$U, Tucker$V)
######lambda_n wil change over time based on the last step
Tucker.Re2 <- TuckerALS_Re2(Train.raw,30,30, 10^-6, "way1", c(10,10), Tucker$U, Tucker$V)
# PC <- BFPCA(Train.raw, 30, 30,0, 10^-5, "way1")
Tucker$t
Tucker.Re1$t
Tucker.Re2$t

Tucker$Time
Tucker.Re1$Time
Tucker.Re2$Time

burn.in <- 5
Tucker.Re2$Obj[length(Tucker.Re2$Obj)]
Tucker.Re1$Obj[length(Tucker.Re1$Obj)]
Tucker$Obj[length(Tucker$Obj)]

Tu <- Tucker$UV.Gap[(burn.in + 1):Tucker$t]
Tu.Re1 <- Tucker.Re1$UV.Gap[(burn.in + 1):Tucker.Re1$t]
Tu.Re2 <- Tucker.Re2$UV.Gap[(burn.in + 1):Tucker.Re2$t]

sample <- c(rep("literature",Tucker$t - burn.in),rep("C-Regulatization", Tucker.Re1$t - burn.in),rep("AdaptiveRegulatization", Tucker.Re2$t - burn.in))
dims <- c((burn.in + 1):Tucker$t,(burn.in + 1):Tucker.Re1$t,(burn.in + 1):Tucker.Re2$t)
y <- c(Tu, Tu.Re1, Tu.Re2)
data_graph <- data.frame(dims, sample, y)
length(sample)
length(dims)

p <- ggplot(data=data_graph,aes(x=dims,y= y,group=sample)) + geom_line(aes(colour=sample)) +
  geom_point(size=1.1,aes(shape=sample,colour=sample)) + xlab("")+ylab("UV.Gap") + 
  scale_y_log10() + theme(panel.background = element_rect( colour = "black", size = 1))
                                                       
p
dev.off()


#####################################################
########### w/O noise is not good maybe try different rank
#####################################################
Data.raw <- NULL
for (i in 1:655) {
  Data.raw[[i]] <- matrix(data1[i,], 80, 80)/255 
}

set.seed(20)
train.index <- 1:300
Train.raw <- Data.raw[train.index ]
Test.raw <- Data.raw[c(1:655)[-train.index ]]

set.seed(205)
######HOOI
Tucker <- TuckerALS(Train.raw, 30, 30, 10^-6, "way1")
######lambda_n is constant
Tucker.Re1 <- TuckerALS_Re(Train.raw,30,30, 10^-6, "way1", c(20,20), Tucker$U, Tucker$V)  ##c(5,5)
######lambda_n wil change over time based on the last step
Tucker.Re2 <- TuckerALS_Re2(Train.raw,30,30, 10^-6, "way1", c(20,20), Tucker$U, Tucker$V)
# PC <- BFPCA(Train.raw, 30, 30,0, 10^-5, "way1")
Tucker$t
Tucker.Re1$t
Tucker.Re2$t
burn.in <- 50
Tucker.Re2$Obj[length(Tucker.Re2$Obj)]
Tucker.Re1$Obj[length(Tucker.Re1$Obj)]
Tucker$Obj[length(Tucker$Obj)]

Tu <- Tucker$UV.Gap[(burn.in + 1):Tucker$t]
Tu.Re1 <- Tucker.Re1$UV.Gap[(burn.in + 1):Tucker.Re1$t]
Tu.Re2 <- Tucker.Re2$UV.Gap[(burn.in + 1):Tucker.Re2$t]

sample <- c(rep("literature",Tucker$t - burn.in),rep("C-Regulatization", Tucker.Re1$t - burn.in),rep("AdaptiveRegulatization", Tucker.Re2$t - burn.in))
dims <- c((burn.in + 1):Tucker$t,(burn.in + 1):Tucker.Re1$t,(burn.in + 1):Tucker.Re2$t)
y <- c(Tu, Tu.Re1, Tu.Re2)
data_graph <- data.frame(dims, sample, y)
length(sample)
length(dims)

p <- ggplot(data=data_graph,aes(x=dims,y= y,group=sample)) + geom_line(aes(colour=sample)) +
  geom_point(size=1.1,aes(shape=sample,colour=sample)) + xlab("")+ylab("UV.Gap")+ 
  scale_y_log10() + theme(panel.background = element_rect( colour = "black", size = 2))
p
dev.off()

########################## This time you should run the BFPCA.R file
####################################################################
################# Numerical study 2: UV gap w/ or w/o noise
####################################################################

#####################################################
########### w/ noise is good 
#####################################################
Data.raw <- NULL
for (i in 1:655) {
  Data.raw[[i]] <- matrix(data1[i,], 80, 80)/255 
}

set.seed(20)
train.index <- 1:300
Train.raw <- Data.raw[train.index ]
Test.raw <- Data.raw[c(1:655)[-train.index ]]

###### add pure noise in the training step
set.seed(10)
for (j in 301:320)   #######1.500 2.700 
  Train.raw[[j]]  <- matrix(runif(80*80), 80, 80)


set.seed(205)
######HOOI
Tucker <- TuckerALS(Train.raw, 30, 30, 10^-5, "way1")
######lambda_n is constant
Tucker.Re1 <- TuckerALS_Re(Train.raw,30,30, 10^-5, "way1", c(20,20), Tucker$U, Tucker$V)
######lambda_n wil change over time based on the last step
Tucker.Re2 <- TuckerALS_Re2(Train.raw,30,30, 10^-5, "way1", c(20,20), Tucker$U, Tucker$V)
# PC <- BFPCA(Train.raw, 30, 30,0, 10^-5, "way1")
Tucker$t
Tucker.Re1$t
Tucker.Re2$t
burn.in <- 20
Tucker.Re2$Obj[length(Tucker.Re2$Obj)]
Tucker.Re1$Obj[length(Tucker.Re1$Obj)]
Tucker$Obj[length(Tucker$Obj)]

Tu <- Tucker$UV.Gap[(burn.in + 1):Tucker$t]
Tu.Re1 <- Tucker.Re1$UV.Gap[(burn.in + 1):Tucker.Re1$t]
Tu.Re2 <- Tucker.Re2$UV.Gap[(burn.in + 1):Tucker.Re2$t]

sample <- c(rep("literature",Tucker$t - burn.in),rep("C-Regulatization", Tucker.Re1$t - burn.in),rep("AdaptiveRegulatization", Tucker.Re2$t - burn.in))
dims <- c((burn.in + 1):Tucker$t,(burn.in + 1):Tucker.Re1$t,(burn.in + 1):Tucker.Re2$t)
y <- c(Tu, Tu.Re1, Tu.Re2)
data_graph <- data.frame(dims, sample, y)
length(sample)
length(dims)

p <- ggplot(data=data_graph,aes(x=dims,y= y,group=sample)) + geom_line(aes(colour=sample)) +
  geom_point(size=1.1,aes(shape=sample,colour=sample)) + xlab("")+ylab("UV.Gap")+ 
  scale_y_log10() + theme(panel.background = element_rect( colour = "black", size = 2))
p
dev.off()



#####################################################
########### w/O noise is not good
#####################################################
Data.raw <- NULL
for (i in 1:655) {
  Data.raw[[i]] <- matrix(data1[i,], 80, 80)/255 
}

set.seed(20)
train.index <- 1:300
Train.raw <- Data.raw[train.index ]
Test.raw <- Data.raw[c(1:655)[-train.index ]]

set.seed(205)
######HOOI
Tucker <- TuckerALS(Train.raw, 30, 30, 10^-5, "way1")
######lambda_n is constant
Tucker.Re1 <- TuckerALS_Re(Train.raw,30,30, 10^-5, "way1", c(20,20), Tucker$U, Tucker$V)
######lambda_n wil change over time based on the last step
Tucker.Re2 <- TuckerALS_Re2(Train.raw,30,30, 10^-5, "way1", c(20,20), Tucker$U, Tucker$V)
# PC <- BFPCA(Train.raw, 30, 30,0, 10^-5, "way1")
Tucker$t
Tucker.Re1$t
Tucker.Re2$t
burn.in <- 20
Tucker.Re2$Obj[length(Tucker.Re2$Obj)]
Tucker.Re1$Obj[length(Tucker.Re1$Obj)]
Tucker$Obj[length(Tucker$Obj)]

Tu <- Tucker$UV.Gap[(burn.in + 1):Tucker$t]
Tu.Re1 <- Tucker.Re1$UV.Gap[(burn.in + 1):Tucker.Re1$t]
Tu.Re2 <- Tucker.Re2$UV.Gap[(burn.in + 1):Tucker.Re2$t]

sample <- c(rep("literature",Tucker$t - burn.in),rep("C-Regulatization", Tucker.Re1$t - burn.in),rep("AdaptiveRegulatization", Tucker.Re2$t - burn.in))
dims <- c((burn.in + 1):Tucker$t,(burn.in + 1):Tucker.Re1$t,(burn.in + 1):Tucker.Re2$t)
y <- c(Tu, Tu.Re1, Tu.Re2)
data_graph <- data.frame(dims, sample, y)
length(sample)
length(dims)

p <- ggplot(data=data_graph,aes(x=dims,y= y,group=sample)) + geom_line(aes(colour=sample)) +
  geom_point(size=1.1,aes(shape=sample,colour=sample)) + xlab("")+ylab("UV.Gap")+ 
  scale_y_log10() + theme(panel.background = element_rect( colour = "black", size = 2))
p
dev.off()
