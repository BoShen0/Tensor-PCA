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
#################UV gap as stopping criteria  ########################
######################################################################

######################################################################
######################################################################
######################   R1 norm realization  ########################
######################################################################
######################################################################
####################
BFPCA <- function(Data, p, q, gamma0, epsilon, init ){  ##Data is the sample list, epsilon is the threldhold 
  timestart<-Sys.time()                                                ## p*q is the dims of feature, gamma0 not necessary
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
  while (round(norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
               norm(V %*% t(V) - V.past %*% t(V.past), type = "F"),20 )> epsilon ) { ## abs(s.past - s)
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
  
  result <- list(U = U, V = V, t = t - 1, Obj = obj.track, UV.Gap = UVGap.track, Time = runningtime)
  return(result)
}

######################################################################
######################################################################
######################   R1 norm realization  ########################
######################     regularization     ########################
######################################################################
######################################################################
BFPCA.Re <- function(Data, p, q, gamma0, epsilon, init, init.reg, lambda){  ##Data is the sample list, epsilon is the threldhold 
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
  while (round(norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
               norm(V %*% t(V) - V.past %*% t(V.past), type = "F"),20 )> epsilon ) { ## abs(s.past - s)
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

CombR <- function(Data, p, q, gamma0, epsilon, init, alpha){  
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
  while (round(norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
               norm(V %*% t(V) - V.past %*% t(V.past), type = "F"),20 )> epsilon ) {
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
  
  result <- list(U = U, V = V, t = t - 1, Obj = obj.track, UV.Gap = UVGap.track, Time = runningtime)
  return(result)
}

######################################################################
##############   R1 norm + F2 norm  realization  #####################
######################     regularization     ########################
######################################################################

CombR.Re <- function(Data, p, q, gamma0, epsilon, init, alpha, init.reg, lambda){  
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
  while (round(norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
               norm(V %*% t(V) - V.past %*% t(V.past), type = "F"),20 )> epsilon ) {
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
##Data is the sample list, epsilon is the threldhold 
## p*q is the dims of feature, init is the way of initialization

TuckerALS <- function(Data, p, q, epsilon, init ){  
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
  while (round(norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
               norm(V %*% t(V) - V.past %*% t(V.past), type = "F"),20 )> epsilon ) {
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

TuckerALS_Re <- function(Data, p, q, epsilon, init, lambda, U_initial, V_initial ){  ##Data is the sample list, epsilon is the threldhold 
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
  while (round(norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
               norm(V %*% t(V) - V.past %*% t(V.past), type = "F"),20 )> epsilon ) {
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

TuckerALS_Re2 <- function(Data, p, q, epsilon, init, lambda, U_initial, V_initial ){  ##Data is the sample list, epsilon is the threldhold 
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
  while (round(norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
               norm(V %*% t(V) - V.past %*% t(V.past), type = "F"),20 )> epsilon ) {
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


# set.seed(10)
# A= NULL
# for (j in 1:200){
#   A[[j]] <- matrix(rnorm(100*100),100,100)
# }
# result1 <- BFPCA(A, 2, 1, 0, 0.000001)
##################### 
########################################################################
#####  EXPERIEMNT1: numeriacal result-ORL face data 
########################################################################
path <- "C:/Users/boshen.HOKIES/Desktop/take home/Research/CCA_coding new/Case study/att_faces"
fileNames <- dir(path)  
sub.path <- "C:/Users/boshen.HOKIES/Desktop/take home/Research/CCA_coding new/Case study/att_faces/s2"
sub.fileNames <- dir(sub.path)
filePath <- data.frame(matrix(0, length(fileNames), 10))
for(i in 1:40)
  for (j in 1:10)
    filePath[i,j] <- paste(path,fileNames[i],sub.fileNames[j],sep='/')

############################
label <- rep(1:40, rep(10, 40))
raw.data <- NULL
for (i in 1:40)
  for(j in 1:10){
    Y <- read.pnm(filePath[i, j])
    raw.data[[10 * (i - 1) + j]] <- as.matrix(Y@grey)
  }

################ - the mean
# Mean <- matrix(0, 112, 92)
# for (j in 1:400)
# Mean <- Mean + 1/400 * raw.data[[j]]
# 
# for (j in 1:400)
#   raw.data[[j]] <- raw.data[[j]] - Mean
set.seed(10)
for (j in 401:420)
raw.data[[j]] <- matrix(runif(112*92), 112, 92)

PC <- BFPCA(raw.data[c(1:30,401,402)],25,25,0, 10^-11, "way1")
Tucker <- TuckerALS(raw.data[c(1:30,401,402)],25,25, 10^-11, "way1")

TT <- CombR(raw.data[c(1:30,401,402)],25,25,0, 10^-11, "way1", 10 )
TTT <- TuckerALS_Re(raw.data,30,30, 10^-11, "way3", c(0.4,0.4), Tucker$U, Tucker$V)

########################################################################
#####          compare with RMSE
########################################################################

RmseF2 <-NULL
RmseR1<- NULL
for (k in 10:30){
PC <- BFPCA(raw.data[c(1:30)],k,k,0, 10^-11, "way1")
Tucker <- TuckerALS(raw.data[c(1:30)],k,k, 10^-11, "way1")

U <- Tucker$U
V <- Tucker$V
z <- 0
for (j in 31:400)
{
  z <-  z + norm(raw.data[[i]] - t(U) %*% U %*% raw.data[[i]] %*% V %*% t(V), type = "F")
}
RmseF2[k-9] <- z/370

U <- PC$U
V <- PC$V
z <- 0
for (j in 31:400)
{
  z <-  z + norm(raw.data[[i]] - t(U) %*% U %*% raw.data[[i]] %*% V %*% t(V), type = "F")
}
RmseR1[k-9] <- z/370

}

# plot(RmseF2,type = "l",ylim=c(0,100),col="red")
# lines(RmseR1,col="blue" )

sample <- c(rep("lierature",21),rep("Mine",21))
dims <- rep(10:30, 2)
y <- c(RmseF2,RmseR1)
data_graph <- data.frame(dims, sample, y)

 
p <- ggplot(data=data_graph,aes(x=dims,y= y,group=sample)) + geom_line(aes(colour=sample)) +
  geom_point(size=2,aes(shape=sample,colour=sample)) + xlab("")+ylab("Reconstruction Error")
p
dev.off()



####################################################################
################# EXPERIEMNT2: image reconstruction
####################################################################
par(mfrow=c(1,3))

U <- Tucker$U
V <- Tucker$V
image(t(t(U) %*% U %*% raw.data[[40]] %*% V %*% t(V)),col=grey.colors(255), xlab = "Literature")

U <- PC$U
V <- PC$V
image(t(t(U) %*% U %*% raw.data[[40]] %*% V %*% t(V)),col=grey.colors(255))

U <- TT$U
V <- TT$V
image(t(t(U) %*% U %*% raw.data[[40]] %*% V %*% t(V)),col=grey.colors(255))

########################################################################
##### EXPERIEMNT3: classification result-ORL face data 
########################################################################

number.repeatition <- 10
train.index <- matrix(0,number.repeatition, 40 * 5) 
#########generate the random sample
set.seed(31)
for(k in 1:number.repeatition)
{ xx <- rep(0,200)
for (j in 1:40) 
{ 
  xx[((j-1)*5 + 1):(5 * j)] <- 10 * (j - 1) + sample(1:10,5) 
  train.index[k,] <- sort(xx)
}
}

label.training <- rep(1:40, rep(5, 40))
label.test <- rep(1:40, rep(5, 40))


accuracy.ORLknn <- rep(0, number.repeatition)
precision.ORLknn <- matrix(0, number.repeatition, 40 )
recall.ORLknn <- matrix(0, number.repeatition, 40 )

for (k in 1:number.repeatition){
  
  ORL.face.train <- raw.data[c(train.index[k,], 401:420)]  #### get train, test data
  ORL.face.test <- raw.data[c(1:400)[-train.index[k,]]]
  
 
  Tucker <- TuckerALS(ORL.face.train, 30, 30, 10^-11, "way1")
    U <- Tucker$U
    V <- Tucker$V
  # PC <- BFPCA(ORL.face.train, 30, 30, 0, 10^-11, "way1")
  # U <- PC$U
  # V <- PC$V

  train_vars <- matrix(0,200,30*30)
  test_vars <- matrix(0,200,30*30)
  for(j in 1:200) {
    train_vars[j,] <- U %*% ORL.face.train[[j]] %*% V
    test_vars[j,] <- U %*% ORL.face.test[[j]] %*% V
  }
  # 
  ###########kNN based on the canonical variates for training and test
  set.seed(100)
  figure.knn <- knn(train_vars, test_vars, label.training,k=1,prob=T,use.all=T)  ###k=? is good question
  accuracy.ORLknn[k] <- 1 - classError(figure.knn,label.test )$errorRate
  predict.label <- as.numeric(figure.knn[1:length(label.test)])
  
  ###########get calculation results
  for (j in 1:40)
  {
    precision.ORLknn[k, j ] <- length(which(predict.label[which(label.test == j)] == j ))/length(which(predict.label == j)) 
    recall.ORLknn[k, j ] <- length(which(predict.label[which(label.test == j)] == j ))/length(which(label.test == j)) 
  }
  
  print(k)
}  #the end 



########## calculate the results
mean(accuracy.ORLknn)
mean(precision.ORLknn)
mean(recall.ORLknn)
mean(2 * precision.ORLknn * recall.ORLknn/ (precision.ORLknn + recall.ORLknn) )

# [1] "number of iteration: 6"
# [1] "U,V gap 5.97343536443577e-06"
# [1] "past objective value: 10139.9532912925"
# [1] "now objective value: 10139.9532912923"
# [1] "objective gap: -5.45696821e-12"

# ####### manually add noise to verify the algorithm
# #####add noise
# image_noise <- raw.data
# zz <- image_noise[[1]]
# for (i in 1:30)
# zz[20+i, 31:50] <- rbinom(20, 1, 0.5)
# image_noise[[1]] <- zz
# # image_noise[20,sample(c(1:92),size=10)] <- rbinom(10, 1, 0.5)
# image(image_noise[[1]],col=grey.colors(255))
# 
# PC <- BFPCA(image_noise,80,80,0, 10^-12)
# dim(PC$U)
# U <- PC$U[1:50,]
# V <- PC$V[,1:50]
# image(t(U) %*% U %*% image_noise[[1]] %*% V %*% t(V),col=grey.colors(255))
# 
# ########
# tt <- matrix(0,112,92)
# for (j in 1:92)
# tt[,j ] <- rbinom(112, 1, 0.5)
# image(tt,col=grey.colors(255))
# image_noise2 <- raw.data
# image_noise2[[401]] <- tt
# 
# image(image_noise2[[401]],col=grey.colors(255))
# 
# PC <- BFPCA(image_noise2,80,80,0, 10^-12)
# dim(PC$U)
# U <- PC$U[1:80,]
# V <- PC$V[,1:80]
# image(image_noise2[[1]],col=grey.colors(255))
# image(t(U) %*% U %*% image_noise2[[1]] %*% V %*% t(V),col=grey.colors(255))
# 
# 
# ######################
# ###################### tensor slice learning i.e. HOSVD
# X.col <- rep(0, 92)
#  for (j in 1:400)
#   X.col <- rbind(X.col, raw.data[[j]])
# X.col <- X.col[-1,]
# 
# G.col <- t(X.col) %*% X.col
# ev.col <- eigen(G.col)
# V <- Re(ev.col$vectors[,1:60])
# 
# X.row <- rep(0, 112)
# for (j in 1:400)
#   X.row <- cbind(X.row, raw.data[[j]])
# X.row <- X.row[,-1]
# 
# G.row <- X.row %*% t(X.row)
# ev.row <- eigen(G.row)
# U <- Re(t(ev.row$vectors[,1:60]))
# 
# 
# s <- 0
# for (i in 1:400)
#   s <- s + norm(raw.data[[i]] - t(U) %*% U %*% raw.data[[i]] %*% V %*% t(V), type = "F")
