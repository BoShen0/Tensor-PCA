######################################################################
######################################################################
#####
##### This code is the Implementation of TPCA-L1 from the below paper
##### Pang, Yanwei, Xuelong Li, and Yuan Yuan. "Robust tensor analysis with L1-norm."
##### IEEE Transactions on Circuits and Systems for Video Technology 20.2 (2010): 172-178.
#####
######################################################################
######################################################################
############ polarity.function
polarity.f <- function(xxxx)
{
  if (xxxx>0) return(1)
  else return(-1)
}

########## to get one u from the data  
########## u is form the left side
TPCAL1.u <- function(DATA, u, epsilon.u, max.u){
  N <- length(DATA)  ### sample size
  m <- dim(DATA[[1]])[1]  ## number of rows
  n <- dim(DATA[[1]])[2]  ## number of columns
  u.past <- -u  ####### 
  t <- 1
  while (round(norm(u - u.past,"2"),20 )> epsilon.u ){
    u.past <- u   ### last iteration vector
    ####### to get a
    a <- rep(0,m)
    #################################################################
    for(jj in 1:N){ 
      Y <- DATA[[jj]]    ####### jj_th sample
      for (ii in 1:n){
        a <- a + polarity.f(sum(u.past*Y[,ii])) * Y[,ii]
      }
    }
    #################################################################
    u <- a/norm(a, "2")
    t <- t + 1
    if(t > max.u) break
  }
  result <- list(u = u, t = t)
  return(result)
} 

########## to get one v from the data  
########## v is form the right side
TPCAL1.v <- function(DATA, v, epsilon.v, max.v){
  N <- length(DATA)  ### sample size
  m <- dim(DATA[[1]])[1]  ## number of rows
  n <- dim(DATA[[1]])[2]  ## number of columns
  v.past <- -v 
  t <- 1
  while (round(norm(v - v.past,"2"),20 ) > epsilon.v ){
    v.past <- v 
    ####### to get a
    a <- rep(0, n)
    #################################################################
    for(jjj in 1:N){
      Y <- DATA[[jjj]]
      for (iii in 1:m){
        a <- a + polarity.f(sum(v.past * Y[iii,])) * Y[iii,]
      }
    }
    #################################################################
    v <- a/norm(a, "2")
    t <- t + 1
    if(t > max.v) break
  }
  result <- list(v = v, t = t)
  return(result)
} 


######### run all above functions first and related packages
######### the main function for TPCAL1
TPCA.L1 <- function(Data, p, q, epsilon.inner, maxit.inner, epsilon.UV, maxiteration){
  ############ Data should be list form, for each element in the list is a matrix
  ############ p, q is the number of feature for U \in R^{p*m} and V \in R^{n*q}
  ############ epsilon.inner and maxit.inner are for inside iterations
  ############ epsilon.UV and maxiteration are for outside iterations 
  timestart<-Sys.time()
  N <- length(Data)  ### sample size
  m <- dim(Data[[1]])[1]  ## number of rows
  n <- dim(Data[[1]])[2]  ## number of columns
  ######### initialize U,V
  ######### initialization with slicing learning i.e. HOSVD
  X.col <- rep(0, n)
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
  U <- Re(t(ev.row$vectors[,1:p]))
  ######### this is to record U,V
  #################################################################
  U.past <- matrix(0, p, m)
  V.past <- matrix(0, n, q)
  #################################################################
  t=0
  UVGap.track <- NULL
  
  #################################################################
  while (round(norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
               norm(V %*% t(V) - V.past %*% t(V.past), type = "F"),20 )> epsilon.UV ) {
    
    print(paste("number of iteration:",t + 1))
    U.past <- U
    V.past <- V
    ########## fix V, try get a new U
    #################################################################
    Data.V <- NULL
    for(i in 1:N)
      Data.V[[i]] <- Data[[i]] %*% V
    
    for (i in 1:p)
    {
      FF.u <-  TPCAL1.u(Data.V, U[i,], epsilon.inner, maxit.inner) 
      u.new <- FF.u$u    #### to get updated u
      U[i,] <- u.new
      u.matrix <- matrix(u.new, m, 1)
      for(i in 1:N)
        Data.V[[i]] <- Data.V[[i]] - u.matrix %*% t(u.matrix) %*% Data.V[[i]]
    }
    
    Data.U <- NULL
    for(i in 1:N)
      Data.U[[i]] <- U %*% Data[[i]]  
    
    for (i in 1:q)
    {
      FF.v <-  TPCAL1.v(Data.U, V[,i], epsilon.inner, maxit.inner) 
      v.new <- FF.v$v
      V[,i] <- v.new
      v.matrix <- matrix(v.new,n,1)
      for(i in 1:N)
        Data.U[[i]] <- Data.U[[i]] - Data.U[[i]] %*% v.matrix %*% t(v.matrix)  
    }
    
    #################################################################
    t <- t + 1
    UVGap.track[t] <- norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
      norm(V %*% t(V) - V.past %*% t(V.past), type = "F")
    print(paste("U,V gap",UVGap.track[t]))
    cat("\n")
    if(t > maxiteration) break
  }
  
  timeend<-Sys.time()
  runningtime<-timeend-timestart
  
  result <- list(U = U, V = V, t = t - 1, UV.Gap = UVGap.track,Time = runningtime)
  return(result)
}
mm <- TPCA.L1(raw.data[1:20],80, 80, 10^-3, 20, 10^-3, 100)

U <- mm$U
V <- mm$V
image(t(U) %*% U %*% raw.data[[5]] %*% V %*% t(V),col=grey.colors(255))
image(raw.data[[5]] ,col=grey.colors(255))
U %*% t(U) 
dim(V)
 