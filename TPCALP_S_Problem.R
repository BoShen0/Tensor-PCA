######################################################################
######################################################################
################ Implementation of TPCALp-S   ########################
######################################################################
######################################################################
 TPCALP.u <- function(Data, u, lambda, gamma, epsilon, PP){
   N <- length(Data)  ### sample size
   m <- dim(Data[[1]])[1]  ## number of rows
   n <- dim(Data[[1]])[2]
   u.past <- -u 
   while (round(norm(u-u.past,"2"),20 )> epsilon ){
     u.past <- u 
     ####### to get a
     a <- rep(0,m)
     #################################################################
     for(jj in 1:N){
       Y <- Data[[jj]]
       for (ii in 1:n){
         if (sum(u*Y[,ii]) == 0) {
           next
         }
         a <- a + sign(sum(u*Y[,ii])) * abs(sum(u*Y[,ii]))^(PP-1) * Y[,ii]
       }
     }
     #################################################################
     b <- abs(u)/(lambda + gamma * abs(u))
     u <- a * b
     u <- u/norm(u,"2")
   }
   result <- list(u = u)
   return(result)
 }   
 
 
 
TPCALP.v <- function(Data, v,lambda, gamma, epsilon, PP){
   N <- length(Data)  ### sample size
   m <- dim(Data[[1]])[1]  ## number of rows
   n <- dim(Data[[1]])[2]
   v.past <- -v 
   while (round(norm(v-v.past,"2"),20 )> epsilon ){
     v.past <- v 
     ####### to get a
     a <- rep(0,n)
     #################################################################
     for(jjj in 1:N){
       Y <- Data[[jjj]]
       for (iii in 1:m){
         if (sum(v*Y[iii,]) == 0) {
           next
         }
         a <- a + sign(sum(v*Y[iii,])) * abs(sum(v*Y[iii,]))^(PP-1) * Y[iii,]
       }
     }
     #################################################################
     b <- abs(v)/(lambda + gamma * abs(v))
     v <- a * b
     v <- v/norm(v,"2")
   }
   result <- list(v = v)
   return(result)
 }   
 
 
 
 
FF.u <-  TPCALP.u(raw.data, rep(sqrt(1/112),112),1,10^-1,10^-3,0.5) 
FF.v <-  TPCALP.v(raw.data, rep(sqrt(1/92),92),1,10^-1,10^-3,0.5) 
FF.u$u
FF.v$v
 
TPCALP.S <- function(Data, p, q, lambda, gamma, epsilon.uv,epsilon.UV, maxiteration,PPP){
   timestart<-Sys.time()
   N <- length(Data)  ### sample size
   m <- dim(Data[[1]])[1]  ## number of rows
   n <- dim(Data[[1]])[2]  ## number of columns
   ######initialize U,V
   #########initialization with slicing learning i.e. HOSVD
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
   
   
   U.past <- matrix(0, p, m)
   V.past <- matrix(0, n, q)
   for(i in 1:p)
     U.past[i,i] <- 1
   
   for(i in 1:q)
     V.past[i,i] <- 1 
   
   t=0
   UVGap.track <- NULL
   
  
   
   while (round(norm(t(U) %*% U - t(U.past) %*% U.past, type = "F") + 
                norm(V %*% t(V) - V.past %*% t(V.past), type = "F"),20 )> epsilon.UV ) {
     print(paste("number of iteration:",t + 1))
     U.past <- U
     V.past <- V
     ##########to get new U
     Data.V <- NULL
     for(i in 1:N)
       Data.V[[i]] <- Data[[i]] %*% V
     
     for (i in 1:p)
     {
       FF.u <-  TPCALP.u(Data.V, U[i,], lambda, gamma, epsilon.uv, PPP) 
       u.new <- FF.u$u
       U[i,] <- u.new
       u.matrix <- matrix(u.new,m,1)
       for(i in 1:N)
         Data.V[[i]] <- Data.V[[i]] -u.matrix %*% t(u.matrix) %*% Data.V[[i]]
     }
     
     Data.U <- NULL
     for(i in 1:N)
       Data.U[[i]] <- U %*% Data[[i]]  
     
     for (i in 1:q)
     {
       FF.v <-  TPCALP.v(Data.U, V[,i], lambda, gamma, epsilon.uv, PPP) 
       v.new <- FF.v$v
       V[,i] <- v.new
       v.matrix <- matrix(v.new,n,1)
       for(i in 1:N)
         Data.U[[i]] <- Data.U[[i]] - Data.U[[i]] %*% v.matrix %*% t(v.matrix)  
     }
     
   
     
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
ABC <- TPCALP.S(raw.data, 25, 25, 10^-2,10^-2, 10^-3,10^-3, 200, 1)
U <- ABC$U
V <- ABC$V
ABC$UV.Gap
image(t(U) %*% U %*% raw.data[[1]] %*% V %*% t(V),col=grey.colors(255))
image(raw.data[[1]],col=grey.colors(255))
t(U) %*% U
