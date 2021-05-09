######################################################################
######################################################################
#####
##### This code is the Implementation of 2d^2 PCA from the below paper
##### Zhang, D., & Zhou, Z. H. (2005). (2D) 2PCA: Two-directional two-dimensional PCA 
##### for efficient face representation and recognition. Neurocomputing, 69(1-3), 224-231.
#####
######################################################################
######################################################################
TwodtwoPCA <- function(Data, p, q){
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

result <- list(U = U, V = V)
return(result)

}

mm <- TwodtwoPCA(raw.data[1:20],20,20)

U <- mm$U
V <- mm$V
image(t(U) %*% U %*% raw.data[[5]] %*% V %*% t(V),col=grey.colors(255))
image(raw.data[[5]] ,col=grey.colors(255))

