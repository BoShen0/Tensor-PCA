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

path <- "C:/Users/boshen.HOKIES/Desktop/take home/Research/CCA_coding new/Case study/att_faces"
fileNames <- dir(path)  
sub.path <- "C:/Users/boshen.HOKIES/Desktop/take home/Research/CCA_coding new/Case study/att_faces/s2"
sub.fileNames <- dir(sub.path)
filePath <- data.frame(matrix(0, length(fileNames), 10))
for(i in 1:40)
  for (j in 1:10)
    filePath[i,j] <- paste(path,fileNames[i],sub.fileNames[j],sep='/')

label <- rep(1:40, rep(10, 40))
raw.data <- NULL
for (i in 1:40)
  for(j in 1:10){
    Y <- read.pnm(filePath[i, j])
    raw.data[[10 * (i - 1) + j]] <- as.matrix(Y@grey) 
  }

twodpca <- function(X, d){ #####dims of output data
  N <- length(X)
  m <- dim(X[[1]])[1]
  n <- dim(X[[1]])[2]
  E.mean <- matrix(0,m,n)
  for(i in 1:N){
    E.mean <- E.mean + X[[i]]
  }
  E.mean <- E.mean/N
  
  E.matrix = matrix(0, m, m)
  for(i in 1:N){
    E.matrix <- E.matrix + (X[[i]] - E.mean) %*% t(X[[i]] - E.mean)
  }
  eig.problem <- eigen(E.matrix)
  data = NULL
  for(i in 1:N){
    data[[i]] <- as.matrix(eig.problem$vectors)[1:d,] %*% X[[i]]
  }
  
  result <- list(eigenValue=eig.problem$values, eignVector=as.matrix(eig.problem$vectors), data = data )
  return(result)
  }
  


twodpca.t <- function(X, d){
  N <- length(X)
  m <- dim(X[[1]])[1]
  n <- dim(X[[1]])[2]
  E.mean <- matrix(0,m,n)
  for(i in 1:N){
    E.mean <- E.mean + X[[i]]
  }
  E.mean <- E.mean/N
  
  E.matrix = matrix(0, n, n)
  for(i in 1:N){
    E.matrix <- E.matrix + t(X[[i]] - E.mean) %*% (X[[i]] - E.mean)
  }
  eig.problem <- eigen(E.matrix)
  
  data = NULL
  for(i in 1:N){
    data[[i]] <- X[[i]] %*% as.matrix(eig.problem$vectors)[,1:d]
  }
  
  result <- list(eigenValue=eig.problem$values, eignVector=as.matrix(eig.problem$vectors), data = data )
  return(result)
}

# eigen.val1 <- matrix(0,100,112)
# eigen.val2 <- matrix(0,100,92)
# 
# for (j in 1:100){
#   result1 <- twodpca(raw.data)
#   eigen.val1[j,] <- result1$eigenValue
#   raw.data <- lapply(raw.data, function(x) result1$eignVector %*% x)
#   
#   result2 <- twodpca.t(raw.data)
#   eigen.val2[j,] <- result2$eigenValue
#   raw.data <- lapply(raw.data, function(x) x %*% result2$eignVector)
# }
# length( eigen.val1[j,])
# dim( eigen.val1)
# 
# lapply(list(1,-1,2,3), function(x) x^2)

re1 <- twodpca(raw.data,27)
re2 <- twodpca.t(raw.data,26)

pca.data <- matrix(0,400, 27*92 + 112*26)
for (i in 1:400){
  pca.data[i,] <- as.vector(c(re1$data[[i]],re2$data[[i]]))
} 


set.seed(100)
figure.knn <- knn(pca.data[train.index[5,],], pca.data[-train.index[5,],], 
                  label.training,k=1,prob=T,use.all=T)  ###k=

1 - classError(figure.knn,label.test )$errorRate





pca.data <- matrix(0,400, 26 * 27)
for (i in 1:400){
  pca.data[i,] <- as.vector(re1$eignVector[1:27,] %*% raw.data[[i]] %*% re2$eignVector[,1:26])
} 

#############################
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

set.seed(100)
figure.knn <- knn(pca.data[train.index[1,],], pca.data[-train.index[1,],], 
                  label.training,k=1,prob=T,use.all=T)  ###k=

1 - classError(figure.knn,label.test )$errorRate



data1 <- NULL
for (i in 1:400){
  data1[[i]] <- re1$eignVector[1:27,] %*% raw.data[[i]]
} 
re3 <- twodpca.t(data1)
re3$eignVector

pca2.data <- matrix(0,400, 26 * 27)
for (i in 1:400){
  pca2.data[i,] <- as.vector(data1[[i]] %*% re3$eignVector[,1:26])
} 

set.seed(100)
figure.knn <- knn(pca2.data[train.index[1,],], pca2.data[-train.index[1,],], 
                  label.training,k=1,prob=T,use.all=T)  ###k=

1 - classError(figure.knn,label.test )$errorRate



# ##########NBC
# set.seed(100)
# NBC.Lscca <- naiveBayes(pca2.data[train.index[1,],],as.factor(label.training))
# predict.label <- as.numeric(predict(NBC.Lscca, pca2.data[-train.index[1,],]))
# 
# 1 - classError(predict.label,label.test )$errorRate
# 
# 
# ##########SVM
# set.seed(100)
# svp.Lscca <- ksvm(pca2.data[train.index[1,],],as.factor(label.training),type="C-svc",kernel='rbfdot',C=1,scaled=c())
# predict.label = as.numeric(predict(svp.Lscca,pca2.data[-train.index[1,],]))
# 1 - classError(predict.label,label.test)$errorRate
