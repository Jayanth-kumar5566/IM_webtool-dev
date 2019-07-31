#Functions to source into the app
library(SNFtool)
library("vegan")
library("reticulate")
use_python("/usr/local/bin/python")
source_python("Python_codes/plot_mat.py")
source_python("Python_codes/sil.py")
#=================Trail data====================
#data1<-read.csv("./../Data/bacteria.csv",header = TRUE,row.names = 1)
#data2<-read.csv("./../Data/fungi.csv",header = TRUE,row.names = 1)
#data3<-read.csv("./../Data/virus.csv",header = TRUE,row.names = 1)
#=====Function for plotting individual biomes===================
biome_plot<-function(data,k){
  dsim=vegdist(data,method='bray',diag=TRUE,upper=TRUE)
  sim=(as.matrix(dsim)-1)*-1
  sim[is.nan(sim)]=1
  labels=spectralClustering(sim,k)
  labels=as.data.frame(labels,row.names = row.names(data))
  sim=as.data.frame(sim)
  #bplot(sim,labels)
  m<-bplot(sim,labels)
  return(heatmap(as.matrix(m),diss,Rowv = NA, Colv = NA, scale = "none",
                 main = "Spectral clustering",xlab = "Sample ID",ylab = "Sample ID",labRow = "",labCol = ""))
}

#=====Function for plotting log of individual biomes===================
biome_log_plot<-function(data,k){
  dsim=vegdist(data,method='bray',diag=TRUE,upper=TRUE)
  sim=(as.matrix(dsim)-1)*-1
  sim[is.nan(sim)]=1
  labels=spectralClustering(sim,k)
  labels=as.data.frame(labels,row.names = row.names(data))
  sim=as.data.frame(sim)
  #bplot(sim,labels)
  sim<-log(sim)
  m<-bplot(sim,labels)
  return(heatmap(as.matrix(m),diss,Rowv = NA, Colv = NA, scale = "none",
                 main = "Spectral clustering",xlab = "Sample ID",ylab = "Sample ID",labRow = "",labCol = ""))
}

#===========Function for giving k based on maximum silhoutee score============
#Given Similarity matrix and labels
max_k<-function(sim){
  k_sil=list()
  for (i in 2:10){
    incProgress(1/10, detail = paste("Calculating", i))
    labels=spectralClustering(sim,i)
    labels=as.data.frame(labels,row.names = row.names(sim))
    k_sil[[i]]<-silhouette_score(sim,labels)}
    k_sil=data.frame(number_of_clusters=2:10,unlist(k_sil))
    colnames(k_sil)<-c("Number of Clusters","Silhouette Score")
    return(k_sil)
  }
#=================Merging Function==================================
merge_snf<-function(x,k,t){
  #x<-list(data1,data2,data3)
  for (i in 1:3){
    incProgress(1/3, detail = paste("Merging", i))
    if (is.null(x[[i]])==FALSE){
      dsim=vegdist(x[[i]],method='bray',diag=TRUE,upper=TRUE)
      dsim[is.nan(dsim)]<-0 
      x[[i]]=(as.matrix(dsim)-1)*-1
      }
  }
  if (sum(sapply(x,is.null))>0){
  x<-x[-which(sapply(x, is.null))]
  }
  W = SNF(x,k,t)
  #print(W)
  return(W)
}
#===============Weighted SNF Merging Function===================
SNF_weighted_iter<-function (Wall, K = 20, t = 20, weight) 
{
  check_wall_names <- function(Wall) {
    name_match <- function(names_A, names_B) {
      return(identical(dimnames(names_A), dimnames(names_B)))
    }
    return(all(unlist(lapply(Wall, FUN = name_match, Wall[[1]]))))
  }
  wall.name.check <- check_wall_names(Wall)
  wall.names <- dimnames(Wall[[1]])
  if (!wall.name.check) {
    warning("Dim names not consistent across all matrices in Wall.\n            Returned matrix will have no dim names.")
  }
  LW <- length(Wall) #Total number of views
  weight=weight/sum(weight)
  weight=weight*(LW)
  normalize <- function(X) {
    row.sum.mdiag <- rowSums(X) - diag(X)
    row.sum.mdiag[row.sum.mdiag == 0] <- 1 #making zeros 1
    X <- X/(2 * (row.sum.mdiag))
    diag(X) <- 0.5
    return(X)
  }
  newW <- vector("list", LW) #Stores S matrices
  nextW <- vector("list", LW) #Stores fused matrices
  #Creates Q matrix
  for (i in 1:LW) {
    Wall[[i]] <- normalize(Wall[[i]])
    Wall[[i]] <- (Wall[[i]] + t(Wall[[i]]))/2
  }
  #Creates S matrix
  for (i in 1:LW) {
    newW[[i]] <- (SNFtool:::.dominateset(Wall[[i]], K))
  }
  #Fusion
  for (i in 1:t) {
    for (j in 1:LW) {
      sumWJ <- matrix(0, dim(Wall[[j]])[1], dim(Wall[[j]])[2])
      tmp<-0
      for (k in 1:LW) {
        if (k != j) {
          sumWJ <- sumWJ + weight[k]*Wall[[k]]  #possible weighting option 1
          tmp<-tmp+weight[k]
        }
      }
      nextW[[j]] <- newW[[j]] %*% (sumWJ/(tmp)) %*% 
        t(newW[[j]])
    }
    for (j in 1:LW) {
      Wall[[j]] <- normalize(nextW[[j]])
      Wall[[j]] <- (Wall[[j]] + t(Wall[[j]]))/2
    }
  }
  W <- matrix(0, nrow(Wall[[1]]), ncol(Wall[[1]]))
  for (i in 1:LW) {
    W <- W + Wall[[i]]  
  }
  W <- W/LW
  W <- normalize(W)
  W <- (W + t(W))/2
  if (wall.name.check) {
    dimnames(W) <- wall.names
  }
  return(W)
}

merge_wsnf<-function(x,k,t,weight){
  #x<-list(data1,data2,data3)
  #weight<-c(weight1,weight2,weight3)
  for (i in 1:3){
    incProgress(1/3, detail = paste("Merging", i))
    if (is.null(x[[i]])==FALSE){
      dsim=vegdist(x[[i]],method='bray',diag=TRUE,upper=TRUE)
      dsim[is.nan(dsim)]<-0  #Dissimilarity is 0 if both patients have no microbes
      x[[i]]=(as.matrix(dsim)-1)*-1
    }
  }
  if (sum(sapply(x,is.null))>0){
    x<-x[-which(sapply(x, is.null))]
  } 
  W = SNF_weighted_iter(x,k,t,weight)
  return(W)
}
#=============Label Creating Function=======================
label_create<-function(sim,k){
  labels=spectralClustering(sim,k)
  labels=as.data.frame(labels,row.names = row.names(sim))
}
#============Function for barplot=====================
bar<-function(data,k){
  summ<-colSums(data)
  summ<-sort(summ,decreasing = TRUE)
  par(mar=c(1,10,1,1)+.1)
  barplot(sort(summ[1:k]),las=2,axes=FALSE, col="darkblue",horiz = TRUE)
}
