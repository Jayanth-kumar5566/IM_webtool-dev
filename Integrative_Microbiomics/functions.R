#Functions to source into the app
library(SNFtool)
library("vegan")
library("reticulate")
use_python("/usr/local/bin/python")
source_python("Python_codes/plot_mat.py")
source_python("Python_codes/sil.py")
#=================Trail data====================
#data<-read.csv("./../Data/bacteria.csv",header = TRUE,row.names = 1)
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
      x[[i]]=(as.matrix(dsim)-1)*-1
      }
  }
  x<-x[-which(sapply(x, is.null))]
  W = SNF(x,k,t)
  #print(W)
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
