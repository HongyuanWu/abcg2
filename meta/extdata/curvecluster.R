

#function linear model
library(fda)
library(MASS)

fpca.test.gwas.group<- function(snp.name, gene.group, snp.position=NULL, cc, data,basis=37,rng=50000){
  data = as.matrix(data)
  if( !is.matrix(data) ) stop('Data cannot be convert to matrix')
  if( length(cc) != dim(data)[1] ) stop('The dimension of Affected Status is not match with sample size')
  rlt<-list()      				#initial the result
  rlt$gene.symbol <- vector()			#record the gene symbol of each test
  rlt$p.value.80<-vector()			#get result with 80% of the principal component variance
  rlt$p.value.90<-vector()			#get result with 90% of the principal component variance
  rlt$test.snp.number<-vector()		#get the tested SNP numbers in each gene region
  
  gene.unique<-unique(gene.group)
  for ( i in 1:length(gene.unique)  ) {  		#Test each gene in gene.csv
    sub.idx<- gene.group==gene.unique[i] 	
    x<-fdata[,snp.name[sub.idx]]			#get the test genotype data
    rlt$gene.symbol[i] <- gene.symbol[i]	#record the tested gene symbol 
    pos<- snp.position[sub.idx]			#get the test SNP position
    if ( length(pos) >= 3 ){
      fpca.rlt<-fpca.genotype(x,cc,pos, nbasis=basis)	#test the p-value
      rlt$test.snp.number[i]<-length(pos)			#record how many snps is in the test
      rlt$p.value.80[i]<-fpca.rlt$pv.all[fpca.rlt$prop>0.8][1]
      rlt$p.value.90[i]<-fpca.rlt$pv.all[fpca.rlt$prop>0.9][1]
    }else{
      rlt$test.snp.number[i]<-0			
      rlt$p.value.80[i]<- -1
      rlt$p.value.90[i]<- -1
    }
  }  
  data.frame(Gene_Symbol=rlt$gene.symbol,SNP_tested=rlt$test.snp.number,FPCA_P_Value_80=rlt$p.value.80,
             FPCA_P_Value_90=rlt$p.value.90,	  Chi_Square_Permutation=rlt$p.value.chi.square, CMC=rlt$p.value.cmc,
             Collaps=rlt$p.value.collaps, VT=rlt$p.value.vt, WE=rlt$p.value.we, T_Square=rlt$p.value.t2)
}

fpca.test.gwas<- function( gene.symbol, gene.chrom,gene.start,gene.end, snp.position,snp.chrom, cc, geno,basis=37,rng=50000){
  geno = as.matrix(geno)
  if( !is.matrix(geno) ) stop('Data cannot be convert to matrix')
  if( length(cc) != dim(geno)[1] ) stop('The dimension of Affected Status is not match with sample size')
  rlt<-data.frame( )						#initial the result
  for ( i in 1:dim(gene.list)[1]  ) {  		#Test each gene in gene.csv
    sub.idx<- snp.chrom == gene.chrom[i]   	#specify the chromosome
    sub.idx<- sub.idx & (snp.position > (gene.start[i] -rng) )  #choose the SNPs from the gene start minus range set before
    sub.idx<- sub.idx & (snp.position < (gene.end[i] + rng) ) #choose the SNPs from the gene end plus range set before
    x<-geno[,sub.idx]			#get the test genotype data
    rlt[i,"Gene_symbol"] <- gene.symbol[i]	#record the tested gene symbol 
    pos<- snp.position[sub.idx]		#get the test SNP position
    x<-as.matrix(x)
    
    if( dim(x)[2] >= 3){
      maf<-colMeans(x,na.rm=TRUE)
      x[,maf>1]=2-x[,maf>1]
      x[is.na(x)]=0
      maf=colMeans(x)
      pos<-pos[maf>0]
      x=x[,maf>0]
    }
    if ( length(pos) >= 3 ){      # no need since exon usually greater than 200bp
      pos2= 0:(length(pos)-1)
      pos2=(pos2-pos2[1])/(pos2[length(pos2)]-pos2[1])
      fpca.rlt<- try( fpca.genotype(x,cc,pos, nbasis=basis)	)#test the p-value
      
      rlt[i,"Snp_Number"]<-length(pos)			#record how many snps is in the test
      rlt[i,"FPCA_80"]<-try( fpca.rlt$pv.all[fpca.rlt$prop>0.8][1] )
      rlt[i,"FPCA_90"]<-try( fpca.rlt$pv.all[fpca.rlt$prop>0.9][1] )
      
    }
    print( i/dim(gene.list)[1])
  }  
  rlt
}

fourier.expansion<- function(x,n_of_basis,pos){
  frange <- c(pos[1], pos[length(pos)])
  rlt=list();
  rlt$fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
  N=dim(x)[2]
  rlt$phi = eval.basis(pos,rlt$fbasis);
  rlt$coef<-ginv(t(rlt$phi)%*%rlt$phi)%*%t(rlt$phi)%*%t(x)
  rlt$x_hat <- t( rlt$phi  %*% rlt$coef );
  return(rlt)	
}

fourier.expansion.smoothed<- function(x,n_of_basis,pos,lambda){
  frange <- c(pos[1], pos[length(pos)])
  rlt=list();
  rlt$fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
  rlt$penalty<-fdPar(rlt$fbasis,2,lambda)
  N=dim(x)[2]	
  
  rlt$phi = eval.basis(pos,rlt$fPar$basis);
  rlt$coef<-ginv(t(rlt$phi)%*%rlt$phi)%*%t(rlt$phi)%*%t(x)
  
  rlt$x_hat <-t( rlt$phi%*%rlt$coefs );
  return(rlt)	
}

fpca.genotype <- function(x,pos=NULL,percentage=0.8,nbasis=37,lambda=NULL){
  nsnps <- dim(x)[2]
  ninds <- dim(x)[1]
  if ( is.null(pos) ){
    pos <- (0:( nsnps-1) )/(nsnps-1)
  }else {
    idx<-order(pos)
    x<-x[,idx]
    pos<-pos[idx]
    pos<- (pos-pos[1])/(pos[nsnps]-pos[1])
  }
  if( is.null(lambda)){
    expanded<-fourier.expansion(x,nbasis,pos)
  }else if( lambda > 0 ){
    expanded<-fourier.expansion.smoothed(x,nbasis,pos,lambda)
  }else if( lambda == 0 ){
  }
  coef<-t(expanded$coef-rowMeans(expanded$coef))/sqrt(ninds) 
  pca.rlt<-prcomp(coef)
  pca.rlt$scores<-coef%*%pca.rlt$rotation
  v0<-diag(var(pca.rlt$scores))
  rlt<-list()
  rlt$prop<-cumsum(v0)/sum(v0)
  rlt$pv<-rlt$pv.all[rlt$prop>percentage][1]
  rlt
}

#simulation for cluster analysis
mydata=matrix(rnorm(3000),60,50)
mydata[12:26,12:26]=rnorm(15*15,3,0.3)
mydata[27:44,27:44]=rnorm(18*18,10,0.3)
mydata[45:50,45:50]=rnorm(6*6,30,0.3)

mydata <- na.omit(mydata) # listwise deletion of missing
mydata <- scale(mydata) # standardize variables 



###One way cluster
# Ward Hierarchical Clustering(uni-cluster for row)  Jan 31st 2013
fhclust<-function(mydata,method="euclidean",cut,bootstrap){
  d <- dist(mydata, method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward")
  plot(fit) # display dendogram
  clusters <- cutree(fit, k=cut) # cut tree into 5 clusters
  # draw dendogram with red borders around the 5 cluster
  rect.hclust(fit, k=cut, border="red") 
  if(bootstrap==TRUE){
    library(pvclust)
    # Ward Hierarchical Clustering with Bootstrapped p values
    fit <- pvclust(mydata, method.hclust="ward",method.dist="euclidean")
    plot(fit) # dendogram with p values
    # add rectangles around groups highly supported by the data
    pvrect(fit, alpha=.95)   
  }else{}
  rlt<-list()
  rlt$group<-clusters;
}

# K-Means Clustering with 5 clusters  Jan 31st 2013
fkmeans<-function(mydata,cut){
  library(cluster)
  fit <- kmeans(mydata, cut)
  # Cluster Plot against 1st 2 principal components
  # vary parameters for most readable graph
  clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE,labels=2, lines=0)
  rlt<-list()
  rlt$group<-fit$cluster
}

# Model Based Clustering   Jan 31st 2013
fmclust<-function(mydata){
  library(mclust)
  fit <- Mclust(mydata)
  plot(fit, "BIC") # plot results
  print(fit) # display the best model 
  rlt<-list()
  rlt$group<-fit$classification
}

# select feature with laso and then do hierarchical cluster to treat spare problem  #Feb 1 2013
sparclfeature<-function(mydata,method="complete",cut=4){
  library("sparcl")
  set.seed(1)
  # Do tuning parameter selection for sparse hierarchical clustering
  perm.out <- HierarchicalSparseCluster.permute(x, wbounds=c(1.5,2:9),nperms=5)
  # Perform sparse hierarchical clustering
  sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists,wbound=perm.out$bestw,method="complete")
  m<-x[,which(sparsehc$ws !=0)]
  #drawHeatmap2(m)
  d <- dist(m, method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward")
  #plot(fit) # display dendogram
  clusters <- cutree(fit, k=cut) # cut tree into 5 clusters
  # draw dendogram with red borders around the 5 cluster
  rect.hclust(fit, k=cut, border="red") 
  #retrun output
  rlt<-list()
  rlt$featurenumber<-sum(sparsehc$ws !=0)
  rlt$group<-clusters;
} 

###Two way cluster(bicluster)
# biclust Clustering
fbiclust<-function(mydata){
  library("biclust")
  set.seed(1)
  bics <- biclust(mydata,BCPlaid(), back.fit = 2, shuffle = 3, fit.model = ~m + a + b,iter.startup = 5, iter.layer = 30, verbose = TRUE)
  rlt<-list()
  r1<-fbicorder(bics, cols=FALSE, rev=FALSE)
  r2<-fbicorder(bics, cols=TRUE, rev=FALSE)
  rlt$row=r1
  rlt$col=r2
  rlt
}

#isa2 cluster analysis 
fisa2<-function(mydata){
  library("isa2")
  isa.result <- isa(mydata)
  
  bc <- isa.biclust(isa.result)
  # drawHeatmap(mydata, bc, 1)
  #bubbleplot(mydata, bc)
  parallelCoordinates(mydata, bc, number=1)
  
  names(isa.result)
  ## Find the best bicluster for each block in the input
  best <- apply(cor(isa.result$rows, data[[2]]), 2, which.max)
  ## Check correlation
  sapply(seq_along(best),function(x) cor(isa.result$rows[,best[x]], data[[2]][,x]))
  ## The same for the columns
  sapply(seq_along(best),function(x) cor(isa.result$columns[,best[x]], data[[3]][,x]))
  ## Plot the data and the modules found
  if (interactive()) {
    layout(rbind(1:2,3:4))
    image(data[[1]], main="In-silico data")
    sapply(best, function(b) image(outer(isa.result$rows[,b],
                                         isa.result$columns[,b]),
                                   main=paste("Module", b)))
  }
}


###Cluster compare
# Centroid Plot against 1st 2 discriminant functions
clusterplot<-function(mydata,group1){
  library(fpc)
  # plotcluster(mydata, fit$cluster) 
}

# compare different cluster reslut
clustercompare<-function(mydata,group1,group2){
  library(fpc)
  d<-dist(mydata)
  cluster.stats(d,group1 ,group2)  
}

#### Function to order variables or objects that appear in a bicluster[function bicorder]
fbicorder<-function(bicResult, cols=TRUE, rev=FALSE)
{
  i<-numeric()
  res<-c()
  order<-vector();
  if(!cols){
    le<-dim(bicResult@RowxNumber)[2]
    for(i in 1:le){
      order[which(bicResult@RowxNumber[,i])[!(which(bicResult@RowxNumber[,i]) %in% res)]]<-i;
      res<-c(res,which(bicResult@RowxNumber[,i])[!(which(bicResult@RowxNumber[,i]) %in% res)])
    }
    count<-1:dim(bicResult@RowxNumber)[1]
    order[count[!(count %in% res)]]<-0;
  }else{
    le<-dim(bicResult@NumberxCol)[1]
    for(i in 1:le){
      order[which(bicResult@NumberxCol[i,])[!(which(bicResult@NumberxCol[i,]) %in% res)]]<-i;
      res<-c(res,which(bicResult@NumberxCol[i,])[!(which(bicResult@NumberxCol[i,]) %in% res)])
    }
    count<-1:dim(bicResult@NumberxCol)[2]
    order[count[!(count %in% res)]]<-0;
  }
  
  if(rev) order<-rev(order)
  order
}




group1<-fhclust(mydata,cut=4,bootstrap=F)
group2<-fkmeans(mydata,cut=4)
clusterplot(mydata,group1)
clustercompare(mydata,group1,group2)

group1
group2
length(group1)
dim(mydata)
