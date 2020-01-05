#Packages
#install.packages("devtools")
#install_github("zdk123/SpiecEasi", force=TRUE)
#source("https://bioconductor.org/biocLite.R")

library(SpiecEasi)
library(igraph)

#find x position in ref
find_pos <- function(x,ref)
{
  res <- as.vector(NULL)
  for(i in 1:length(x))
  {
    res<-c(res, grep(x[i],ref))
  }
  return(res)
}
#Sample Division
sample_division <- function(samples,table)
{
  Sample_list <- list(NULL)
  length(Sample_list) <- length(Sample.Division)
  for(i in 1:length(samples))
  {
    sample.names <- colnames(table)
    temp <- find_pos(samples[i],sample.names)
    Sample_list[[i]] <- mOTU.table[,temp]
  }
  names(Sample_list) <- samples
  return(Sample_list)
}
#filtering otus
filtering_otu <- function(otutable,otu_obs=0.5,Transpose=FALSE)
{
  if(Transpose) {otutable <- t(otutable)}
  if(otu_obs >0)
  {
    sample_num <- nrow(otutable)
    sample_obs_th <- otu_obs*sample_num
    otu_absence <- apply(otutable,2,function(x){length(which(x!=0))})
    otu_keep <- which(otu_absence >= sample_obs_th)
    otutable <- otutable[,otu_keep]
  }
  return(otutable)
}
#Network construction
spiec_easi_network <- function(otutable,Transpose=FALSE,mb_th=0.3,lambda.min.ratio=1e-2,nlambda=20,rep.num=20)
{
  otutable <- as.matrix(otutable)
  if(Transpose) {otutable <- t(otutable)}
  se.mb.amgut <- spiec.easi(otutable,method="mb",lambda.min.ratio=lambda.min.ratio,nlambda=nlambda
                            ,pulsar.params=list(rep.num=rep.num))
  adj_mb=as.matrix(symBeta(getOptBeta(se.mb.amgut), mode='maxabs'))
  colnames(adj_mb) <- colnames(otutable)
  rownames(adj_mb) <- colnames(otutable)
  adj_mb <- abs(adj_mb) > mb_th
  diag(adj_mb) <- 0
  adj_mb_v <- adj_mb*as.matrix(symBeta(getOptBeta(se.mb.amgut), mode='maxabs'))
  colnames(adj_mb_v) <- colnames(otutable)
  rownames(adj_mb_v) <- colnames(otutable)
  ig.mb     <- graph_from_adjacency_matrix(adj_mb,mode="undirected",diag=FALSE)
  return(list(res.net=ig.mb,res.spiec.easi=se.mb.amgut,otutable=otutable))
}

sparcc_network <- function(otutable,Transpose=FALSE,iter=20,inner_iter=10,th=0.3)
{
  otutable <- as.matrix(otutable)
  if(Transpose) {otutable <- t(otutable)}
  sparcc.amgut <- sparcc(otutable,iter=iter,inner_iter = inner_iter,th=th)
  adj_sparcc <- abs(sparcc.amgut$Cor) >= sparcc_th
  diag(adj_sparcc) <- 0
  adj_sparcc <- as.matrix(Matrix(adj_sparcc, sparse=TRUE))
  adj_sparcc_v <- adj_sparcc*sparcc.amgut$Cor
  colnames(adj_sparcc_v) <- colnames(OTU)
  rownames(adj_sparcc_v) <- colnames(OTU)
  colnames(adj_sparcc) <- colnames(OTU)
  rownames(adj_sparcc) <- colnames(OTU)
  ig.sparcc <- graph_from_adjacency_matrix(adj_sparcc,mode="undirected",diag=FALSE)
  return(list(res.net=ig.sparcc,res.sparcc=sparcc.amgut,otutable=otutable))
}
construct_spiecEasi_network <- function(se.mb.amgut,mb_th=0.3)
{
  adj_mb=as.matrix(symBeta(getOptBeta(se.mb.amgut), mode='maxabs'))
  colnames(adj_mb) <- colnames(otutable)
  rownames(adj_mb) <- colnames(otutable)
  adj_mb <- abs(adj_mb) > mb_th
  diag(adj_mb) <- 0
  adj_mb_v <- adj_mb*as.matrix(symBeta(getOptBeta(se.mb.amgut), mode='maxabs'))
  colnames(adj_mb_v) <- colnames(otutable)
  rownames(adj_mb_v) <- colnames(otutable)
  ig.mb     <- graph_from_adjacency_matrix(adj_mb,mode="undirected",diag=FALSE)
  return(list(res.net=ig.mb,res.spiec.easi=se.mb.amgut))
}

construct_sparcc_network <- function(sparcc.res,th)
{
  C <- sparcc.res$Cor
  C[abs(C) < th] <- 0
  diag(C) <- 0
  C.adj <- C != 0
  C.graph.with_d0 <- graph_from_adjacency_matrix(C.adj,mode="undirected",diag=FALSE)
  Cs <- colSums(C)
  del <- which(Cs==0)
  if(length(del) !=0) C <- C[-del,-del]
  C.adj <- C != 0
  C.graph.a <- graph_from_adjacency_matrix(C.adj,mode="undirected",diag=FALSE)
  C.graph.a <- delete_d0(C.graph.a)
  C.adj <- (C > 0)
  C.graph.p <- graph_from_adjacency_matrix(C.adj,mode="undirected",diag=FALSE)
  C.graph.p <- delete_d0(C.graph.p)
  C.adj <- (C < 0)
  C.graph.n <- graph_from_adjacency_matrix(C.adj,mode="undirected",diag=FALSE)
  C.graph.n <- delete_d0(C.graph.n)
  return(list(g.with_d0=C.graph.with_d0,g.a=C.graph.a,g.p = C.graph.p,g.n= C.graph.n,Cor=C))
}
delete_d0 <- function(g)
{
  d0 <- names(which(degree(g) == 0))
  if(length(d0) != 0) g <- delete_vertices(g,d0)
  return(g)
}

#network path and phylogenetic distance
net_phy_dis <- function(net,phy_tre)
{
  net_dis.mat <- distances(net,V(net))
  ids <- rownames(net_dis,mat)
  phy_dis.mat <- cophenetic(phy_tre)
  phy_dis.mat <- phy_dis.mat[ids,ids]
  net_dis.edgelist <- sym_matrix_to_list(net_dis.mat)
  net_dis.edgelist$phy_dis <- extract_target_from_matrix(phy_dis.mat,net_dis.edgelist,1,2)[,1]
  return(list(net_dis.mat=net_dis.mat,phy_dis.mat=phy_dis.mat,path2.list=net_dis.edgelist))
}

#extract targets from matrix
extract_target_from_matrix <- function(matrix1,extract_target,n=1,m=2)
{
  target_numbers <- nrow(extract_target)
  rownames_matrix <- rownames(matrix1)
  colnames_matrix <- colnames(matrix1)
  target_vector <- data.frame(NULL)
  for(i in 1:target_numbers)
  {
    row_rank <- which(rownames_matrix == extract_target[i,n])
    col_rank <- which(colnames_matrix == extract_target[i,m])
    target_vector[i,1] <- matrix1[row_rank,col_rank]
    rownames(target_vector)[i] <- paste(extract_target[i,n],extract_target[i,m],sep="-")
  }
  return(target_vector)
}

#convert symmetric matrix to list
sym_matrix_to_list <- function(Symmetric_matrix)
{
  k <- nrow(Symmetric_matrix)
  mat.names <- colnames(Symmetric_matrix)
  target <- Symmetric_matrix[lower.tri(Symmetric_matrix)]
  data1 <- data.frame(NULL)
  l <- 1
  for(i in 1:(k-1))
  {
    for(j in (i+1):k)
    {
      data1[l,1] <- mat.names[i]
      data1[l,2] <- mat.names[j]
      l <- l+1
    }
  }
  cbind(data1,target)
  return(data1)
}

three_dis_mats_2_list <- function(N,P,H,cn,z_scores_4_P=FALSE)
{
  ob <- colnames(N)
  sub.P <- P[ob,ob]
  sub.H <- H[ob,ob]
  lower.N <- N[lower.tri(N)]
  lower.sub.P <- sub.P[lower.tri(sub.P)]
  lower.sub.H <- sub.H[lower.tri(sub.H)]
  sub.p.vector <- as.vector(lower.sub.P)
  res <- cbind.data.frame(as.vector(lower.N),sub.p.vector,as.vector(lower.sub.H))
  n <- nrow(N)
  N1 <- NULL
  N2 <- NULL
  for(i in 1:(n-1))
  {
    N1 <- c(N1,rep(ob[i],n-i))
    N2 <- c(N2,ob[(i+1):n])
  }
  res <- cbind.data.frame(N1,N2,res,stringsAsFactors =FALSE)
  res <- res[!is.infinite(res[,3]),]
  if(z_scores_4_P) {res$sub.p.vector <- (res$sub.p.vector-min(res$sub.p.vector))/(max(res$sub.p.vector)-min(res$sub.p.vector))}
  colnames(res) <- cn
  return(res)
}
two_dis_mats_2_list <- function(P,S,cn)
{
  ob <- colnames(P)
  sub.S <- S[ob,ob]
  lower.P <- P[lower.tri(P)]
  lower.sub.S <- sub.S[lower.tri(sub.S)]
  res <- cbind.data.frame(as.vector(lower.P),as.vector(lower.sub.S))
  colnames(res) <- cn
  return(res)
}
#otunames
extract_otunames <- function(otutable,n=1)
{
  if(n==1)
  {otu_names <- rownames(otutable)} 
  if(n==2)
  {otu_names <- colnames(otutable)}
  df <- data.frame("Otu_names"=otu_names,'Seq_names'=otu_names)
  return(df)
}
###optima###
cal.env.optima = function(otu, meta.used, variables = c("elevation"),z_scores=FALSE) {
  
  # inputs:
  # 1, otu: otu table
  # 2, meta.used: meta matrix
  # 3, variables: variable selected
  otu[is.na(otu)] <- 0
  
  if (length(variables) == 1) { variables = c(variables, variables)  } 
  
  meta.used = meta.used[, variables]
  
  #*******************************************************************************
  # 1, calculate mean habitat value  #####
  #*******************************************************************************
  
  otu.names = colnames(otu); # otu.names[1:3]; # just the names of otus
  
  env.vars = colnames(meta.used)[1:ncol(meta.used)]; env.vars; ## a vector of names for environmental variables
  
  dim(otu); dim(meta.used); 
  otu.meta = merge(otu, meta.used, by=0); # otu.meta[1:5,c(1:5,(ncol(otu.meta)-5):ncol(otu.meta))]; dim(otu.meta); 
  ## merging otu table with meta data
  otu.meta <- na.omit(otu.meta)
  
  otu.mean.habitats = as.data.frame(matrix(c(-999), ncol=length(env.vars), nrow=length(otu.names)));  
  colnames(otu.mean.habitats) = env.vars; # otu.mean.habitats[1:2,1:2]; 
  # setting up an 'empty' matrix to fill with habitat values
  
  temporary = cbind(rep(1,nrow(otu.mean.habitats)),seq(1:nrow(otu.mean.habitats))); head(temporary)
  rownames(otu.mean.habitats) = paste(temporary[,1],temporary[,2],sep="_"); 
  # just putting in dummy row names that will be changed in the for loop that goes across otus
  
  for (i in 1:length(otu.names)) { 
    # for loop to calculate abundance weighted mean habitat values for each otu across each environmental variable
    
    otu.mean.habitats[i,] = apply(otu.meta[,otu.names[i]]*otu.meta[,env.vars],2,FUN='sum')/sum(otu.meta[,otu.names[i]]); 
    # the function multiplies the abundance vector to each environmental vector, 
    # takes the sum and then divides the sum by the total abundance of the otu
    
    rownames(otu.mean.habitats)[i] = otu.names[i]; 
    # changing the rowname just to be sure the correct otu is matched with the habitat values
    
    if(is.element(i,seq(0,length(otu.names),1000))==T) {print(i)} else{}; 
    # track progress of the for loop
  }; rm(i) ;
  del <- which(is.nan(rowSums(otu.mean.habitats)))
  if(length(del) != 0 ) otu.mean.habitats <- otu.mean.habitats[-del,]
  if(z_scores)
  {
    otu.mean.habitats.z_scores <- as.data.frame(apply(otu.mean.habitats,2,FUN=function(x){res=(x-min(x))/(max(x)-min(x))}))
    colnames(otu.mean.habitats.z_scores) <- variables
    return(otu.mean.habitats.z_scores)  
  }
  return(otu.mean.habitats)
}


#Mean neighbor XX distance
mnxd <- function(net,dis,delete.aloned_nodes=FALSE)
{
  if(delete.aloned_nodes){
    node_ids <- as_ids(V(net))
    no_adj_nodes <- unlist(lapply(adjacent_vertices(net,V(net)),length))
    delete.nodes <-node_ids[which(no_adj_nodes == 0)]
    net <- delete.vertices(net,delete.nodes)
  }
  nodes <- as_ids(V(net))
  res <- NULL
  adj_nodes <- adjacent_vertices(net,nodes)
  for(i in 1:length(nodes))
  {
    adj_nodes.temp <-  as_ids(adj_nodes[[i]])
    if(length(adj_nodes.temp) != 0) 
    {
      res[i] <- mean(dis[nodes[i],adj_nodes.temp])
    }else{ res[i] <- NA}
  }
  names(res) <- nodes
  return(res)
}
library(foreach)
library(doParallel)
mnxd.ses <- function(net,dis,per=999,delete.aloned_nodes=FALSE,cpu=2)
{
  if(delete.aloned_nodes){
    node_ids <- as_ids(V(net))
    no_adj_nodes <- unlist(lapply(adjacent_vertices(net,V(net)),length))
    delete.nodes <-node_ids[which(no_adj_nodes == 0)]
    if(length(delete.nodes)) net <- delete.vertices(net,delete.nodes)
  }
  cl <- makeCluster(cpu)
  clusterExport(cl,c("mnxd"))
  registerDoParallel(cl)
  n <- seq(1,vcount(net))
  names(n) <- as_ids(V(net))
  dis2 <- dis[names(n),names(n)]
  colnames(dis2) <- n[names(n)]
  rownames(dis2) <- n[names(n)]
  rd_mnxd <- foreach(i = 1:per,.packages=c("igraph","parallel","doParallel"),
                     .combine = cbind) %dopar%
                     {
                       rd_net <- sample_gnm(vcount(net),ecount(net))
                       mnxd.temp <- mnxd(rd_net,dis2)
                       return(mnxd.temp)
                     }
  stopCluster(cl)
  rd_mean <- apply(as.matrix(rd_mnxd),1,mean,na.rm=TRUE)
  rd_st <- apply(as.matrix(rd_mnxd),1,sd,na.rm=TRUE)
  mnxd.temp.obs <- mnxd(net,dis)
  rd_mean <- rd_mean[n[names(mnxd.temp.obs)]]
  rd_st <- rd_st[n[names(mnxd.temp.obs)]]
  nxd.res <- as.data.frame(as.matrix((mnxd.temp.obs-rd_mean)/rd_st,ncol=1))
  rownames(nxd.res) <- names(mnxd.temp.obs)
  return(nxd.res)
}


#convert symmetric matrix to list
sym_matrix_to_list <- function(Symmetric_matrix)
{
  k <- nrow(Symmetric_matrix)
  mat.names <- colnames(Symmetric_matrix)
  target <- Symmetric_matrix[lower.tri(Symmetric_matrix)]
  data1 <- data.frame(NULL)
  l <- 1
  for(i in 1:(k-1))
  {
    for(j in (i+1):k)
    {
      data1[l,1] <- mat.names[i]
      data1[l,2] <- mat.names[j]
      l <- l+1
    }
  }
  cbind(data1,target)
  return(data1)
}
#convert abundance table to binary table
convert_abund_to_binary <- function(abund_table)
{
  abund_table <- as.matrix(abund_table)
  abund_table[which(abund_table !=0)] <- 1
  return(abund_table)
}
#randomforest explanation
rf.func <- function(res)
{
  tdm2l <-  three_dis_mats_2_list(distances(graph_from_adjacency_matrix(res$stable_adj)),
                                  res$phy_dis,res$habitat_dis,
                                  cn=c("node1","node2","net_dis","phy_dis","habitat_dis"))
  rt.tdm2l <- randomForest(as.factor(net_dis)~phy_dis+habitat_dis,
                           data=tdm2l,importance=TRUE)
  return(rt.tdm2l)
}
rf_explanation <- function(rf_res)
{
  error.rate <- rf_res$confusion[,"class.error"]
  explaination <- 1-error.rate
  n <- ncol(rf_res$importance)
  res <- t(rf_res$importance[,-c(n-1,n)])
  res <- ifelse(res >0,res,0)
  res.rowsum <- rowSums(res)
  res <- res/res.rowsum*explaination
  res[is.nan(res)] <- 0
  res <- as.data.frame(res)
  res$net_dis <- rownames(res)
  return(res)
}
rf.func2 <- function(res)
{
  rf.res <- rf.func(res)
  rf.res <- rf_explanation(rf.res)
  return(rf.res)
}
rf.func3 <- function(res)
{
  res.tmp <- lapply(res,rf.func2)
  res.tmp <- melt(res.tmp,id.vars=c("net_dis"))
  return(res.tmp)
}

in_which_module <- function(Ver,Modules)
{
  n <- length(Modules)
  for(i in 1:n)
  {
    if(Ver %in% Modules[[i]]) return(i)
  }
}

lm.p <- function(x,y){
  res.lm <- lm(y~x)
  res.lm.sum <- summary(res.lm)
  n.r <- nrow(res.lm.sum$coefficients)
  if(n.r == 1) return(1)
  if(n.r != 1) return(res.lm.sum$coefficients[2,4])}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

position.nmxd <- function(MotifPos,nmxd)
{
  tmp <- merge(MotifPos,nmxd,by=0,all.x=T)
  del <- which( is.na(tmp$degree) == T)
  if(length(del) != 0) tmp <- tmp[-del,]
  tmp1 <- tmp[,2:149]
  tmp2 <- tmp[,150:151]
  tmp.nmpd <- colSums(tmp1*tmp2$mnpd.ses)
  tmp.nmhd <- colSums(tmp1*tmp2$mnhd.ses)
  res.tmp <-  as.data.frame(cbind(tmp.nmpd,tmp.nmhd))
  colnames(res.tmp) <- c("mnpd","mnhd")
  res.tmp$Position <- seq(1,148)
  del <- which(((res.tmp$mnpd==0) + (res.tmp$mnhd==0))!=0)
  if(length(del) != 0) res.tmp <- res.tmp[-del,]
  return(res.tmp)
}

match_two_mat <- function(x,a)
{
  an <- colnames(a)
  x.tmp <- x[an,an]
  return(x.tmp)
}

BICg_VB <- function(Net,d=2,maxg=5,lcc=T)
{
  data1 <- data.frame(NULL)
  if(maxg < 1) {maxg=1;print("maxp is less than 1, it is setted to be 1 automatically")}
  for(i in 1:maxg)
  {
    print(paste("start caculating with dimension of",i,sep=" "))
    samp.fit <- vblpcmstart(Net,G=i,d=d,lcc=lcc,LSTEPS = 1000,CLUST=0)
    data1[i,1] <- i
    data1[i,2] <- samp.fit$BIC$overall
    data1[i,3] <- samp.fit$BIC$Y[1,1]
    data1[i,4] <- samp.fit$BIC$MBC
  }
  colnames(data1) <- c("No.Groups","Overall BIC","Likelihood BIC","Latent space/clustering BIC")
  return(data1)
}
motif_table <- function(x)
{
  N <- length(x)
  res <- matrix(0,nrow=N,ncol=44)
  colnames(res) <- seq(1,44)
  rownames(res) <- names(x)
  for(i in 1:N)
  {
    res[i,] <- x[[i]]$normalise_sum
  }
  return(res)
}

ls_z <- function(ls.res)
{
  Z <- ls.res$V_z
  rownames(Z) <- network.vertex.names(ls.res$net)
  colnames(Z) <- c("z1","z2")
  return(Z)
}

ls_Z_dist <- function(x,Met="eucli")
{
  tmp <- as.matrix(vegdist(x,method=Met))
  tmp <- (tmp-min(tmp))/(max(tmp)-min(tmp))
  return(tmp)
}

Z_dist <- function(x)
{
  N <- ncol(x)
  z.dist <- list()
  z.dist[[1]] <- ls_Z_dist(x,Met="eucli")
  for(i in 1:N)
  {
    z.dist[[i+1]] <- ls_Z_dist(x[,i],Met="manhattan")
  }
  names(z.dist) <- c("Z.overall.eucli",paste("z",seq(1:N),".manhattan",sep=""))
  return(z.dist)
}

ls_G <- function(x)
{
  tmp <- t(x$V_lambda)
  rownames(tmp) <- network.vertex.names(x$net)
  return(tmp)
}

z_scores_4_c <- function(x,vs)
{
  if(length(vs)==0) stop("Specifying the column index")
  N <- length(vs)
  for(i in 1:N)
  {
    tmp <- vs[i]
    tmp.x <- x[,tmp]
    x[,tmp] <- (tmp.x-min(tmp.x))/(max(tmp.x)-min(tmp.x))
  }
  return(x)
}
cor.test.p <- function(x,y)
{
  if(length(x)>=10)
  {
    tmp <- cor.test(x,y)
    tmp$p.value
  }else{
    return(NA)
  }
}
glm.res <- function(y,x1,x2){
  if(length(x1)>=10) 
  {
    res <- glm(y~x1+x2)
    res <- summary(res)$coefficients
    return(res)
  }else{
    return(matrix(NA,nrow=3,ncol=4))
  }
}

#sum.j(Fij * Pij)/sum.j(Pij)
mat.dot <- function(mat.f,mat.p,col.index)
{
  mat.names <- rownames(mat.f)
  mat.p.match <- t(as.matrix(mat.p[mat.names,col.index]))
  res <- mat.p.match %*% mat.f /colSums(mat.f)
  return(res)
}

#classify for MNXD.ses
classify_mnxd.1 <- function(x)
{
  if(x >=2) return(1)
  if(x >-2 && x < 2) return(2)
  if(x<= -2) return(3)
}
classify_mnxd <- function(mnpd,mnhd)
{
  index.mnpd <- 4-sapply(mnpd,classify_mnxd.1)
  index.mnhd <- sapply(mnhd,classify_mnxd.1)
  m <- t(matrix(seq(1,9),3,3))
  n <- length(mnpd)
  tmp <- NULL
  for(i in 1:n)
  {
    tmp <- c(tmp,m[index.mnhd[i],index.mnpd[i]])
  }
  return(tmp)
}

pos_pro=function(mat,trn=F)
{
  if(trn) mat=t(mat)
  Nc=nrow(mat)
  for(i in 1:Nc)
  {
    tmp=mat[i,]
    if(sum(tmp) != 0) tmp=tmp/sum(tmp)
    mat[i,]=tmp
  }
  if(trn) mat=t(mat)
  return(mat)
}