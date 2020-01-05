library(igraph)
pro_df <- function(data)
{
  pro_dis.round <- unlist(lapply(data,function(x){x$pro_dis}))
  TL <- lapply(data,function(x){x$TL})
  adj_mat <- lapply(data,function(x){x$adj_mat})
  gr.res <- NULL
  lr.res <- NULL
  mr.res <- NULL
  md.res <- NULL
  for(i in 1:length(pro_dis.round))
  {
    pro_dis <- pro_dis.round[i]
    if(!is.na(pro_dis))
    {
      s1 <- pro_dis-1
      s2 <- length(TL[[i]])
      adj1 <- adj_mat[[i]][[s1]]
      adj2 <- adj_mat[[i]][[s2]]
      gr.tmp <- gain_ratio(adj1,adj2)
      lr.tmp <- loss_ratio(adj1,adj2)
      mr.tmp <- module_changed_ratio(adj1,adj2)
      md.tmp <- md(adj1,adj2)
      gr.res <- c(gr.res,gr.tmp)
      lr.res <- c(lr.res,lr.tmp)
      mr.res <- c(mr.res,mr.tmp)
      md.res <- c(md.res,md.tmp)
    }
  }
  res <- cbind.data.frame(gr.res,lr.res,mr.res,md.res)
  colnames(res) <- c("Gain","Loss","M2/M1","Modularity difference")
  return(res)
}

gain_ratio <- function(adj1,adj2)
{
  tmp.g <- sum((adj2-adj1) == 1)/2
  tmp.a <- sum(adj1)/2
  return(tmp.g/tmp.a)
}
loss_ratio <- function(adj1,adj2)
{
  tmp.l <- sum((adj2-adj1) == -1)/2
  tmp.a <- sum(adj1)/2
  return(tmp.l/tmp.a)
}

module_changed_ratio <- function(adj1,adj2)
{
  m1 <- modularity(cluster_fast_greedy(graph_from_adjacency_matrix(adj1,mode="undirected")))
  m2 <- modularity(cluster_fast_greedy(graph_from_adjacency_matrix(adj2,mode="undirected")))
  return(m2/m1)
}

#variation coefficient of modularity difference
md <- function(adj1,adj2)
{
  m1 <- modularity(cluster_fast_greedy(graph_from_adjacency_matrix(adj1,mode="undirected")))
  m2 <- modularity(cluster_fast_greedy(graph_from_adjacency_matrix(adj2,mode="undirected")))
  return(m2-m1)
}

my_aov <- function(y,x,n)
{
  x <- factor(x)
  tmp.aov <- summary(aov(y~x))
  res <- c(tmp.aov[[1]][1,4],tmp.aov[[1]][1,5])
  return(res[n])
}

cal_vc <- function(x){sd(x)/mean(x)}

my_bootstrap=function(df,v,l,m,theta,N)
{
  vs=unique(df[,v])
  res <- NULL
  lr <- NULL
  Nr <- NULL
  for(i in 1:length(vs))
  {
    df.tmp <- df[which(df[,v] == vs[i]),]
    l1 <- unique(df.tmp[,l])
    for(j in 1:length(l1))
      {
         df.tmp2 <- df[which(df.tmp[,l] == l1[j]),]
         sam <- df.tmp2[,m]
         res.tmp <- bootstrap(sam,N,theta)$thetastar
         res <- c(res,res.tmp)
         lr <- c(lr,rep(as.numeric(l1[j]),N))
         Nr <- c(Nr,rep(vs[i],N))
      }
  }
  res[is.nan(res)] <- 0
  re.res=cbind.data.frame(Nr,lr,res)
  return(re.res)
}