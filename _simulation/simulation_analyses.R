pro_df <- function(data)
{
  pro_dis.round <- unlist(lapply(data,function(x){x$pro_dis}))
  TL <- lapply(data,function(x){x$TL})
  adj_mat <- lapply(data,function(x){x$adj_mat})
  gr.res <- NULL
  lr.res <- NULL
  mr.res <- NULL
  for(i in 1:length(pro_dis.round))
  {
    pro_dis <- pro_dis.round[i]
    if(!is.na(p
              ro_dis))
    {
      s1 <- pro_dis-1
      s2 <- length(TL[[i]])
      adj1 <- adj_mat[[i]][[s1]]
      adj2 <- adj_mat[[i]][[s2]]
      gr.tmp <- gain_ratio(adj1,adj2)
      lr.tmp <- loss_ratio(adj1,adj2)
      mr.tmp <- module_changed_ratio(adj1,adj2)
      gr.res <- c(gr.res,gr.tmp)
      lr.res <- c(lr.res,lr.tmp)
      mr.res <- c(mr.res,mr.tmp)
    }
  }
  res <- cbind.data.frame(gr.res,lr.res,mr.res)
  colnames(res) <- c("Gain","Loss","M2/M1")
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

