#1. species pools
#The weight across network path
#Wn = importance*(1+D)/(1+exp(rate*n-C))
#n is network path distance
get_weighted_across_net_dis <- function(rate=1,importance=1,D=0,n=6,C=5)
{
  res <- NULL
  for(i in 1:n)
  {
    res[i] <- importance*(1+D)/(1+exp(rate*(i-C)))
  }
  return(res)
}
#Generate community phylogeny with N species
#And scale with formula: phy_dis = (phy_dis-min(phy_dis))/(max(phy_dis)-min(phy_dis))
get_phylogeny_dis <- function(N,tree.br)
{
  res <- cophenetic(rtree(N,br=tree.br))
  tn <- paste("t",seq(1,N),sep="")
  res <- res[tn,tn]
  res <- (res-min(res))/(max(res)-min(res))
  return(res)
}
#Generate habitat distance according to community phylogeny with phylogeny signal
#       if phylogeny distance < phy_sig: habiat distance  = phy_sig*rate + e; 
#       else: habiat distance  = phy_sig*rate + U(lower_limit*phy_sig*rate, upper_limit*phy_sig*rate) + e;
#And scale with formula: habitat_dis = (habitat_dis-min(habitat_dis))/(max(habitat_dis)-min(habitat_dis))
get_habitat_dis <- function(phy_adj,phy_sig=0.3,rate=0.5,C=0.1,b1=0.9,b2=1.1,e1=0.01,e2=0.01)
{
  res <- phy_adj
  for(i in 1:nrow(phy_adj))
  {
    for(j in 1:ncol(phy_adj))
    {
      temp <- phy_adj[i,j]
      if(temp <= phy_sig && i!=j) res[i,j] <- rate*phy_adj[i,j] + rnorm(1,0,e1)+C
      if(temp > phy_sig && i!=j) res[i,j] <- runif(1,b1*(rate*phy_sig+C),b2*(rate*phy_sig+C))+rnorm(1,0,e2)
    }
  }
  res <- (res-min(res))/(max(res)-min(res))
  return(res)
}
#Generate phylogeny/habitat preference matrix
#phylogeny/habitat preference = 0.5*(preference.i+preference.j)/distance
get_preference_adj <- function(dis,unif_max=1,unif_min=0,Scales=2)
{
  N <- nrow(dis)
  preference <- runif(N,max=unif_max,min=unif_min)
  names(preference) <- rownames(dis)
  res <- dis
  for(i in 1:N)
  {
    for(j in 1:N)
    {
      if(i != j)
      {
        res[i,j] <- 0.5*(preference[i]+preference[j])/dis[i,j]
      }
    }
  }
  res <- (res-min(res))/(max(res)-min(res))*Scales
  return(res)
}
#Get netowrk path distance
get_net_dis <- function(adj_matrix)
{
  g <- graph_from_adjacency_matrix(adj_matrix,mode="undirected")
  dis_mat <- distances(g)
  return(dis_mat)
}
#Get indirected preference
#Indirected_preference (path=n) =  0.5*(sum_ki(preference(i,j,ki,n)+sum_kj(preference(j,i,kj,n))
# ki/kj are the neighbors of node i/j with n path distance
get_indirected_preference <- function(ni,nj,preference_mat,net_dis.mat,path_N=5)
{
  pre1 <- NULL
  pre2 <- NULL
  pre <- NULL
  ni_pos <- which(colnames(net_dis.mat)==ni)
  nj_pos <- which(colnames(net_dis.mat)==nj)
  if(path_N == 0) return(0)
  for(i in 1:path_N)
  {
    adj.nodes.ni.temp <- get_node_ids(net_dis.mat[ni,-nj_pos],path_N)
    adj.nodes.nj.temp <- get_node_ids(net_dis.mat[nj,-ni_pos],path_N)
    pre1[i] <- mean(preference_mat[ni,adj.nodes.nj.temp])
    pre2[i] <- mean(preference_mat[nj,adj.nodes.ni.temp])
    if(is.nan(pre1[i])) pre1[i]=0
    if(is.nan(pre2[i])) pre2[i]=0
    pre[i] <- 0.5*(pre1[i]+pre2[i])
  }
  return(pre)
}
#Get directed preference
get_directed_preference <- function(ni,nj,prefernce_mat)
{
  return(prefernce_mat[ni,nj])
}
#scale free properties
node_degree <- function(n,adj_mat)
{
  g <- graph_from_adjacency_matrix(adj_mat)
  N <- length(adjacent_vertices(g,n)[[1]])
  return(N)
}
#Get node ids
get_node_ids <- function(V,path_N)
{
  names(which(V==path_N))
}
#Get probability matrix of connectivity
#Connetivity strength = 
#     Wp/sum(Wp+Wh)*(degree.i+degree.j+2)^(-2+Wp*Phylogeny_preference) + 
#     Wh/sum(Wp+Wh)*(degree.i+degree.j+2)^(-2+Wh*habitat_preference)
connetivity_strength <- function(ni,nj,adj_mat,phy_preference.mat,
                                 habitat_preference.mat,net_dis.mat,
                                 Wp,Wh,path_N,a,ori.st=2)
{
  Ppre.directed.ij <- get_directed_preference(ni,nj,prefernce_mat =phy_preference.mat)
  Hpre.directed.ij <- get_directed_preference(ni,nj,prefernce_mat =habitat_preference.mat)
  Ppre.indirected.ij <- get_indirected_preference(ni,nj,preference_mat = phy_preference.mat,net_dis.mat,path_N)
  Hpre.indirected.ij <- get_indirected_preference(ni,nj,preference_mat =habitat_preference.mat,net_dis.mat,path_N)
  Ppre.ij <- c(Ppre.directed.ij,Ppre.indirected.ij)
  Hpre.ij <- c(Hpre.directed.ij,Hpre.indirected.ij)
  Wp.scale <- Wp/sum(Wp+Wh)
  Wh.scale <- Wh/sum(Wp+Wh)
  Wp.scale.b <- Wp/sum(Wp)
  Wh.scale.b <- Wh/sum(Wh)
  Wp.Ppre.ij <- Wp.scale.b*Ppre.ij
  Wh.Hpre.ij <- Wh.scale.b*Hpre.ij
  node.degree.i <- node_degree(ni,adj_mat)
  node.degree.j <- node_degree(nj,adj_mat)
  degr <- (node.degree.i+node.degree.j+ori.st)/2
  P.connected.ij <- sum(Wp.scale*degr^(-a+Wp.Ppre.ij)) + sum(Wh.scale*degr^(-a+Wh.Hpre.ij))
  return(P.connected.ij)
}

#Initial network generated from Binomial experiments with a probability
initial_network <- function(adj_mat,pro)
{
  N <- nrow(adj_mat)
  res <- adj_mat
  while(sum(res)/2 < pro)
  {
    temp <- sample(seq(1,N),2)
    if(temp[1] != temp[2]) {res[temp[1],temp[2]] = 1;res[temp[2],temp[1]] = 1}
  }
  return(res)
}
#Get unconnected nodes of node i.
get_unconnected_nodes <- function(ni,adj_mat)
{
  g <- graph_from_adjacency_matrix(adj_mat)
  V_ids <- as_ids(V(g))
  connected_nodes <- which(V_ids %in% as_ids(adjacent_vertices(g,ni)[[1]]))
  if( length(connected_nodes) !=0) {unconnected_nodes <- V_ids[-connected_nodes]}else{
    unconnected_nodes <- connected_nodes
  }
  return(unconnected_nodes)
}
#Get connected nodes
get_connected_nodes <- function(adj_mat)
{
  g <- graph_from_adjacency_matrix(adj_mat)
  V_ids <- as_ids(V(g))
  V_g <- degree(g)
  connected_nodes <- V_ids[which(V_g!=0)]
  return(connected_nodes)
} 
#Get connectivity strength matrix
get_connective_strength_mat <- function(adj_mat,phy_preference.mat,
                                        habitat_preference.mat,Wp,Wh,path_N,a,connect_p,ori.st)
{
  N <- ncol(adj_mat)
  Species <- rownames(adj_mat)
  res <- matrix(0,nrow=N,ncol=N)
  colnames(res) <- Species
  rownames(res) <- Species
  for(i in 1:(N-1))
  {
    for(j in (i+1):N)
    {
      ni <- Species[i]
      nj <- Species[j]
      if(adj_mat[i,j] ==1) {
        res[i,j] <- connetivity_strength(ni,nj,adj_mat,phy_preference.mat,
                                         habitat_preference.mat,get_net_dis(adj_mat),Wp,Wh,path_N,a,ori.st)
        res[j,i] <- res[i,j]
      }
    }
  }
  return(res)
}
#Get binary matrix with Binomial experiment 
get_binary_matrix <- function(mat)
{
  res <- mat
  n <- nrow(mat)
  m <- ncol(mat)
  for(i in 1:n)
  {
    for(j in 1:m)
    {
      res[i,j] <- rbinom(1,1,mat[i,j])
      res[j,i] <- res[i,j]
    }
  }
  return(res)
}

#main function
#1. Construt species pool
#   1.1 community phylogeny
#   1.2 phylogeny/habitat preferece matrix
#2. Get wight of Phylogeny/habitat across network path. Wp,Wh
#3. Initialize network: return adjacency matrix
#4. Loop body
#  while
#     4.1 if pro_disturbance 
#           Get new Wp and Wh
#     4.2 calculate the all connectivity strength of unconnected links between species
#     4.3 keep the unconnected links with top of connectivity strength ready to connect
#     4.4 get the connection matrix with Binomial experiment
#     4.5 transitive adjacency matrix = adjacency matrix + connection matrix
#     4.6 transitive connectivity strenght matrix with transitive adjacency matrix
#     4.7 the links with low connectivity strenght would disconnect with a probability
#         and get disconnection matrix
#     4.8 adjacency matrix = transitive adjacency matrix - disconnection matrix
stimulation_main <- function(Nspecies,rounds=50,phy_sig=0,tree.br,
                             rate=c(1.8,2),importance=c(2,1),pre_D=0,n=6,a=3,
                             pro=0.01,th_ratio=0.95,connect_th,dis_th=5e-2,pro_disturbance=FALSE,
                             pro_disturbance_round=10,pro_D=1,connect_p=0.1,dis_p=0.05,ori.st=2,C,Samples=NA,plot.graph=TRUE)
{
  done <- 0
  oop <- 0
  pro_dis <- NA
  if(length(Samples) ==1) {ifelse(is.na(Samples),Samples = rounds,Samples=Samples)}
  adj_mat.list <- list(NULL)
  length(adj_mat.list) <- length(Samples)
  list.n=1
  names(adj_mat.list) <- paste("round",Samples,sep="_")
  GL <- NULL
  LL <- NULL
  TL <- NULL
  #1. Construt species pool
  # 1.1 ====
  phy_dis.mat <- get_phylogeny_dis(Nspecies,tree.br)
  Species <- rownames(phy_dis.mat)
  # 1.2 ====
  habitat_dis.mat <- get_habitat_dis(phy_dis.mat,phy_sig)
  phy_preference.mat <- get_preference_adj(dis = phy_dis.mat)
  habitat_preference.mat <- get_preference_adj(dis = habitat_dis.mat)
  #2. Get wight of Phylogeny/habitat across network path. Wp,Wh
  Wp <- get_weighted_across_net_dis(rate = rate[1],importance = importance[1],D = 0,n = n,C)
  Wh <- get_weighted_across_net_dis(rate = rate[2],importance = importance[2],D = pre_D,n = n,C)
  t_rounds=rounds
  #3. Initialize network: return adjacency matrix
  net_adj.mat <- matrix(0,nrow=Nspecies,ncol=Nspecies)
  rownames(net_adj.mat) <- Species
  colnames(net_adj.mat) <- Species
  adj_mat <- initial_network(net_adj.mat,pro=pro)
  # 4. Loop body
  while(rounds > 0)
  {
    # 4.1 Pro disturbance
    print(paste("==== Round: ",t_rounds-rounds+1," ====",sep=""))
    
    if(pro_disturbance && done==1 && oop==0)
    {
      print("Oops! Disturbance coming..")
      oop <- oop+1
      pro_dis <- t_rounds-rounds+1
      Wp <- get_weighted_across_net_dis(rate = rate[1],importance = importance[1],D = 0,n = n,C)
      Wh <- get_weighted_across_net_dis(rate = rate[2],importance = importance[2],D = pro_D,n = n,C)
    }
    # 4.2 calculate the all connectivity strength of unconnected links between species
    connected_nodes <- get_connected_nodes(adj_mat)
    connection.mat <- matrix(0,nrow=Nspecies,ncol=Nspecies)
    rownames(connection.mat) <- Species
    colnames(connection.mat) <- Species
    for(i in 1:length(connected_nodes))
    {
      node.i <- connected_nodes[i]
      unconnected_nodes.i <- get_unconnected_nodes(node.i,adj_mat)
      C.temp <- NULL
      for(j in 1:length(unconnected_nodes.i))
      {
        node.j <- unconnected_nodes.i[j]
        if(node.i != node.j){
          connection.mat[node.i,node.j] <- connetivity_strength(node.i,node.j,adj_mat,phy_preference.mat,
                                                                habitat_preference.mat,get_net_dis(adj_mat),Wp,Wh,n-1,a,ori.st=ori.st)
          connection.mat[node.j,node.i] <- connection.mat[node.i,node.j]
        }
      }
    }
    # 4.3 keep the unconnected links with top of connectivity strength ready to connect
    connection.mat.na <- connection.mat
    connection.mat.na <- connection.mat.na[connection.mat.na!=0]
    th <- quantile(connection.mat.na,th_ratio,na.rm = TRUE)
    connection.mat <- ifelse(connection.mat >=1,1,connection.mat)
    connection.mat <- ifelse(connection.mat >=max(th,connect_th),max(connection.mat,connect_p),0)
    print(paste("Max connectivity strength: ",max(connection.mat.na),sep=""))
    # 4.4 get the connection matrix with Binomial experiment
    connection.mat <- get_binary_matrix(connection.mat)
    # 4.5 adjacency matrix = adjacency matrix + connection matrix
    adj_mat.temp <- adj_mat + connection.mat
    # 4.6 transitive connectivity strenght matrix with adjacency matrix
    pro_connection_strength.mat <- get_connective_strength_mat(adj_mat.temp,phy_preference.mat,habitat_preference.mat,Wp,Wh,n-1,a,ori.st=0)
    pro_connection_mat <- adj_mat*pro_connection_strength.mat
    disconnection.mat <- pro_connection_mat
    # 4.7 the links with low connectivity strenght would disconnect with a probability
    #         and get disconnection matrix
    disconnection.mat <- ifelse((disconnection.mat !=0)*(disconnection.mat < dis_th) == 1,dis_p,0)
    PLL <- sum(disconnection.mat !=0)/2
    print(paste("Potential loss Links: ",PLL,sep=""))
    disconnection.mat <- get_binary_matrix(disconnection.mat)
    #4.8 adjacency matrix = adjacency matrix - disconnection matrix
    adj_mat <- adj_mat.temp - disconnection.mat
    rounds=rounds-1
    #=============== basic dynamic =================#
    get_links <- sum(connection.mat)/2
    loss_links <- sum(disconnection.mat)/2
    print(paste("Wp/(Wp+Wh): ",sum(Wp)/(sum(Wh)+sum(Wp)),sep=""))
    print(paste("Get: ",round(get_links,5)," links."," Loss: ",round(loss_links,5)," links.",sep=""))
    if(max(connection.mat.na) < connect_th && PLL == 0) {print("Assamblage done.");stable_adj=adj_mat;stable_rounds=t_rounds-rounds;done <- done +1}
    if(!pro_disturbance && done == 1) break
    if(pro_disturbance && done == 2) break
    #=============== plot function =================#
    if(plot.graph){
      g <- graph_from_adjacency_matrix(adj_mat,mode="undirected")
      m <- round(modularity(fastgreedy.community(g)),2)
      plot(g,main=paste("Round: ",t_rounds-rounds,"; Link: ",sum(adj_mat)/2,"; Modu: ",m,sep=""),vertex.label="")
    }
    #============== Sample the rounds ==============#
    if((t_rounds-rounds) %in% Samples) {adj_mat.list[[list.n]] = adj_mat;list.n=list.n+1;GL <- c(GL,get_links);LL <- c(LL,loss_links);TL <- c(TL,sum(adj_mat)/2)}
  }
  return(list(pro_dis=pro_dis,GL=GL,LL=LL,TL=TL,phy_dis = phy_dis.mat,habitat_dis = habitat_dis.mat,phy_prefernece=phy_preference.mat,habitat_preference=habitat_preference.mat,stable_adj=adj_mat,stable_rounds=t_rounds-rounds,adj_mat=adj_mat.list))
}