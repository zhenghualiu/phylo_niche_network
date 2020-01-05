library(igraph)
library(doParallel)
library(ape)

source("simulation_function.R")
#===============================================#
rounds <- 200      #Iterations
phy_sig <- 0.15     #phylogeny sinal eg. [0,1] 
pro <- 3        #the links for initailing network, e.g. 1,2,3,4 
pro_disturbance <- FALSE  #Add pro-disturbance or not (TRUE/FALSE)
pre_D <- 0         #Add pre-distance firstly in formula (2):
# Wn = importance*(1+pre_D)/(1+exp(rate*(n-C)))
pro_disturbance_round <- 101 #Add pro-disturbance in which round
ori.st <- 2               #The constant in formula (1) [(di+dj+ori.st)] 
rate <- c(0.5,0.5)          #The constant in formula (2). 
C <- 0
importance <- c(0.2,0.1)      #The constant in formula (2). 
Samples <- seq(1,rounds,100)  #The rounds are sampled. a vector. 
plot.graph <- FALSE        #plot the network or not
tree.br <- runif    #an R function used to generate the branch lengths 

#Sensitive analysis
set.seed(111)
Per <- 100
cpu=40
cl <- makeCluster(cpu)
registerDoParallel(cl)
#==========
#N <- seq(50,100,1)
#n <- seq(1,5,1)
#connection <- seq(0.08,0.12,0.01)
#disconnection <- seq(0.06,0.08,0.01)
#a <- seq(2,3,0.1)
#th <- seq(0.95,0.99,0.01)
#phy_sig <-seq(0.12,0.2,0.02)
#num_obs <- 100
#vars <- cbind(sample(N,num_obs,replace=TRUE),
#              sample(n,num_obs,replace=TRUE),
#              sample(connection,num_obs,replace=TRUE),
#              sample(disconnection,num_obs,replace=TRUE),
#              sample(a,num_obs,replace=TRUE),
#              sample(th,num_obs,replace=TRUE),
#              sample(phy_sig,num_obs,replace=TRUE)
#)
#save(vars,file="sensitive_anlysis_meta.Rdata")
load("sensitive_anlysis_meta.Rdata")
vars <- vars[81:90,]
print(paste("Len of vars: ",nrow(vars),sep=""))
variable_Nspecies <- list()
for (j in 1:nrow(vars)) {
  Nspecies <- vars[j,1]
  n <- vars[j,2]
  connect_p <- vars[j,3]
  connect_th <- vars[j,3]
  dis_p <- vars[j,4]
  dis_th <- vars[j,4]
  a <- vars[j,5]
  th_ratio <- vars[j,6]
  phy_sig <- vars[j,7]
  print(j)
  variable_Nspecies[[j]] <- foreach(i = 1:Per,.packages=c("igraph","parallel","doParallel","ape")) %dopar%
  {
    stimulation_main(Nspecies,rounds,phy_sig,tree.br,rate,importance,pre_D,n,a,
                     pro,th_ratio,connect_th,dis_th,pro_disturbance,pro_disturbance_round,
                     pro_D,connect_p,dis_p,ori.st,C,Samples,plot.graph)
  }
}
stopCluster(cl)
save(variable_Nspecies,vars,file="sensitive_analysis_81_90.Rdata")