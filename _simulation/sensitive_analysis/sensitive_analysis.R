library(igraph)
library(doParallel)
library(ape)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./sensitive_analysis/simulation_function.R")
#===============================================#
Nspecies <- 50     #Number of species
rounds <- 200      #Iterations
phy_sig <- 0.15     #phylogeny sinal eg. [0,1]
l <- 1             #Distance of network path between interacting speices, e.g. 1,2,3,4.. 
# if l = 1, there merely are directed effects on speices interaction
a <- 2             #Constant in formula (1): 
#    pij = Wp/(Wp+Wh)*((di+dj+2)/2)^(-a+Wp*Ppre.ij) + Wh/(Wp+Wh)*((di+dj+2)/2)^(-a+Wh*Hpre.ij)
#about 2 
e.t0 <- 3        #the links for initailing network, e.g. 1,2,3,4 
pc.th <- 0.99   #Top percentage of preparing connectivity strength
pc.con.min <- 0.12  #The minimum connectivity strength
pc.dis.min <- 0.06     #Disconnected threshold of pij
p.con.min <- pc.con.min  #The minimum probability of connection
p.dis.min <- 0.06      #The probability of disconnection
pro_Eh_yn <- FALSE  #Add pro-enviromental heterogeneity or not (TRUE/FALSE)
Eh <- 0         #Add pre-enviromental heterogeneity firstly in formula (2):
# Wn = A*(1+Eh)/(1+exp(B*(n-C)))
pro_Eh_round <- 101 #Add pro-disturbance in which round
pro_Eh <- 1    #Add pro-enviromental heterogeneity at round pro_Eh_round
ori.k <- 2               #The constant in formula (1) [(di+dj+ori.k)/2] 
A <- c(0.2,0.1)      #The constant in formula (2). 
B <- c(0.5,0.5)          #The constant in formula (2). 
C <- 0
Samples <- seq(1,rounds,100)  #The rounds are sampled. a vector. 
plot.graph <- FALSE        #plot the network or not
tree.br <- runif    #an R function used to generate the branch lengths 

#=================
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

load("./sensitive_anlysis_meta.Rdata")

#=================
set.seed(111)
Per <- 100
cpu=2
cl <- makeCluster(cpu)
registerDoParallel(cl)
variable_Nspecies <- list()

for (j in 1:nrow(vars)) {
  Nspecies <- vars[j,1]
  l <- vars[j,2]
  pc.con.min <- vars[j,3]
  p.con.min <- vars[j,3]
  pc.dis.min <- vars[j,4]
  p.con.min <- vars[j,4]
  a <- vars[j,5]
  pc.th <- vars[j,6]
  phy_sig <- vars[j,7]
  print(j)
  variable_Nspecies[[j]] <- foreach(i = 1:Per,.packages=c("igraph","parallel","doParallel","ape")) %dopar%
  {
    stimulation_main(Nspecies,rounds,phy_sig,tree.br,B,A,pre_D,l,a,
                     e.t0,pc.th,pc.con.min,pc.dis.min,pro_Eh_yn,pro_Eh_round,
                     pro_Eh,p.con.min,p.dis.min,ori.k,C,Samples,plot.graph)
  }
}
stopCluster(cl)
#=====================================
