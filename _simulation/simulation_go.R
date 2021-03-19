library(igraph)
library(doParallel)
library(ape)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./simulation_function.R")
#===============================================#
Nspecies <- 50     #Number of species
rounds <- 200     #Iterations
phy_sig <- 0.15    #phylogeny sinal eg. [0,1]
l <- 3             #Distance of network path between interacting speices, e.g. 1,2,3,4.. 
                   # if l = 1, there merely are directed effects on speices interaction
                   #Commpared with reference, l in here has a lag of 1.
a <- 2             #Constant in formula (1): 
                   #cij = Wp/(Wp+Wh)*((di+dj+2)/2)^(-a+Wp*Ppre.ij) + Wh/(Wp+Wh)*((di+dj+2)/2)^(-a+Wh*Hpre.ij)
                   #about 2 
e.t0 <- 3               #the links for initailing network, e.g. 1,2,3,4 
pc.th <- 0.99           #Top percentage of preparing connectivity strength
pc.con.min <- 0.12      #The minimum connectivity strength
pc.dis.min <- 0.06      #Disconnected threshold of pij
p.con.min <- pc.con.min #The minimum probability of connection
p.dis.min <- 0.06       #The probability of disconnection
pro_Eh_yn <- FALSE      #Add pro-enviromental heterogeneity or not (TRUE/FALSE)
Eh <- 0                 #Add pre-enviromental heterogeneity firstly in formula (2):
                        # Wn = A*(1+Eh)/(1+exp(B*(n-C)))
pro_Eh_round <- NA      #Add pro-disturbance in which round
pro_Eh <- 1    #Add pro-enviromental heterogeneity at round pro_Eh_round
ori.k <- 2          #The constant in formula (1) [(di+dj+ori.k)/2] 
A <- c(0.2,0.1)     #The constant in formula (2). 
B <- c(0.5,0.5)     #The constant in formula (2). 
C <- 0              #The constant in formula (2).
Samples <- seq(1,rounds,1)  #The rounds are sampled. a vector. 
plot.graph <- TRUE           #plot the network or not
tree.br <- runif    #an R function used to generate the branch lengths 

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
  a <- vars[j,3]
  print(j)
  variable_Nspecies[[j]] <- foreach(i = 1:Per,.packages=c("igraph","parallel","doParallel","ape")) %dopar%
  {
    stimulation_main(Nspecies,rounds,phy_sig,tree.br,B,A,Eh,l,a,
                     e.t0,pc.th,pc.con.min,pc.dis.min,pro_Eh_yn,pro_Eh_round,
                     pro_Eh,p.con.min,p.dis.min,ori.k,C,Samples,plot.graph)
  }
}
stopCluster(cl)
#=====================================
rf_res.ttt <- lapply(variable_Nspecies[[1]],rf_res)
melt(rf_res.ttt,id.vars=c("variable","net_dis")) %>% ggplot(aes(x=net_dis)) +geom_smooth(aes(y=value,colour=variable),se=FALSE)+
  geom_point(aes(y=value,colour=variable))+guides(colour=FALSE)

#=====================================
#=====Sensitive analysis=====#
n <- vars[j,2]
connect_p <- vars[j,3]
connect_th <- vars[j,3]
dis_p <- vars[j,4]
dis_th <- vars[j,4]
a <- vars[j,5]
th_ratio <- vars[j,6]
phy_sig <- vars[j,7]
rf.res30 <- list()
ttt.2 <- lapply(ttt,func.4)
ttt.2.melt <- melt(ttt.2)
ttt.2.melt$Nspecies <- vars[ttt.2.melt$L1,1]
ttt.2.melt$n <- vars[ttt.2.melt$L1,2]
ttt.2.melt$connect_p <- vars[ttt.2.melt$L1,3]
ttt.2.melt$connect_th <- vars[ttt.2.melt$L1,3]
ttt.2.melt$dis_p <- vars[ttt.2.melt$L1,4]
ttt.2.melt$dis_th <- vars[ttt.2.melt$L1,4]
ttt.2.melt$a <- vars[ttt.2.melt$L1,5]
ttt.2.melt$th_ratio <- vars[ttt.2.melt$L1,6]
ttt.2.melt$phy_sig <- vars[ttt.2.melt$L1,7]
tmp <- ttt.2.melt %>% filter(variable=="phy_dis",net_dis==9);summary(aov(value~Nspecies+n+connect_p+dis_p+a+th_ratio+phy_sig,data=tmp))