rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./_Rfunction/Rfunctions.R")
source("./_Rfunction/MRM.R")
#Packages
#install.packages("devtools")
#install_github("zdk123/SpiecEasi", force=TRUE)
#source("https://bioconductor.org/biocLite.R")
#library(SpiecEasi)
library(vegan)
library(plyr)
library(bmotif)
library(dplyr)
library(ggplot2)
library(ape)
library(latentnet) #Latent Model
library(VBLPCM)
library(igraph)
library(viridis)
library(reshape)
library(picante)
library(randomForest)
library(SpiecEasi)
library(gridExtra)
library(ggpubr)
library(cowplot)
.id <- as.character(unlist(as.data.frame(strsplit(dir("./_phylo.tree"),split=".tree"))))
Levels <- c("SRF_Free","SRF_Part","DCM_Free","DCM_Part","MES_Free","MES_Part","AMD","Topsoil")
Levels.f <- c("SRF_Free","DCM_Free","MES_Free","Topsoil","SRF_Part","DCM_Part","MES_Part","AMD")
phylo.tree <- lapply(paste('./_phylo.tree/',dir("./_phylo.tree"),sep=""),read.tree)
otutable <- lapply(paste('./_otutable/',dir("./_otutable"),sep=""),read.table,header=T,row.names=1,sep="\t")
optima.env <- lapply(paste('./_optima.cal/env/',dir("./_optima.cal/env/"),sep=""),read.table,header=T,row.names=1,sep="\t")
optima.otutable <- lapply(paste('./_optima.cal/otutable/',dir("./_optima.cal/otutable/"),sep=""),read.table,header=T,row.names=1,sep="\t")
optima.id <- as.character(unlist(as.data.frame(strsplit(dir("./_optima.cal/otutable/"),split=".txt"))))
optima.env.null <- lapply(optima.env,function(x){x <- as.matrix(x)
x[x!=1] <- 1;x[is.na(x)] <- 1;return(x)})

otutable <- lapply(otutable,filtering_otu,otu_obs=0.5,Transpose=T)
N <- length(phylo.tree)
phylo.comm <- list()
for(i in 1:N){
  phylo.comm[[i]] <-  match.phylo.comm(phylo.tree[[i]],otutable[[i]])
}
otutable.match <- lapply(phylo.comm,function(x){x$comm})
phylo.match <- lapply(phylo.comm,function(x){x$phy})


Phylo_dist <- lapply(phylo.match,cophenetic)
Phylo_dist <- lapply(Phylo_dist,function(x){(x-min(x))/(max(x)-min(x))})

Optima_dist <- list()
Optima_dist.null <- list()
for(i in 1:length(optima.env))
{
  ef <- colnames(optima.env[[i]])
  Optima_dist[[i]] <- cal.env.optima(t(optima.otutable[[i]]),optima.env[[i]],ef,z_scores = T)
  Optima_dist.null[[i]] <- cal.env.optima(t(optima.otutable[[i]]),optima.env.null[[i]],ef,z_scores = T)
}
names(Optima_dist) <- optima.id
names(Optima_dist.null) <- optima.id

ij <- c(1,2,3,2,3,2,3,4)
Optima_dist2 <- list()
Optima_dist2.null <- list()
for(i in 1:N)
{
  otu.tmp <- colnames(otutable.match[[i]])
  Optima_dist.tmp <- Optima_dist[[ij[i]]]
  Optima_dist.null.tmp <- Optima_dist.null[[ij[i]]]
  Optima_dist2[[i]] <-  Optima_dist.tmp[otu.tmp,]
  Optima_dist2.null[[i]] <-  Optima_dist.null.tmp[otu.tmp,]
}

Optima_dist <- Optima_dist2
Optima_dist.null <- Optima_dist2.null
Optima_dist <- lapply(Optima_dist,function(x){tmp <- as.matrix(vegdist(x,method="euclid"));res <- (tmp-min(tmp))/(max(tmp)-min(tmp));return(res)})
Optima_dist.null <- lapply(Optima_dist.null,function(x){tmp <- as.matrix(vegdist(x,method="euclid"));res <- (tmp-min(tmp))/(max(tmp)-min(tmp));return(res)})

#save(otutable.match,file="./_Rdata/otutable.match.Rdata")
#sparcc.res <- lapply(otutable.match,sparcc,iter=100,inner_iter = 100,th=0.3)
#names(Optima_dist) <- id
#names(Phylo_dist) <- id
#names(sparcc.res) <- id
#save(Optima_dist,file="./_Rdata/Optima_dist.Rdata")
#save(Phylo_dist,file="./_Rdata/Phylo_dist.Rdata")
#save(sparcc.res,file="./_Rdata/sparcc_res.Rdata")

load("./_Rdata/Optima_dist.Rdata")
load("./_Rdata/Phylo_dist.Rdata")
load("./_Rdata/Sparcc_res.Rdata")

#****
N <- length(Phylo_dist)
sparcc.igraph.C <- lapply(sparcc.res,construct_sparcc_network,th=0.65)
sparcc.igraph.Cor <- lapply(sparcc.igraph.C,function(x){x$Cor})
Negative_Cor_propotion <- lapply(sparcc.igraph.Cor,function(x){sum(x <0)/sum(x!=0)})
sparcc.igraph <- lapply(sparcc.igraph.C,function(x){x$g.a})
sparcc.igraph.path_dist <- lapply(sparcc.igraph,function(x){distances(x)})
sparcc.igraph.adj_mat <- lapply(sparcc.igraph,function(x){tmp <- get.adjacency(x);as.matrix(tmp)})
sparcc.igraph.modules <- lapply(sparcc.igraph,cluster_fast_greedy)

Optima_dist.match <- list()
Phylo_dist.match <- list()
for(i in 1:N){
  Phylo_dist.match[[i]] <- match_two_mat(Phylo_dist[[i]],sparcc.igraph.adj_mat[[i]])
  Optima_dist.match[[i]] <- match_two_mat(Optima_dist[[i]],sparcc.igraph.adj_mat[[i]])
}

#Mantel.test
#Mantel.test.res.phylo <- list()
#Mantel.test.res.habitat <- list()
#for(i in 1:N)
#{
#  Mantel.test.res.phylo[[i]] <- mantel.test(Phylo_dist.match[[i]],sparcc.igraph.adj_mat[[i]])
#  Mantel.test.res.habitat[[i]] <- mantel.test(Optima_dist.match[[i]],sparcc.igraph.adj_mat[[i]])
#}
#names(Mantel.test.res.phylo) <- names(Phylo_dist)
#names(Mantel.test.res.habitat) <- names(Optima_dist)
#Mantel.test.res.phylo <- ldply(lapply(Mantel.test.res.phylo,function(x){tmp <- c(x$z.stat,x$p);names(tmp) <- c("z.stat","p");return(tmp)}))
#Mantel.test.res.habitat <- ldply(lapply(Mantel.test.res.habitat,function(x){tmp <- c(x$z.stat,x$p);names(tmp) <- c("z.stat","p");return(tmp)}))
#Mantel.test.res <- list(phylo=Mantel.test.res.phylo,habitat=Mantel.test.res.habitat)
#save(Mantel.test.res,file="./_Rdata/Mantel.test.res.Rdata")
load("./_Rdata/Mantel.test.res.Rdata")
Mantel.test.res <- merge(Mantel.test.res$phylo,Mantel.test.res$habitat,by=".id")[,2:5]
colnames(Mantel.test.res) <- c("phylo_adj.z","phylo_adj.p","habitat_adj.z","habitat_adj.p")
rownames(Mantel.test.res) <- .id
write.table(format(Mantel.test.res,digits=3),"./_table/Mantel.test.HP_vs_adj_mat.txt",sep="\t",col.names = NA,quote = F)
#MRM
#library(ecodist)
#MRM.res <- list()
#for(i in 1:N)
#{
#  MRM.res[[i]] <- MRM(as.dist(sparcc.igraph.adj_mat[[i]])~as.dist(Phylo_dist.match[[i]])+as.dist(Optima_dist.match[[i]]),
#                      nperm=999,method = "logistic")
#}
#names(MRM.res) <- .id
#save(MRM.res,file="./_Rdata/MRM.res.Rdata")
#detach("package:ecodist", unload=TRUE)
load("./_Rdata/MRM.res.Rdata")
MRM.res.coef <- lapply(MRM.res,function(x){x$coef})
MRM.res.coef <- ldply(MRM.res.coef)
MRM.res.coef$Variable <- c("Int","phylo","habitat")
MRM.res.coef <- MRM.res.coef[-which(MRM.res.coef$Variable == "Int"),]
MRM.res.coef$beta <- -MRM.res.coef$b
write.table(format(MRM.res.coef,digits=3),"./_table/MRM.adj.res.coef.txt",quote=F,row.names = F,sep="\t")
MRM.res.coef <- within(MRM.res.coef,beta <- ifelse(pval>=0.1,0,beta))
S <- MRM.res.coef %>% group_by(.id) %>% summarise(S=sum(beta))
MRM.res.coef <- merge(MRM.res.coef,S,by=".id")
MRM.res.coef$'Contribution' <- MRM.res.coef$beta/MRM.res.coef$S
MRM.res.coef <- within(MRM.res.coef,.id <- factor(.id,levels = Levels))
pdf(file="./_figure/p1.mrm.relative.adj_phylo_optima.pdf",width = 5,height = 4)
p1 <- MRM.res.coef %>% ggplot(aes(x=.id,y=Contribution,fill=Variable)) + geom_bar(stat="identity")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+labs(x="",y="Relative Contribution")
print(p1)
dev.off()
pdf(file="./_figure/p2.mrm.Beta.adj_phylo_optima.pdf",width = 5,height = 4)
p2 <- MRM.res.coef %>% ggplot(aes(x=.id,y=-b)) + geom_bar(aes(fill=Variable),stat="identity",position="dodge")+
  labs(y="-Beta",x="")+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+guides(fill=F)
print(p2)
dev.off()
#Randomforest
#Meta data generation
#meta_data.list <- list()
#for(i in 1:N)
#{
#  print(i)
#  meta_data.list[[i]] <- three_dis_mats_2_list(sparcc.igraph.path_dist[[i]],
#                                               Phylo_dist[[i]],Optima_dist[[i]],
#                                               cn=c("node1","node2","net_dis","phylo_dis","habitat_dis"))
#  Nr <- nrow(meta_data.list[[i]])
#  meta_data.list[[i]]$No.module <- "InTwo"
#  meta_data.list[[i]]$In.same_modules <- "No"
#  for(j in 1:Nr)
#  {
#    m1 <- in_which_module(meta_data.list[[i]]$node1[j],sparcc.igraph.modules[[i]])
#    m2 <- in_which_module(meta_data.list[[i]]$node2[j],sparcc.igraph.modules[[i]])
#    if(m1 == m2) {meta_data.list[[i]]$No.module[j] <- as.character(m1);meta_data.list[[i]]$In.same_modules[j] <- "Yes"}
#  }
#}
#names(meta_data.list) <- .id
#save(meta_data.list,file="./_Rdata/meta_data.list.Rdata")
load("./_Rdata/meta_data.list.Rdata")
meta_data.df <- ldply(meta_data.list)

#************
lm.test.for.phy_net_dis <- meta_data.df %>% group_by(.id,No.module) %>% do(data.frame(lm.phy.p=lm.p(.$net_dis,.$phylo_dis)))
lm.test.for.habitat_net_dis <- meta_data.df %>% group_by(.id,No.module) %>% do(data.frame(lm.phy.p=lm.p(.$net_dis,.$habitat_dis)))
lm.test.hp.n <- lm.test.for.phy_net_dis
lm.test.hp.n$lm.habitat.p <- lm.test.for.habitat_net_dis$lm.phy.p

meta_data.df$groups <- paste(meta_data.df$.id,meta_data.df$No.module,sep="_")
lm.test.hp.n$groups <- paste(lm.test.hp.n$.id,lm.test.hp.n$No.module,sep="_")
meta_data2 <- merge.data.frame(meta_data.df,lm.test.hp.n,by.x="groups",by.y="groups",all.x=T)
meta_data2 <- within(meta_data2,phy_sig <- ifelse(lm.phy.p<0.05,"sig","no sig"))
meta_data2 <- within(meta_data2,habitat_sig <- ifelse(lm.habitat.p<0.05,"sig","no sig"))
meta_data.df <- meta_data2
meta_data.df <- within(meta_data.df,phy_sig <- factor(phy_sig,levels = c("sig","no sig")))
meta_data.df <- within(meta_data.df,habitat_sig <- factor(habitat_sig,levels = c("sig","no sig")))
meta_data.df$id <- meta_data.df$.id.x


meta_data.df <- within(meta_data.df,id <- factor(id,levels = Levels.f))
#save(meta_data.df,file="./_Rdata/meta_data.df.Rdata")

#******
pdf(file="./_figure/p3.lm.net_dis_phylo_dis.pdf",width = 6,height = 4)
p3 <- meta_data.df %>%ggplot(aes(x=net_dis))+geom_smooth(aes(y=phylo_dis,colour=In.same_modules,linetype=phy_sig,group=No.module.x),method="lm",se=FALSE,size=0.5)+facet_wrap(~id,nrow=2)+guides(colour=FALSE,linetype=FALSE)+theme_bw()+scale_x_continuous(limits = c(0,12),breaks = seq(0,12,3))+scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))+labs(x="Network path distance",y="Phylogenetic distance")+theme(axis.text = element_text(color="black"),axis.title = element_text(color="black"))
print(p3)
dev.off()
pdf(file="./_figure/p4.lm.net_dis_habitat_dis.pdf",width = 6,height = 4)
p4 <- meta_data.df %>%ggplot(aes(x=net_dis))+geom_smooth(aes(y=habitat_dis,colour=In.same_modules,linetype=habitat_sig,group=No.module.x),method="lm",se=FALSE,size=0.5)+facet_wrap(~id,nrow=2)+guides(colour=FALSE,linetype=FALSE)+theme_bw()+scale_x_continuous(limits = c(0,12),breaks = seq(0,12,3))+scale_y_continuous(limits = c(0,0.75),breaks = seq(0,0.75,0.25))+labs(x="Network path distance",y="Habitat distance")+theme(axis.text = element_text(color="black"),axis.title = element_text(color="black"))
print(p4)
dev.off()
pdf(file="./_figure/p5.smooth.phylo_dis_habitat_dis.pdf",width = 6,height = 4)
p5 <- meta_data.df %>%ggplot(aes(x=phylo_dis,y=habitat_dis))+geom_smooth(se=F)+labs(x="Phylogenetic distance",y="Habitat distance")+theme_bw()+scale_x_continuous(breaks = seq(0,1,0.2))+facet_wrap(~id,nrow=2)
print(p5)
dev.off()
pdf(file="./_figure/p6.smooth.phylo_dis_net_dis.pdf",width = 6,height = 4)
p6 <- meta_data.df %>%ggplot(aes(x=jitter(net_dis),y=phylo_dis))+geom_smooth(se=F)+labs(x="Network distance",y="Phylogenetic distance")+theme_bw()+facet_wrap(~id,nrow=2)+scale_x_continuous(breaks = seq(1,15,3))
print(p6)
dev.off()
pdf(file="./_figure/p7.smooth.habitat_dis_net_dis.pdf",width = 6,height = 4)
p7 <- meta_data.df %>%ggplot(aes(x=jitter(net_dis),y=habitat_dis))+geom_smooth(se=F)+labs(x="Network distance",y="Habitat distance")+theme_bw()+facet_wrap(~id,nrow=2)+scale_x_continuous(breaks = seq(1,15,3))
print(p7)
dev.off()
#rf
#randomforest.res <- lapply(meta_data.list,function(x){randomForest(as.factor(net_dis)~phylo_dis+habitat_dis,data=x,importance=TRUE)})
#randomforest.res.explaination <- lapply(randomforest.res,rf_explanation)
#names(randomforest.res.explaination) <- .id
#randomforest.res.explaination <- melt(randomforest.res.explaination,id.vars = c("net_dis"))

#randomforest.res.explaination <- within(randomforest.res.explaination,L1 <- factor(L1,levels = Levels.f))

#randomforest.res.explaination <- within(randomforest.res.explaination,variable <- factor(variable,levels = c("habitat_dis","phylo_dis")))
#save(randomforest.res.explaination,file="./_Rdata/randomforest.res.explaination.Rdata")
load("./_Rdata/randomforest.res.explaination.Rdata")
randomforest.res.explaination %>% filter( L1 != "AMD" & L1 != "Topsoil") %>% group_by(variable) %>% do(data.frame(max=max(.$value)))

pdf(file="./_figure/p8.randomforest.HP_explaination.pdf",width = 6,height = 4)
p8 <- randomforest.res.explaination %>% ggplot(aes(x=as.numeric(net_dis),y=value)) + geom_smooth(aes(colour=variable),se=FALSE)+
  theme_bw()+facet_wrap(~L1,nrow=2,scales = "free_x")+guides(colour=FALSE)+geom_point(aes(colour=variable))+scale_x_continuous(limits = c(0,15),breaks = seq(0,15,3))+labs(x="Network distance",y="Explaination")
print(p8)
dev.off()
##Part2
#***MNPD.ses,MNHD.ses***
#sparcc.igraph <- lapply(sparcc.igraph.C,function(x){x$g.with_d0})
#sparcc.igraph <- lapply(sparcc.igraph.C,function(x){x$g.a})

#mnpd.ses.list <- list()
#mnhd.ses.list <- list()
#node.degre.list <- list()
#for(i in 1:N)
#{
#  mnpd.ses.list[[i]] <- mnxd.ses(sparcc.igraph[[i]],Phylo_dist[[i]],delete.aloned_nodes=T,cpu=2)
#  mnhd.ses.list[[i]] <- mnxd.ses(sparcc.igraph[[i]],Optima_dist[[i]],delete.aloned_nodes=T,cpu=2)
#}
#node.degre.list <- lapply(sparcc.igraph,function(x){as.matrix(degree(x),ncol=1)})
#names(mnpd.ses.list) <- .id
#names(mnhd.ses.list) <- .id
#names(node.degre.list) <- .id

#save(mnpd.ses.list,mnhd.ses.list,node.degre.list,file="./_Rdata/node_mnxd_degree.Rdata")
load("./_Rdata/node_mnxd_degree.Rdata")
#Node_mnpd_mnhd_degre <- list()
#for(i in 1:N)
#{
#  tmp <- merge(mnpd.ses.list[[i]],mnhd.ses.list[[i]],by=0)
#  rownames(tmp) <- tmp$Row.names
#  tmp <- tmp[,-1]
#  tmp <- merge(tmp,node.degre.list[[i]],by=0)
#  rownames(tmp) <- tmp$Row.names
#  tmp <- tmp[,-1]
#  colnames(tmp) <- c("mnpd.ses","mnhd.ses","degree")
#  Node_mnpd_mnhd_degre[[i]] <- tmp
#}
#names(Node_mnpd_mnhd_degre) <- .id
#save(Node_mnpd_mnhd_degre,file="./_Rdata/Node_mnpd_mnhd_degre.Rdata")
load("./_Rdata/Node_mnpd_mnhd_degre.Rdata")
#Node_mnpd_mnhd_degre.df <- ldply(Node_mnpd_mnhd_degre)

#Node_mnpd_mnhd_degre.df <- Node_mnpd_mnhd_degre.df %>% group_by(.id) %>% do(data.frame(density=get_density(.$mnpd.ses,.$mnhd.ses),degree=.$degree,mnpd.ses=.$mnpd.ses,mnhd.ses=.$mnhd.ses))
#Node_mnpd_mnhd_degre.df <- within(Node_mnpd_mnhd_degre.df,.id <- factor(.id,levels = Levels.f))                 

#save(Node_mnpd_mnhd_degre.df,file="./_Rdata/Node_mnpd_mnhd_degre.df.Rdata")
load("./_Rdata/Node_mnpd_mnhd_degre.df.Rdata")
n1=Node_mnpd_mnhd_degre.df %>% group_by(.id) %>% summarise(n1=count(.id))
n2=Node_mnpd_mnhd_degre.df %>% filter(mnpd.ses < -2)  %>% group_by(.id) %>% summarise(n1=count(.id))
n3=Node_mnpd_mnhd_degre.df %>% filter(mnpd.ses < 0)  %>% group_by(.id) %>% summarise(n1=count(.id))
n4=Node_mnpd_mnhd_degre.df %>% filter(mnhd.ses < -2)  %>% group_by(.id) %>% summarise(n1=count(.id))
n5=Node_mnpd_mnhd_degre.df %>% filter(mnhd.ses < 0)  %>% group_by(.id) %>% summarise(n1=count(.id))
n6=Node_mnpd_mnhd_degre.df %>% filter(mnhd.ses < 0 | mnpd.ses <0)  %>% group_by(.id) %>% summarise(n1=count(.id))
n7=Node_mnpd_mnhd_degre.df %>% filter(mnhd.ses < -2 | mnpd.ses < -2)  %>% group_by(.id) %>% summarise(n1=count(.id))
n8=Node_mnpd_mnhd_degre.df %>% filter(mnhd.ses < 0 & mnpd.ses <0)  %>% group_by(.id) %>% summarise(n1=count(.id))
n9=Node_mnpd_mnhd_degre.df %>% filter(mnhd.ses < -2 & mnpd.ses < -2)  %>% group_by(.id) %>% summarise(n1=count(.id))


n1$mnpd_less_n2=n2$n1$freq/n1$n1$freq
n1$mnpd_less_0=n3$n1$freq/n1$n1$freq
n1$mnhd_less_n2=n4$n1$freq/n1$n1$freq
n1$mnhd_less_0=n5$n1$freq/n1$n1$freq
n1$either_less_0=n6$n1$freq/n1$n1$freq
n1$either_less_n2=n7$n1$freq/n1$n1$freq
n1$both_less_0=n8$n1$freq/n1$n1$freq
n1$both_less_n2=n9$n1$freq/n1$n1$freq

pdf(file="./_figure/p9.node.mnhd_mnpd.ses.pdf",width = 8,height = 4.8)
p9 <- Node_mnpd_mnhd_degre.df %>% ggplot(aes(x=mnpd.ses,y=mnhd.ses,colour=density))+geom_point(size=0.2)+
  geom_hline(yintercept=c(-2,2))+geom_vline(xintercept=c(-2,2))+scale_color_viridis()+facet_wrap(~.id,nrow=2,scales = "free")+scale_x_continuous(expand = c(0.1,0.1,0.1,0.1))+theme_bw()+labs(y="MNHD.ses",x="MNPD.ses",colour="Density")
print(p9)
dev.off()
pdf(file="./_figure/p10.node.degree_mnpd.ses.pdf",width = 6,height = 4)
p10 <- Node_mnpd_mnhd_degre.df %>% ggplot(aes(x=mnpd.ses,y=degree))+geom_point(size=0.2)+geom_vline(xintercept = c(-2,2),color="red",size=0.2)+facet_wrap(~.id,nrow=2,scales = "free")+labs(x="MNPD.ses",y="Degree")+theme_bw()+scale_x_continuous(expand = c(0.1,0.1,0.1,0.1))
print(p10)
dev.off()
pdf(file="./_figure/p11.node.degree_mnhd.ses.pdf",width = 6,height = 4)
p11 <- Node_mnpd_mnhd_degre.df %>% ggplot(aes(x=mnhd.ses,y=degree))+geom_point(size=0.2)+geom_vline(xintercept = c(-2,2),color="red",size=0.2)+facet_wrap(~.id,nrow=2,scales = "free")+labs(x="MNHD.ses",y="Degree")+theme_bw()+scale_x_continuous(expand = c(0.1,0.1,0.1,0.1))
print(p11)
dev.off()

Node_mnpd_mnhd_degre.glm <- lapply(Node_mnpd_mnhd_degre,function(x){
    glm.temp <- summary(glm(degree~mnpd.ses+mnhd.ses,data=x))
    return(glm.temp$coefficients)
}) 
Node_mnpd_mnhd_degre.glm.df <- ldply(Node_mnpd_mnhd_degre.glm)
Node_mnpd_mnhd_degre.glm.df$Para <- c("Intercept","MNPD.ses","MNHD.ses")
Node_mnpd_mnhd_degre.glm.df <- Node_mnpd_mnhd_degre.glm.df %>% select(.id,Para,Estimate,`Std. Error`,`t value`,`Pr(>|t|)`)
write.table(format(Node_mnpd_mnhd_degre.glm.df,digits=3),"./_table/glm.degree~MNPD+MNHD.df.txt",quote=F,row.names = F,sep="\t",qmethod="double")
Node_mnpd_mnhd_degre.rf <- lapply(Node_mnpd_mnhd_degre,function(x){randomForest(degree~mnpd.ses+mnhd.ses,
                                                                                data=x,importance=TRUE)})

#***Motif & Position***
sparcc.igraph.adj_mat <- lapply(sparcc.igraph.adj_mat,as.matrix)
#motif.list <- lapply(sparcc.igraph.adj_mat,mcount,six_node = TRUE,normalise = TRUE)
#position.list <- lapply(sparcc.igraph.adj_mat,function(x){
#  x <- as.matrix(x)
#  tmp <- positions(x,six_node = TRUE, level = "rows", normalisation ="none")
#  tmp.cs <- colSums(tmp)
#  tmp.freq <- t(apply(tmp,1,function(x,y){x/y},y=tmp.cs))
#  return(tmp.freq)
#  }
#)
#names(motif.list) <- .id
#names(position.list) <- .id
#save(motif.list,position.list,file="./_Rdata/Motif_position.Rdata")
load("./_Rdata/Motif_position.Rdata")
motif.df <- ldply(motif.list)
motif.df <- within(motif.df,.id<-factor(.id,levels = Levels))
pdf(file="./_figure/p12.motif.frequency.pdf",width = 6,height = 4)
p12 <- motif.df %>% ggplot(aes(x=motif))+geom_line(aes(y=normalise_sum,colour=.id),size=0.3)+theme_bw()+labs(x="Motif",y="Normalise frequency",colour="")
print(p12)
dev.off()

motif.matrix <- motif_table(motif.list)
motif.nmds <- as.data.frame(metaMDS(motif.matrix)$points)
motif.nmds$.id<- rownames(motif.nmds)
pdf(file="./_figure/p13.motif.NMDS.pdf",width = 6,height = 4)
p13 <-  motif.nmds %>% ggplot(aes(x=MDS1,y=MDS2,colour=.id))+geom_point()+theme_bw()+labs(colour="")
print(p13)
dev.off()

#motif.matrix.tara <- motif.matrix[2:7,]
#tara.g <- rep(c("Free","Part"),3)
#mrpp(motif.matrix.tara,tara.g,distance = "horn")

Position_mnxd.list <- list()
for(i in 1:N)
{
  Position_mnxd.list[[i]] <- position.nmxd(position.list[[i]],Node_mnpd_mnhd_degre[[i]])
}
names(Position_mnxd.list) <- .id
Position_mnxd.list.df <- ldply(Position_mnxd.list)
Position_mnxd.list.df <- within(Position_mnxd.list.df,.id <- factor(.id,levels = Levels.f))                 

pdf(file="./_figure/p14.position.mnxd.pdf",width = 6,height = 4)
p14 <- Position_mnxd.list.df%>% ggplot(aes(x=mnpd,y=mnhd,colour=Position))+geom_point(size=0.2)+
  geom_hline(yintercept=c(-2))+geom_vline(xintercept=c(-2))+scale_color_viridis()+geom_abline(intercept = 0, slope = 1,color="red",linetype="dashed")+facet_wrap(~.id,nrow=2)+labs(x="MNPD.ses",y="MNHD.ses")+theme_bw()
print(p14)
dev.off()

position.list <- lapply(position.list,pos_pro)
Position_dist.list <- lapply(position.list,function(x){tmp <- as.matrix(vegdist(x,method="bray"))
tmp <- (tmp - min(tmp))/(max(tmp)-min(tmp));return(tmp)})

Mantel.test.pos_phylo <- list()
Mantel.test.pos_habitat <- list()
for(i in 1:N)
{
  Mantel.test.pos_phylo[[i]] <- mantel.test(Phylo_dist.match[[i]],Position_dist.list[[i]])
  Mantel.test.pos_habitat[[i]] <- mantel.test(Optima_dist.match[[i]],Position_dist.list[[i]])
}
names(Mantel.test.pos_phylo) <- names(Node_mnpd_mnhd_degre)
names(Mantel.test.pos_habitat) <- names(Node_mnpd_mnhd_degre)
Mantel.test.pos_phylo.res <- merge(ldply(lapply(Mantel.test.pos_phylo,function(x){x$z.stat})),ldply(lapply(Mantel.test.pos_phylo,function(x){x$p})),by=".id")
Mantel.test.pos_habitat.res <- merge(ldply(lapply(Mantel.test.pos_habitat,function(x){x$z.stat})),ldply(lapply(Mantel.test.pos_habitat,function(x){x$p})),by=".id")
colnames(Mantel.test.pos_phylo.res) <- c("id","phylo_pos.z","phylo_position.p")
colnames(Mantel.test.pos_habitat.res) <- c("id","habitat_pos.z","phylo_position.p")
Mantel.test.pos_hp <- cbind(Mantel.test.pos_phylo.res,Mantel.test.pos_phylo.res[,-1])
write.table(format(Mantel.test.pos_hp,digits=3),"./_table/Mantel.test.pos_hp.txt",quote = F,sep="\t",row.names = F)

library(ecodist)
MRM.pos_hp <- list()
for(i in 1:N)
{
  MRM.pos_hp[[i]] <- MRM(as.dist(Position_dist.list[[i]])~as.dist(Phylo_dist.match[[i]])+as.dist(Optima_dist.match[[i]]),nperm = 999)
}
detach("package:ecodist", unload=TRUE)
names(MRM.pos_hp) <- .id
MRM.pos.res.coef <- lapply(MRM.pos_hp,function(x){x$coef})
MRM.pos.res.coef <- ldply(MRM.pos.res.coef)
MRM.pos.res.coef$Variable <- c("Int","phylo","habitat")
MRM.pos.res.coef <- MRM.pos.res.coef[-which(MRM.pos.res.coef$Variable == "Int"),]


colnames(MRM.pos.res.coef) <- c(".id","Estimate","pval","Variable")
write.table(format(MRM.pos.res.coef,digits=3),"./_table/MRM.pos.res.coef.txt",row.names = F,quote=F,sep="\t")
MRM.pos.res.coef <- within(MRM.pos.res.coef,Estimate <- ifelse(pval>=0.1,0,Estimate))
S <- MRM.pos.res.coef %>% group_by(.id) %>% do(data.frame(S=sum(abs(.$Estimate))))
MRM.pos.res.coef <- merge(MRM.pos.res.coef,S,by=".id")
MRM.pos.res.coef$'Contribution' <- MRM.pos.res.coef$Estimate/MRM.pos.res.coef$S
MRM.pos.res.coef <- within(MRM.pos.res.coef,.id <- factor(.id,levels = Levels))



pdf(file="./_figure/p30.mrm.relative.pos_phylo_optima.pdf",width = 5,height = 4)
p30 <- MRM.pos.res.coef %>% ggplot(aes(x=.id,y=Contribution,fill=Variable)) + geom_bar(stat="identity")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+labs(x="",y="Relative Contribution")
print(p30)
dev.off()
pdf(file="./_figure/p31.mrm.Estimate.pos_phylo_optima.pdf",width = 5,height = 4)
p31 <- MRM.pos.res.coef %>% ggplot(aes(x=.id,y=Estimate)) + geom_bar(aes(fill=Variable),stat="identity",position="dodge")+
  labs(y="Estimate",x="")+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+guides(fill=F)
print(p31)
dev.off()


#Latent space
sparcc.network <- lapply(sparcc.igraph.adj_mat,network,directed=F)
#sparcc.ls.res <- list()
#for(i in 1:N)
#{
#  sparcc.ls.res[[i]] <- BICg_VB(sparcc.network[[i]],d=2,maxg = (max(sparcc.igraph.modules[[i]]$membership)+5))
#}
#names(sparcc.ls.res) <- .id
#save(sparcc.ls.res,file="./_Rdata/sparcc.ls.res.Rdata")
load("./_Rdata/sparcc.ls.res.Rdata")
sparcc.ls.bic <- ldply(sparcc.ls.res) %>%  within({.id <- factor(.id,levels = Levels.f)})
write.table(format(sparcc.ls.bic,digits = 3),"./_table/sparcc.ls.bic.txt",row.names = F,quote = F,sep="\t")
pdf("./_figure/p15.ls.BIC.pdf",width = 6,height = 4)
p15 <- sparcc.ls.bic %>%
    ggplot(aes(x=No.Groups,y=`Overall BIC`)) + geom_point(size=1) + geom_line()+scale_x_continuous(limits = c(0,30),breaks = seq(0,30,5))+facet_wrap(~.id,nrow=2,scales="free")+theme_bw()
print(p15)
dev.off()

bestG <- lapply(sparcc.ls.res,function(x){x$No.Groups[which(x$`Overall BIC`==min(x$`Overall BIC`))]})
bestG <- ldply(bestG)
colnames(bestG) <- c("id","Best.G")
write.table(bestG,"./_table/bestG.txt",row.names = F,quote = F,sep="\t")

#sparcc.lsp.mcmc.lcc <- list()
#sparcc.lsp.mcmc.lcc.igm <- list()
#for( i in 1:N)
#{
#  tmp <- bestG[i,2]
#  sparcc.lsp.mcmc.lcc[[i]] <- vblpcmstart(sparcc.network[[i]],d=2,G=tmp)
#  tmp <- 1
#  sparcc.lsp.mcmc.lcc.igm[[i]] <- vblpcmstart(sparcc.network[[i]],d=2,G=tmp)
#}
#names(sparcc.lsp.mcmc.lcc) <- .id
#names(sparcc.lsp.mcmc.lcc.igm) <- .id

#save(sparcc.lsp.mcmc.lcc,sparcc.lsp.mcmc.lcc.igm,file="./_Rdata/sparcc.lsp.mcmc.lcc.res.Rdata")
load("./_Rdata/sparcc.lsp.mcmc.lcc.res.Rdata")

sparcc.lsp.mcmc.lcc.Z <- lapply(sparcc.lsp.mcmc.lcc,ls_z)
sparcc.lsp.mcmc.lcc.igm.Z <- lapply(sparcc.lsp.mcmc.lcc.igm,ls_z)

sparcc.lsp.mcmc.lcc.Z.dist <- lapply(sparcc.lsp.mcmc.lcc.Z,Z_dist)
sparcc.lsp.mcmc.lcc.igm.Z.dist <- lapply(sparcc.lsp.mcmc.lcc.igm.Z,Z_dist)

sparcc.lsp.mcmc.lcc.Z.dist.v <- lapply(sparcc.lsp.mcmc.lcc.Z.dist,function(x){three_dis_mats_2_list(x[[1]],x[[2]],x[[3]],cn=c("n1","n2","Z_eucli.dis","z1_manha.dis","z2_manha.dis"))})
sparcc.lsp.mcmc.lcc.igm.Z.dist.v <- lapply(sparcc.lsp.mcmc.lcc.igm.Z.dist,function(x){three_dis_mats_2_list(x[[1]],x[[2]],x[[3]],cn=c("n1","n2","Z_eucli.dis.igm","z1_manha.dis.igm","z2_manha.dis.igm"))})

sparcc.lsp.mcmc.lcc.Z.dist.meta <- list()
for(i in 1:N)
{
  tmp.x <- sparcc.lsp.mcmc.lcc.Z.dist.v[[i]]
  tmp.y <- sparcc.lsp.mcmc.lcc.igm.Z.dist.v[[i]]
  rownames(tmp.x) <- paste(tmp.x$n1,tmp.x$n2)
  rownames(tmp.y) <- paste(tmp.y$n1,tmp.y$n2)
  tmp.y <- tmp.y[,-c(1,2)]
  tmp.merge <- merge(tmp.x,tmp.y,by=0)[,-1]
  sparcc.lsp.mcmc.lcc.Z.dist.meta[[i]] <- tmp.merge
}
sparcc.lsp.mcmc.lcc.Z.dist.meta <- lapply(sparcc.lsp.mcmc.lcc.Z.dist.meta,z_scores_4_c,vs=seq(3,8))
net_phylo_habitat_dis <- list()
for(i in 1:N)
{
  net_phylo_habitat_dis[[i]] <- three_dis_mats_2_list(sparcc.igraph.path_dist[[i]],Phylo_dist.match[[i]],Optima_dist.match[[i]],cn=c("n1","n2","net_dis","phylo_dis","habitat_dis"))
}
sparcc.lsp.mcmc.lcc.meta <- list()
for( i in 1:N)
{
  tmp.x <- sparcc.lsp.mcmc.lcc.Z.dist.meta[[i]]
  tmp.y <- net_phylo_habitat_dis[[i]]
  rownames(tmp.x) <- paste(tmp.x$n1,tmp.x$n2)
  rownames(tmp.y) <- paste(tmp.y$n1,tmp.y$n2)
  tmp.y <- tmp.y[,-c(1,2)]
  tmp.merge <- merge(tmp.x,tmp.y,by=0)[,-1]
  sparcc.lsp.mcmc.lcc.meta[[i]] <- tmp.merge
}

names(sparcc.lsp.mcmc.lcc.meta) <- .id
sparcc.lsp.mcmc.lcc.meta.df <- ldply(sparcc.lsp.mcmc.lcc.meta)
sparcc.lsp.mcmc.lcc.meta.df  <- within(sparcc.lsp.mcmc.lcc.meta.df,.id <- factor(.id,levels = Levels.f))     

pdf("./_figure/p16.ls.Z_net_dis.pdf",width = 6,height = 4)
p16 <- ggplot(data=sparcc.lsp.mcmc.lcc.meta.df,aes(x=jitter(net_dis),y=Z_eucli.dis))+geom_point()+geom_smooth(method="lm")+labs(x="Network path distance",y="Z distance")+facet_wrap(~.id,nrow=2)+theme_bw()+theme(axis.text = element_text(color="black"))
print(p16)
dev.off()

pdf("./_figure/p17.ls.Z1_net_dis.pdf",width = 6,height = 4)
p17 <- ggplot(data=sparcc.lsp.mcmc.lcc.meta.df,aes(x=jitter(net_dis),y=z1_manha.dis))+geom_point()+geom_smooth(method="lm")+labs(x="Network path distance",y="Z1 distance")+facet_wrap(~.id,nrow=2)+theme_bw()+theme(axis.text = element_text(color="black"))
print(p17)
dev.off()

pdf("./_figure/p18.ls.Z2_net_dis.pdf",width = 6,height = 4)
p18 <- ggplot(data=sparcc.lsp.mcmc.lcc.meta.df,aes(x=jitter(net_dis),y=z2_manha.dis))+geom_point()+geom_smooth(method="lm")+labs(x="Network path distance",y="Z2 distance")+facet_wrap(~.id,nrow=2)+theme_bw()+theme(axis.text = element_text(color="black"))
print(p18)
dev.off()


pdf("./_figure/p19.ls.igm.Z_net_dis.pdf",width = 6,height = 4)
p19 <- ggplot(data=sparcc.lsp.mcmc.lcc.meta.df,aes(x=jitter(net_dis),y=Z_eucli.dis.igm))+geom_point()+geom_smooth(method="lm")+labs(x="Network path distance",y="Z distance")+facet_wrap(~.id,nrow=2)+theme_bw()+theme(axis.text = element_text(color="black"))
print(p19)
dev.off()

pdf("./_figure/p20.ls.igm.Z1_net_dis.pdf",width = 6,height = 4)
p21 <- ggplot(data=sparcc.lsp.mcmc.lcc.meta.df,aes(x=jitter(net_dis),y=z1_manha.dis.igm))+geom_point()+geom_smooth(method="lm")+labs(x="Network path distance",y="Z1 distance")+facet_wrap(~.id,nrow=2)+theme_bw()+theme(axis.text = element_text(color="black"))
print(p21)
dev.off()

pdf("./_figure/p21.ls.igm.Z2_net_dis.pdf",width = 6,height = 4)
p21 <- ggplot(data=sparcc.lsp.mcmc.lcc.meta.df,aes(x=jitter(net_dis),y=z2_manha.dis.igm))+geom_point()+geom_smooth(method="lm")+labs(x="Network path distance",y="Z2 distance")+facet_wrap(~.id,nrow=2)+theme_bw()+theme(axis.text = element_text(color="black"))
print(p21)
dev.off()

sparcc.lsp.pearson.cor <- sparcc.lsp.mcmc.lcc.meta.df %>% group_by(.id,net_dis) %>%
  do(data.frame(Cor.Z.phylo=cor(.$Z_eucli.dis,.$phylo_dis),
                Cor.Z.phylo.p=cor.test.p(.$Z_eucli.dis,.$phylo_dis),
                Cor.Z.habitat=cor(.$Z_eucli.dis,.$habitat_dis),
                Cor.Z.habitat.p=cor.test.p(.$Z_eucli.dis,.$habitat_dis),
                Cor.z1.phylo=cor(.$z1_manha.dis,.$phylo_dis),
                Cor.z1.phylo.p=cor.test.p(.$z1_manha.dis,.$phylo_dis),
                Cor.z1.habitat=cor(.$z1_manha.dis,.$habitat_dis),
                Cor.z1.habitat.p=cor.test.p(.$z1_manha.dis,.$habitat_dis),
                Cor.z2.phylo=cor(.$z2_manha.dis,.$phylo_dis),
                Cor.z2.phylo.p=cor.test.p(.$z2_manha.dis,.$phylo_dis),
                Cor.z2.habitat=cor(.$z2_manha.dis,.$habitat_dis),
                Cor.z2.habitat.p=cor.test.p(.$z2_manha.dis,.$habitat_dis),
                #igm
                Cor.igm.Z.phylo=cor(.$Z_eucli.dis.igm,.$phylo_dis),
                Cor.igm.Z.phylo.p=cor.test.p(.$Z_eucli.dis.igm,.$phylo_dis),
                Cor.igm.Z.habitat=cor(.$Z_eucli.dis.igm,.$habitat_dis),
                Cor.igm.Z.habitat.p=cor.test.p(.$Z_eucli.dis.igm,.$habitat_dis),
                Cor.igm.z1.phylo=cor(.$z1_manha.dis.igm,.$phylo_dis),
                Cor.igm.z1.phylo.p=cor.test.p(.$z1_manha.dis.igm,.$phylo_dis),
                Cor.igm.z1.habitat=cor(.$z1_manha.dis.igm,.$habitat_dis),
                Cor.igm.z1.habitat.p=cor.test.p(.$z1_manha.dis.igm,.$habitat_dis),
                Cor.igm.z2.phylo=cor(.$z2_manha.dis.igm,.$phylo_dis),
                Cor.igm.z2.phylo.p=cor.test.p(.$z2_manha.dis.igm,.$phylo_dis),
                Cor.igm.z2.habitat=cor(.$z2_manha.dis.igm,.$habitat_dis),
                Cor.igm.z2.habitat.p=cor.test.p(.$z2_manha.dis.igm,.$habitat_dis)
  ))

sparcc.lsp.pearson.cor <- as.data.frame(sparcc.lsp.pearson.cor)

sparcc.lsp.pearson.cor <- sparcc.lsp.pearson.cor %>% filter(net_dis <=6) %>% 
  within(
    Cor.Z.phylo <- ifelse(Cor.Z.phylo.p>0.05,0,Cor.Z.phylo),
    Cor.Z.habitat <- ifelse(Cor.Z.habitat.p>0.05,0,Cor.Z.habitat),
    Cor.z1.phylo <- ifelse(Cor.z1.phylo.p>0.05,0,Cor.z1.phylo),
    Cor.z1.habitat <- ifelse(Cor.z1.habitat.p>0.05,0,Cor.z1.habitat),
    Cor.z2.phylo <- ifelse(Cor.z2.phylo.p>0.05,0,Cor.z2.phylo),
    Cor.z2.habitat <- ifelse(Cor.z2.habitat.p>0.05,0,Cor.z2.habitat),
    Cor.igm.Z.phylo <- ifelse(Cor.igm.Z.phylo.p>0.05,0,Cor.igm.Z.phylo),
    Cor.igm.Z.habitat <- ifelse(Cor.igm.Z.habitat.p>0.05,0,Cor.igm.Z.habitat),
    Cor.igm.z1.phylo <- ifelse(Cor.igm.z1.phylo.p>0.05,0,Cor.igm.z1.phylo),
    Cor.igm.z1.habitat <- ifelse(Cor.igm.z1.habitat.p>0.05,0,Cor.igm.z1.habitat),
    Cor.igm.z2.phylo <- ifelse(Cor.igm.z2.phylo.p>0.05,0,Cor.igm.z2.phylo),
    Cor.igm.z2.habitat <- ifelse(Cor.igm.z2.habitat.p>0.05,0,Cor.igm.z2.habitat)
)

write.table(format(sparcc.lsp.pearson.cor,digits=3),"./_table/sparcc.lsp.pearson.cor.txt",row.names = F,quote = F,sep="\t")
#save(sparcc.lsp.pearson.cor,file="./_Rdata/sparcc.lsp.pearson.cor.Rdata")
load("./_Rdata/sparcc.lsp.pearson.cor.Rdata")
pdf("./_figure/p22.lsp.Z.pearson.cor.pdf",width = 6,height = 4)
p22 <- sparcc.lsp.pearson.cor %>% select(.id,net_dis,Cor.Z.phylo,Cor.Z.habitat) %>%
  melt.data.frame(measure.vars = c("Cor.Z.phylo","Cor.Z.habitat")) %>% 
  within({variable <- factor(variable,levels=c("Cor.Z.habitat","Cor.Z.phylo"));.id <- factor(.id,levels = Levels.f)}) %>%
  ggplot(aes(x=net_dis,y=value)) + geom_bar(aes(group=variable,fill=variable),stat="identity",position="dodge")+labs(x="Network path distance",y="Pearson' r")+facet_wrap(~.id,nrow=2)+guides(fill=F)+theme_bw()+scale_x_continuous(breaks = seq(1,6))
print(p22)
dev.off()

pdf("./_figure/p23.lsp.Z1.pearson.cor.pdf",width = 6,height = 4)
p23 <- sparcc.lsp.pearson.cor %>% select(.id,net_dis,Cor.z1.phylo,Cor.z1.habitat) %>%
  melt.data.frame(measure.vars = c("Cor.z1.phylo","Cor.z1.habitat")) %>% 
  within({variable <- factor(variable,levels=c("Cor.z1.habitat","Cor.z1.phylo"));.id <- factor(.id,levels = Levels.f)}) %>%
  ggplot(aes(x=net_dis,y=value)) + geom_bar(aes(group=variable,fill=variable),stat="identity",position="dodge")+labs(x="Network path distance",y="Pearson' r")+facet_wrap(~.id,nrow=2)+guides(fill=F)+theme_bw()+scale_x_continuous(breaks = seq(1,6))
print(p23)
dev.off()

pdf("./_figure/p24.lsp.Z2.pearson.cor.pdf",width = 6,height = 4)
p24 <- sparcc.lsp.pearson.cor %>% select(.id,net_dis,Cor.z2.phylo,Cor.z2.habitat) %>%
  melt.data.frame(measure.vars = c("Cor.z2.phylo","Cor.z2.habitat")) %>% 
  within({variable <- factor(variable,levels=c("Cor.z2.habitat","Cor.z2.phylo"));.id <- factor(.id,levels = Levels.f)}) %>%
  ggplot(aes(x=net_dis,y=value)) + geom_bar(aes(group=variable,fill=variable),stat="identity",position="dodge")+labs(x="Network path distance",y="Pearson' r")+facet_wrap(~.id,nrow=2)+guides(fill=F)+theme_bw()+scale_x_continuous(breaks = seq(1,6))
print(p24)
dev.off()

pdf("./_figure/p25.lsp.igm.Z.pearson.cor.pdf",width = 6,height = 4)
p25 <- sparcc.lsp.pearson.cor %>% select(.id,net_dis,Cor.igm.Z.phylo,Cor.igm.Z.habitat) %>%
  melt.data.frame(measure.vars = c("Cor.igm.Z.phylo","Cor.igm.Z.habitat")) %>% 
  within({variable <- factor(variable,levels=c("Cor.igm.Z.habitat","Cor.igm.Z.phylo"));.id <- factor(.id,levels = Levels.f)}) %>%
  ggplot(aes(x=net_dis,y=value)) + geom_bar(aes(group=variable,fill=variable),stat="identity",position="dodge")+labs(x="Network path distance",y="Pearson' r")+facet_wrap(~.id,nrow=2)+guides(fill=F)+theme_bw()+scale_x_continuous(breaks = seq(1,6))
print(p25)
dev.off()

pdf("./_figure/p26.lsp.igm.Z1.pearson.cor.pdf",width = 6,height = 4)
p26 <- sparcc.lsp.pearson.cor %>% select(.id,net_dis,Cor.igm.z1.phylo,Cor.igm.z1.habitat) %>%
  melt.data.frame(measure.vars = c("Cor.igm.z1.phylo","Cor.igm.z1.habitat")) %>% 
  within({variable <- factor(variable,levels=c("Cor.igm.z1.habitat","Cor.igm.z1.phylo"));.id <- factor(.id,levels = Levels.f)}) %>%
  ggplot(aes(x=net_dis,y=value)) + geom_bar(aes(group=variable,fill=variable),stat="identity",position="dodge")+labs(x="Network path distance",y="Pearson' r")+facet_wrap(~.id,nrow=2)+guides(fill=F)+theme_bw()+scale_x_continuous(breaks = seq(1,6))
print(p26)
dev.off()

pdf("./_figure/p27.lsp.igm.Z2.pearson.cor.pdf",width = 6,height = 4)
p27 <- sparcc.lsp.pearson.cor %>% select(.id,net_dis,Cor.igm.z2.phylo,Cor.igm.z2.habitat) %>%
  melt.data.frame(measure.vars = c("Cor.igm.z2.phylo","Cor.igm.z2.habitat")) %>% 
  within({variable <- factor(variable,levels=c("Cor.igm.z2.habitat","Cor.igm.z2.phylo"));.id <- factor(.id,levels = Levels.f)}) %>%
  ggplot(aes(x=net_dis,y=value)) + geom_bar(aes(group=variable,fill=variable),stat="identity",position="dodge")+labs(x="Network path distance",y="Pearson' r")+facet_wrap(~.id,nrow=2)+guides(fill=F)+theme_bw()+scale_x_continuous(breaks = seq(1,6))
print(p27)
dev.off()

sparcc.lsp.G <- lapply(sparcc.lsp.mcmc.lcc,ls_G)
sparcc.lsp.G.dist <- lapply(sparcc.lsp.G,ls_Z_dist,Met="bray")
mrm.G.hp.dist <- list()
mantel.G.phylo.dist <- list()
mantel.G.habitat.dist <- list()
library(ecodist)
for(i in 1:N)
{
  tmp.G <- sparcc.lsp.G.dist[[i]]
  tmp.p <- Phylo_dist.match[[i]]
  tmp.h <- Optima_dist.match[[i]]
  tmp.p.match <- match_two_mat(tmp.p,tmp.G)
  tmp.h.match <- match_two_mat(tmp.h,tmp.G)
  mrm.G.hp.dist[[i]] <- MRM(as.dist(tmp.G)~as.dist(tmp.p.match)+as.dist(tmp.h.match),nperm=999)
}
detach("package:ecodist", unload=TRUE)
names(mrm.G.hp.dist) <- .id
for(i in 1:N)
{
  tmp.G <- sparcc.lsp.G.dist[[i]]
  tmp.p <- Phylo_dist.match[[i]]
  tmp.h <- Optima_dist.match[[i]]
  tmp.p.match <- match_two_mat(tmp.p,tmp.G)
  tmp.h.match <- match_two_mat(tmp.h,tmp.G)
  mantel.G.phylo.dist[[i]] <- mantel.test(tmp.G,tmp.p.match)
  mantel.G.habitat.dist[[i]] <- mantel.test(tmp.G,tmp.h.match)
}
names(mantel.G.phylo.dist) <- .id
names(mantel.G.habitat.dist) <- .id
Mantel.test.G.phylo <- ldply(lapply(mantel.G.phylo.dist,function(x){tmp <- c(x$z.stat,x$p);names(tmp) <- c("z.stat","p");return(tmp)}))
Mantel.test.G.habitat <- ldply(lapply(mantel.G.habitat.dist,function(x){tmp <- c(x$z.stat,x$p);names(tmp) <- c("z.stat","p");return(tmp)}))
colnames(Mantel.test.G.phylo) <- c("id","phylo.z","phylo.p")
colnames(Mantel.test.G.habitat) <- c("id","habitat.z","habitat.p")
Mantel.test.G <- merge(Mantel.test.G.phylo,Mantel.test.G.habitat,by="id")
write.table(format(Mantel.test.G,digits=3),"./_table/Mantel.test.G.txt",quote=F,row.names = F,sep="\t")


MRM.G.res.coef <- lapply(mrm.G.hp.dist,function(x){x$coef})
MRM.G.res.coef <- ldply(MRM.G.res.coef)
MRM.G.res.coef$Variable <- c("Int","phylo","habitat")
MRM.G.res.coef <- MRM.G.res.coef[-which(MRM.G.res.coef$Variable == "Int"),]
colnames(MRM.G.res.coef) <- c("id","Estimate","pval","Variable")
write.table(MRM.G.res.coef,"./_table/MRM.G.res.coef.txt",row.names = F,sep="\t",quote = F)
MRM.G.res.coef <- within(MRM.G.res.coef,Estimate <- ifelse(pval>=0.1,0,Estimate))
S <- MRM.G.res.coef %>% group_by(id) %>% summarise(S=sum(abs(Estimate)))
MRM.G.res.coef <- merge(MRM.G.res.coef,S,by="id")
MRM.G.res.coef$'Contribution' <- MRM.G.res.coef$Estimate/MRM.G.res.coef$S
MRM.G.res.coef$Contribution[is.nan(MRM.G.res.coef$Contribution)] <- 0
MRM.G.res.coef <- within(MRM.G.res.coef,id <- factor(id,levels = Levels))

pdf(file="./_figure/p28.mrm.relative.G_phylo_habitat.pdf",width = 5,height = 4)
p28 <- MRM.G.res.coef %>% ggplot(aes(x=id,y=Contribution,fill=Variable)) + geom_bar(stat="identity")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+labs(x="",y="Relative Contribution")+guides(fill=F)
print(p28)
dev.off()
pdf(file="./_figure/p29.mrm.Estimate.G_phylo_optima.pdf",width = 5,height = 4)
p29 <- MRM.G.res.coef %>% ggplot(aes(x=id,y=Estimate)) + geom_bar(aes(fill=Variable),stat="identity",position="dodge")+
  labs(y="Estimate",x="")+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+guides(fill=F)
print(p29)
dev.off()

pdf(file="./_figure/p32.smooth.hp_Z_dis.pdf",width = 5.5,height = 3.5)
p32 <- sparcc.lsp.mcmc.lcc.meta.df %>% select(.id,Z_eucli.dis,phylo_dis,habitat_dis) %>% melt.data.frame(measure.vars = c("phylo_dis","habitat_dis")) %>% within(variable <- factor(variable,levels = c("habitat_dis","phylo_dis"))) %>%
  ggplot(aes(x=value,y=Z_eucli.dis)) + geom_smooth(aes(colour=variable))+
  labs(y="Z distance",x="H/P Distance")+theme_bw()+scale_x_continuous(breaks = seq(0,1,0.2))+facet_wrap(~.id,nrow=2)+guides(colour=F)
print(p32)
dev.off()
pdf(file="./_figure/p33.smooth.hp_Z1_dis.pdf",width = 5.5,height = 3.5)
p33 <- sparcc.lsp.mcmc.lcc.meta.df %>% select(.id,z1_manha.dis,phylo_dis,habitat_dis) %>% melt.data.frame(measure.vars = c("phylo_dis","habitat_dis")) %>% within(variable <- factor(variable,levels = c("habitat_dis","phylo_dis"))) %>%
  ggplot(aes(x=value,y=z1_manha.dis)) + geom_smooth(aes(colour=variable))+
  labs(y="Z1 distance",x="H/P Distance")+theme_bw()+scale_x_continuous(breaks = seq(0,1,0.2))+facet_wrap(~.id,nrow=2)+guides(colour=F)
print(p33)
dev.off()
pdf(file="./_figure/p34.smooth.hp_Z2_dis.pdf",width = 5.5,height = 3.5)
p34 <- sparcc.lsp.mcmc.lcc.meta.df %>% select(.id,z2_manha.dis,phylo_dis,habitat_dis) %>% melt.data.frame(measure.vars = c("phylo_dis","habitat_dis")) %>% within(variable <- factor(variable,levels = c("habitat_dis","phylo_dis"))) %>%
  ggplot(aes(x=value,y=z2_manha.dis)) + geom_smooth(aes(colour=variable))+
  labs(y="Z2 distance",x="H/P Distance")+theme_bw()+scale_x_continuous(breaks = seq(0,1,0.2))+facet_wrap(~.id,nrow=2)+guides(colour=F)
print(p34)
dev.off()

pdf(file="./_figure/p35.smooth.igm.hp_Z_dis.pdf",width = 5.5,height = 3.5)
p35 <- sparcc.lsp.mcmc.lcc.meta.df %>% select(.id,Z_eucli.dis.igm,phylo_dis,habitat_dis) %>% melt.data.frame(measure.vars = c("phylo_dis","habitat_dis")) %>% within(variable <- factor(variable,levels = c("habitat_dis","phylo_dis"))) %>%
  ggplot(aes(x=value,y=Z_eucli.dis.igm)) + geom_smooth(aes(colour=variable))+
  labs(y="Z distance",x="H/P Distance")+theme_bw()+scale_x_continuous(breaks = seq(0,1,0.2))+facet_wrap(~.id,nrow=2)+guides(colour=F)
print(p35)
dev.off()
pdf(file="./_figure/p36.smooth.igm.hp_Z1_dis.pdf",width = 5.5,height = 3.5)
p36 <- sparcc.lsp.mcmc.lcc.meta.df %>% select(.id,z1_manha.dis.igm,phylo_dis,habitat_dis) %>% melt.data.frame(measure.vars = c("phylo_dis","habitat_dis")) %>% within(variable <- factor(variable,levels = c("habitat_dis","phylo_dis"))) %>%
  ggplot(aes(x=value,y=z1_manha.dis.igm)) + geom_smooth(aes(colour=variable))+
  labs(y="Z1 distance",x="H/P Distance")+theme_bw()+scale_x_continuous(breaks = seq(0,1,0.2))+facet_wrap(~.id,nrow=2)+guides(colour=F)
print(p36)
dev.off()
pdf(file="./_figure/p37.smooth.igm.hp_Z2_dis.pdf",width = 5.5,height = 3.5)
p37 <- sparcc.lsp.mcmc.lcc.meta.df %>% select(.id,z2_manha.dis.igm,phylo_dis,habitat_dis) %>% melt.data.frame(measure.vars = c("phylo_dis","habitat_dis")) %>% within(variable <- factor(variable,levels = c("habitat_dis","phylo_dis"))) %>%
  ggplot(aes(x=value,y=z2_manha.dis.igm)) + geom_smooth(aes(colour=variable))+
  labs(y="Z2 distance",x="H/P Distance")+theme_bw()+scale_x_continuous(breaks = seq(0,1,0.2))+facet_wrap(~.id,nrow=2)+guides(colour=F)
print(p37)
dev.off()

#mnpd.ses & mnhd.ses of Groups
G_mnpd_mnhd <- list()
for(i in 1:N)
{
  tmp.f <-  sparcc.lsp.G[[i]]
  tmp.p <- Node_mnpd_mnhd_degre[[i]]
  G.mnpd <- mat.dot(tmp.f,tmp.p,1)
  G.mnhd <- mat.dot(tmp.f,tmp.p,2)
  G.mnpd.mhpd <- t(rbind(G.mnpd,G.mnhd))
  colnames(G.mnpd.mhpd) <- c("mnpd.ses","mnhd.ses")
  G_mnpd_mnhd[[i]] <- G.mnpd.mhpd
}

names(G_mnpd_mnhd) <- .id


n.p.g.nmxd <- list()
tmp.n <- ldply(Node_mnpd_mnhd_degre)[,-4]
tmp.p <- ldply(Position_mnxd.list)[,-4]
tmp.g <- ldply(G_mnpd_mnhd)
colnames(tmp.n) <- c("id1","mnpd.ses","mnhd.ses")
colnames(tmp.p) <- c("id1","mnpd.ses","mnhd.ses")
colnames(tmp.g) <- c("id1","mnpd.ses","mnhd.ses")
n.p.g.nmxd$`Node-level` <- tmp.n
n.p.g.nmxd$`Motif-level` <- tmp.p
n.p.g.nmxd$`Module-level` <- tmp.g

n.p.g.nmxd.df <- ldply(n.p.g.nmxd) %>% within({id1 <- factor(id1,levels = Levels.f);.id <- factor(.id,levels = c("Module-level","Motif-level","Node-level"))})
#save(n.p.g.nmxd.df,file="./_Rdata/n.p.g.nmxd.df.Rdata")
load("./_Rdata/n.p.g.nmxd.df.Rdata")
pdf("./_figure/p38.node_motif_module.mnxd.pdf",width =6.67,height = 4)
p38 <- n.p.g.nmxd.df %>%
  ggplot(aes(x=mnpd.ses,y=mnhd.ses,color=.id))+geom_point(size=0.2,alpha=0.5)+stat_ellipse(size=0.2)+facet_wrap(~id1,nrow=2) +
  theme_bw()+labs(x="MNPD.ses",y="MNHD.ses",color="")+guides(color=guide_legend(direction = "horizontal"))+theme(legend.position="bottom")+scale_color_discrete(breaks=c("Node-level","Motif-level","Module-level"))
print(p38)
dev.off()

node_number <- n.p.g.nmxd.df %>% filter(.id=="Node-level") %>% group_by(id1) %>% count()
node.preference.in.similarity <- n.p.g.nmxd.df %>% filter(.id=="Node-level") %>% group_by(id1) %>% filter(mnpd.ses < 0 & mnhd.ses < 0) %>% count()
node.preference.in.similarity$n/node_number$n 
#[1] 0.76037 0.8000 0.84103 0.8571 0.7949 0.9128 0.74208 0.9000

#======Fig. 1==========================#
library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggpubr)
#=====Fig. 1A======
#MRM.res.coef.whole =MRM.res.coef
#MRM.res.coef.whole <- MRM.res.coef.whole[,-5]
#MRM.res.coef.module=MRM.G.res.coef
#MRM.res.coef.motif=MRM.pos.res.coef
#MRM.res.coef.whole$level="Whole network"
#MRM.res.coef.module$level="Module"
#MRM.res.coef.motif$level="Motif"
#colnames(MRM.res.coef.whole)=colnames(MRM.res.coef.motif)
#colnames(MRM.res.coef.module)=colnames(MRM.res.coef.motif)
#MRM.three_levels <- rbind.data.frame(MRM.res.coef.whole,MRM.res.coef.module,MRM.res.coef.motif)
#save(MRM.three_levels,file="./_Rdata/MRM.three_levels.Rdata")
#load("./_Rdata/MRM.three_levels.Rdata")
#MRM.three_levels.2 <-  MRM.three_levels %>% filter(.id != "Topsoil" & .id != "AMD")
#MRM.res <- cbind(MRM.three_levels.2,ldply(strsplit(as.character(MRM.res$.id),split=" ")))
#MRM.res$shape <- 21
#MRM.res$shape[MRM.res$pval > 0.05] <- 4
#MRM.res$shape[MRM.res$Variable == "habitat"] <- 24
#MRM.res$.id <- paste(MRM.res$V1,MRM.res$V2,sep=" ")
#MRM.res$bg <- "#4DFF00"
#MRM.res$bg[MRM.res$V2=="Part"] <- "#FFB300"
#MRM.res$Variable[MRM.res$Variable=="habitat"] <- "Niche"
#MRM.res$Variable[MRM.res$Variable=="phylo"] <- "Phylo."
#MRM.res <- MRM.res %>% within({V1 <- factor(V1,levels = c("SRF","DCM","MES"));
#level <- factor(level,levels=c("Motif","Module","Whole network"));
#V2 <- factor(V2,levels=c("Part","Free"));
#Variable <- factor(Variable,levels = c("Phylo","Niche"))})
#save(MRM.res,file="./_Rdata/MRM.res.Rdata")
load("./_Rdata/MRM.res.Rdata")
MRM.res$Variable <-  as.character(MRM.res$Variable)
MRM.res$level <-  as.character(MRM.res$level)
MRM.res$level[MRM.res$level=="Whole network"] <- "Whole"
MRM.res$level <- factor(MRM.res$level,levels = c("Motif","Module","Whole"))
MRM.res$Variable[MRM.res$Variable == "Phylo."] <- "Phylo"
MRM.res$Variable <- factor(MRM.res$Variable,levels = c("Phylo","Niche"))
fig.1a.shape <- ifelse(MRM.res$shape==4,NA,21)
text_size=7
Fig.1A <- MRM.res %>% ggplot(aes(x=Variable,y=V2))+geom_tile(aes(width=1), fill=MRM.res$bg,colour = "grey50")+
  geom_point(aes(fill=Contribution,size=abs(Contribution)*0.5),pch=fig.1a.shape,color="black")+theme_bw()+
  facet_grid(level~V1)+labs(x="",y="")+theme(axis.ticks = element_blank(),strip.background = element_blank())+
  scale_fill_gradient2(limits=c(-1,1),low="#FF3300",high="#3300FF",guide="colorbar")+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
  theme(axis.text = element_text(color="black",size=text_size),legend.position = "right",
        strip.text = element_text(size=text_size))+
  guides(fill=guide_colorbar(title="Relative \ncontribution",
                             title.theme = element_text(size=text_size),
                             label.theme = element_text(size=text_size),
                             barwidth = 0.7),
         size=F)
#===========================
load("./_Rdata/meta_data.df.Rdata")
#meta_data.df.2 <- meta_data.df %>% filter(id != "AMD" & id != "Topsoil")
#meta_data.df.2 <- cbind.data.frame(meta_data.df.2,ldply(strsplit(meta_data.df.2$.id.y,split="_")))
#meta_data.df.3 <- melt.data.frame(meta_data.df.2,id.vars = c("net_dis","In.same_modules","phy_sig","No.module.x","V2","V1"),measure.vars = c("phylo_dis","habitat_dis"))
#meta_data.df.3$variable <- as.character(meta_data.df.3$variable)
#meta_data.df.3$variable[meta_data.df.3$variable=="phylo_dis"] <- "Phylo."
#meta_data.df.3$variable[meta_data.df.3$variable=="habitat_dis"] <- "Niche"
#meta_data.df.3 <- meta_data.df.3 %>% within({
#  V1 <- factor(V1,levels=c("SRF","DCM","MES"));
#  variable <- factor(variable,levels=c("Phylo.","Niche"))
#})
#save(meta_data.df.3,file="./_Rdata/meta_data.df.3.Rdata")
load("./_Rdata/meta_data.df.3.Rdata")
meta_data.df.3$variable <-  as.character(meta_data.df.3$variable)
meta_data.df.3$variable[meta_data.df.3$variable == "Phylo."] <- "Phylo"
meta_data.df.3$variable <- factor(meta_data.df.3$variable,levels = c("Phylo","Niche"))

#========Fig. 1B
Fig.1B <- meta_data.df.3 %>% ggplot(aes(x=jitter(net_dis)))+
  geom_smooth(aes(y=value,colour=V1,linetype=V2),size=0.5,fill="black",alpha=0.2,fullrange=F)+
  facet_wrap(~variable)+theme_bw()+labs(x="Network path distance",y="Distance")+
  theme(axis.text = element_text(color="black",size=text_size),axis.title = element_text(color="black",size=text_size),
        strip.background = element_blank(),strip.text = element_text(size=text_size),panel.grid = element_blank()
        )+
#  guides(colour=F,linetype=F)+
    guides(colour=guide_legend(title="",title.theme = element_text(size=text_size),label.theme = element_text(size=text_size),
                            keywidth = 1,keyheight = 0.7),
         linetype=guide_legend(title="",title.theme = element_text(size=text_size),label.theme = element_text(size=text_size),
                              keywidth = 1,keyheight = 0.7))+
  scale_x_continuous(limits = c(0,14),breaks = seq(3,15,3))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))
#=======Fig. 1C=============
Fig.1C <- meta_data.df.3 %>%ggplot(aes(x=net_dis))+
  geom_smooth(aes(y=value,colour=In.same_modules,linetype=phy_sig,group=No.module.x),method="lm",se=F,size=0.4)+
  scale_colour_manual(values=c("#FF3366","#3366FF"),labels=c("Across modules","Within modules"))+
  guides(colour=guide_legend(title="",byrow=T),linetype=guide_legend(title="",byrow=T))+theme_bw()+
  scale_linetype_discrete(labels=c("P < 0.05","P > 0.05"))+
  scale_x_continuous(limits = c(0,12),breaks = seq(3,12,3))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,0.8,0.4),expand = c(0.1,0))+
  labs(x="Network path distance",y="Distance")+
  theme(axis.text = element_text(color="black",size=text_size),axis.title = element_text(color="black",size=text_size),
        panel.grid.minor = element_blank(),strip.text = element_text(size=text_size),legend.position = "bottom")+
  facet_grid(V2+variable~V1)+theme(strip.background = element_blank())+
  guides(colour=F,linetype=F)
#    guides(colour=guide_legend(title="",title.theme = element_text(size=text_size),label.theme = element_text(size=text_size)),
#         linetype=guide_legend(title="",title.theme = element_text(size=text_size),label.theme = element_text(size=text_size)))
#===Fig.1
Fig.1 <- arrangeGrob(Fig.1A,Fig.1B,Fig.1C,ncol=2,nrow=5,layout_matrix = cbind(c(1,1,1,2,2), c(3,3,3,3,3)))

Fig.1 <- as_ggplot(Fig.1) +                          
  draw_plot_label(label = c("A", "B", "C"), size = 9,
                  x = c(0, 0, 0.5), y = c(1, 0.4, 1))

pdf("./_figure/Fig.1.pdf",height = 4, width = 6.6)
print(Fig.1)
dev.off()

tiff("./_figure/Fig.1.tiff",height = 4, width = 6.6,units = "in",res=300)
print(Fig.1)
dev.off()

#============Fig. 2=================
#====Fig.2A.bc
#test.df <- data.frame(null=rnorm(1000,0.5,0.01))
#data <- data.frame(x=c(-2,-2),y=c(0,0.2))
#p <- test.df %>% ggplot(aes(null)) + geom_density(fill="#00E6FF")+
#  scale_y_continuous(expand = c(0,0,0,0.02))+
#  theme(axis.text.y= element_blank(),axis.ticks = element_blank(),axis.line.y = element_blank())+
#  labs(x="Null distribution of MNPD (or MNHD) values",y="")+geom_line(aes(x=c(0.48,0.48),y=c(0,16)),data=as.data.frame(data),linetype="dotdash")+
#  annotate(geom="text",x=0.48,y=18,label=c("obs.MNPD\n(or obs.MNND)"),color="red",size=3)
#pdf("./_figure/Fig.2A.2.pdf",height = 5,width=6)
#print(p)
#dev.off()
#test.df2 <- data.frame(null=(0.48-test.df$null)/sd(test.df$null))
#data2 <- data.frame(x=c(-2,-2),y=c(0,0.2))
#p2 <- test.df2 %>% ggplot(aes(null)) + geom_density(fill="#E6FF00")+
#  scale_y_continuous(expand = c(0,0,0,0.02))+
#  theme(axis.text.y= element_blank(),axis.ticks = element_blank(),axis.line.y = element_blank())+
#  labs(x="Null distribution of MNPD (or MNHD) value",y="")+geom_line(aes(x=c(-2,-2),y=c(0,0.40)),data=as.data.frame(data2),linetype="dotdash")+
#  annotate(geom="text",x=-2,y=0.45,label=c("ses.MNPD\n(or ses.MNND)"),color="red",size=3)
#pdf("./_figure/Fig.2A.3.pdf",height = 5,width=6)
#print(p2)
#dev.off()
#====Fig.2A
Fig.2A <- ggdraw() + draw_image("./_figure/Fig.2.A.tif", scale = 0.9)
#====Fig.2B
load("./_Rdata/n.p.g.nmxd.df.Rdata")
n.p.g.nmxd.df2 <- n.p.g.nmxd.df %>% filter(id1 != "AMD" & id1 !="Topsoil")
n.p.g.nmxd.df3 <- cbind(n.p.g.nmxd.df2,ldply(strsplit(as.character(n.p.g.nmxd.df2$id1),split="_")))
n.p.g.nmxd.df3$.id2 <- as.character(n.p.g.nmxd.df3$.id)
n.p.g.nmxd.df3$.id2[n.p.g.nmxd.df3$.id2=="Node-level"] <- "Node"
n.p.g.nmxd.df3$.id2[n.p.g.nmxd.df3$.id2=="Motif-level"] <- "Motif"
n.p.g.nmxd.df3$.id2[n.p.g.nmxd.df3$.id2=="Module-level"] <- "Module"
n.p.g.nmxd.df3 <- n.p.g.nmxd.df3 %>% within({
  V1 <- factor(V1,levels = c("SRF","DCM","MES"));
  V2 <- factor(V2,levels = c("Free","Part"));
  .id2 <- factor(.id2,levels=c("Node","Motif","Module"))
})
Fig.2B <- n.p.g.nmxd.df3 %>%
  ggplot(aes(x=mnpd.ses,y=mnhd.ses,color=.id2))+geom_point(size=1,pch=18,alpha=0.2)+stat_ellipse(size=0.4)+facet_wrap(~id1,nrow=2) +
  theme_bw()+labs(x="ses.MNPD",y="ses.MNHD",color="")+
  guides(color=guide_legend(direction = "vertical",
                            label.position="bottom"
                            ))+
  theme(legend.position="right",panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(color="black",size=text_size),
        axis.text = element_text(color="black",size=text_size),
        axis.ticks = element_line(color="black",size=0.3),
        legend.text = element_text(color="black",size=text_size))+
  scale_color_discrete(breaks=c("Node","Motif","Module")) + 
  facet_grid(V2~V1)

#=======Fig.2C
#===Fig.2C.a
color_blocks <- data.frame(x=rep(c(-6,-2,-2,2,2,6),3),
                           ymax=c(rep(6,6),rep(2,6),rep(-2,6)),
                           ymin=c(rep(2,6),rep(-2,6),rep(-6,6)),
                           col=factor(as.vector(rbind(seq(1,9),seq(1,9)))))
Fig.2C.a <- color_blocks %>% ggplot(aes(x=x))+geom_ribbon(aes(ymin=ymin,ymax=ymax,fill=factor(col)))+
  scale_fill_manual(values = c("#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00"))+
  scale_x_continuous(expand = c(0,0),breaks = c(-2,2))+
  scale_y_continuous(expand = c(0,0),breaks = c(-2,2))+
  labs(x="ses.MNPD",y="ses.MNND")+guides(fill=FALSE)+
  theme(axis.text = element_text(color="black",size=text_size),
        axis.ticks = element_line(color="black"),
        axis.title = element_text(color="black",size=text_size))

#====Fig.2C.b
sep_ses=function(v1,v2)
{
  a=2;b=-2
  res <- rep(0,length(v1))
  res[(v1 < b & v2 > a)] <- 1;res[(v1 > b & v2 > a)] <- 2;res[(v1 > a & v2 > a)] <- 3;
  res[(v1 < b & v2 < a)] <- 4;res[(v1 > b & v2 < a)] <- 5;res[(v1 > a & v2 < a)] <- 6;
  res[(v1 < b & v2 < b)] <- 7;res[(v1 > b & v2 < b)] <- 8;res[(v1 > a & v2 < b)] <- 9;
  return(res)
}
n.p.g.nmxd.df3$type=sep_ses(n.p.g.nmxd.df3$mnpd.ses,n.p.g.nmxd.df3$mnhd.ses)
n.p.g.nmxd.df4 <- n.p.g.nmxd.df3 %>% group_by(.id,V1,V2) %>% dplyr::count(type) %>% 
  group_by(.id,V1,V2) %>% do(data.frame(.id=.$.id,type=.$type,n=.$n,V1=.$V1,V2=.$V2,S=sum(.$n)))
n.p.g.nmxd.df4$V3 <- paste(n.p.g.nmxd.df4$V1,n.p.g.nmxd.df4$V2,sep=" ")
n.p.g.nmxd.df4$.id <- as.character(n.p.g.nmxd.df4$.id)
n.p.g.nmxd.df4$proportions <- n.p.g.nmxd.df4$n/n.p.g.nmxd.df4$S
n.p.g.nmxd.df4$.id[n.p.g.nmxd.df4$.id=="Motif-level"] <- "Motif"
n.p.g.nmxd.df4$.id[n.p.g.nmxd.df4$.id=="Module-level"] <- "Module"
n.p.g.nmxd.df4$.id[n.p.g.nmxd.df4$.id=="Node-level"] <- "Node"
n.p.g.nmxd.df4 <- n.p.g.nmxd.df4 %>% within({
  .id <- factor(.id,levels = c("Node","Motif","Module"));
  V3 <- factor(V3,levels = c("SRF Free","SRF Part","DCM Free","DCM Part","MES Free","MES Part"));
  type <- factor(type,levels=seq(1,9))
})
Fig.2C.b <- n.p.g.nmxd.df4 %>% ggplot(aes(x=V3,y=proportions,fill=type))+
  geom_bar(stat="identity")+
  labs(x="",y="Proportions")+theme_bw()+
  theme(axis.text.x = element_text(color="black",angle = 45,vjust=1,hjust=1,size=text_size),
        strip.background = element_blank(),
        axis.text.y = element_text(color="black",size=text_size),
        axis.title.y = element_text(color="black",size=text_size),
        axis.ticks = element_line(color="black",size=0.3))+
  scale_x_discrete(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c("#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00"))+
  facet_wrap(~.id)+guides(fill=FALSE)

#=======Fig.2C
Fig.2C <- arrangeGrob(Fig.2C.a,Fig.2C.b,ncol=4,nrow=7,layout_matrix = cbind(c(3,1,1,1,1,1,3),
                                                                            c(2,2,2,2,2,2,2),
                                                                            rep(2,7),rep(2,7)))
#=======Fig.2
Fig.2 <- arrangeGrob(Fig.2A,Fig.2B,Fig.2C,ncol=5,nrow = 2,layout_matrix = rbind(c(1,1,2,2,2),c(1,1,3,3,3)))
Fig.2 <- as_ggplot(Fig.2)+draw_plot_label(label=LETTERS[1:3],x=c(0,0.4,0.4),y=c(1,1,0.5),size=9)
pdf("./_figure/Fig.2.pdf",height = 4,width = 6.6)
print(Fig.2)
dev.off()

tiff("./_figure/Fig.2.tif",units = "in",height = 4,width = 6.6,res=300)
print(Fig.2)
dev.off()

#=======Fig.3==========
#===Fig.3a
load("./_Rdata/randomforest.res.explaination.Rdata")
randomforest.res.explaination2 <- randomforest.res.explaination %>% filter(L1 != "AMD" & L1 != "Topsoil")
randomforest.res.explaination2$L1 <- as.character(randomforest.res.explaination2$L1)
randomforest.res.explaination2 <- cbind(randomforest.res.explaination2,
                                        ldply(strsplit(randomforest.res.explaination2$L1,split="_")))
randomforest.res.explaination2 <- randomforest.res.explaination2 %>% within({
V1 <- factor(V1,levels = c("SRF","DCM","MES"));
variable <- factor(variable,levels=c("phylo_dis","habitat_dis"))
})
a_size=1.3
Fig.3a <- 
  randomforest.res.explaination2 %>% ggplot(aes(x=as.numeric(net_dis),y=value),color="black") + 
#geom_line(aes(colour=variable))+
    geom_smooth(aes(colour=variable),se=FALSE)+
  theme_bw()+
  theme(axis.text = element_text(color="black",size=text_size*a_size),
        axis.title = element_text(color="black",size=text_size*a_size),
        strip.text.x = element_text(color = "white",size=text_size*a_size),
        strip.text.y = element_text(color="black",size=text_size*a_size),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  scale_color_manual(values = c("#FFB300","#4DFF00"))+
  facet_grid(V2~V1)+
  guides(colour=FALSE,fill=FALSE)+
  geom_point(aes(fill=variable),pch=21,color="black")+
  scale_fill_manual(values = c("#FFB300","#4DFF00"))+
  scale_x_continuous(limits = c(0,15),breaks = seq(3,15,3))+
  labs(x="",y="Relative contribution")

#====Fig.3b
load("./_Rdata/sparcc.lsp.pearson.cor.Rdata")
sparcc.lsp.pearson.cor2 <- sparcc.lsp.pearson.cor %>% filter(.id != "AMD" & .id != "Topsoil")
sparcc.lsp.pearson.cor2$.id <-as.character(sparcc.lsp.pearson.cor2$.id)
sparcc.lsp.pearson.cor2 <- cbind(sparcc.lsp.pearson.cor2,ldply(strsplit(sparcc.lsp.pearson.cor2$.id,split = "_")))
sparcc.lsp.pearson.cor2 <- sparcc.lsp.pearson.cor2 %>% select(.id,net_dis,Cor.igm.Z.phylo,Cor.igm.Z.habitat,V1,V2) %>%
  melt.data.frame(measure.vars = c("Cor.igm.Z.phylo","Cor.igm.Z.habitat")) %>% 
  within({variable <- factor(variable,levels=c("Cor.igm.Z.habitat","Cor.igm.Z.phylo"));
  .id <- factor(.id,levels = Levels.f);
  V1 <- factor(V1,levels = c("SRF","DCM","MES"));
  variable <- factor(variable,levels=c("Cor.igm.Z.phylo","Cor.igm.Z.habitat"))})

Fig.3b <- sparcc.lsp.pearson.cor2 %>% 
  ggplot(aes(x=net_dis,y=value)) + 
  geom_bar(aes(group=variable,fill=variable),color="black",stat="identity",position="dodge")+
  labs(x="Network path distance",y="Pearson' r")+
  facet_grid(V2~V1)+
  theme_bw()+
  scale_x_continuous(breaks = seq(1,6))+
  theme(axis.text = element_text(color="black",size=text_size*a_size),
        axis.title = element_text(color="black",size=text_size*a_size),
        strip.background = element_blank(),
        strip.text.x = element_text(color = "white",size=text_size*a_size),
        strip.text.y = element_text(color="black",size=text_size*a_size),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")+
#geom_line(aes(colour=variable))+
    geom_smooth(aes(x=net_dis,y=value,colour=variable))+
  guides(colour=FALSE,fill=FALSE)+
  scale_colour_manual(values=c("#FFB300","#4DFF00"))+
  scale_fill_manual(values=c("#FFB300","#4DFF00"),labels=c("Phylo","Niche"))

#==Fig.3
#Fig.3 <- arrangeGrob(Fig.3a,Fig.3b,ncol=1,nrow = 7,layout_matrix = cbind(c(1,1,1,2,2,2,2)))
Fig.3 <- arrangeGrob(Fig.3a,Fig.3b)
Fig.3 <- as_ggplot(Fig.3)+draw_plot_label(label=LETTERS[1:2],x=c(0,0),y=c(1,0.57),size=9)
pdf("./_figure/Fig.3.pdf",height = 5,width = 6)
print(Fig.3)
dev.off()

tiff("./_figure/Fig.3.tif",units = "in",height = 6,width = 6.6,res=300)
print(Fig.3)
dev.off()

#=====Fig.S1==========
#======Fig.S1a
load("./_Rdata/sparcc.ls.res.Rdata")
sparcc.ls.bic <- ldply(sparcc.ls.res) %>%  within({.id <- factor(.id,levels = Levels.f)})
sparcc.ls.bic2 <- sparcc.ls.bic %>% filter(.id != "AMD" & .id != "Topsoil")
sparcc.ls.bic2$.id <- as.character(sparcc.ls.bic2$.id)
sparcc.ls.bic2 <- cbind(sparcc.ls.bic2,ldply(strsplit(sparcc.ls.bic2$.id,split = "_")))
sparcc.ls.bic2$V3 <- paste(sparcc.ls.bic2$V1,sparcc.ls.bic2$V2,sep=" ")
sparcc.ls.bic2 <- sparcc.ls.bic2 %>% within({
  V3 <- factor(V3,levels = c("SRF Free","DCM Free","MES Free","SRF Part","DCM Part","MES Part"));
  
})

Fig.S1a <- sparcc.ls.bic2 %>%
  ggplot(aes(x=No.Groups,y=`Overall BIC`)) + 
  geom_point(size=1) + 
  geom_line()+labs(x="The number of modules")+
  scale_x_continuous(limits = c(0,30),breaks = seq(0,30,5))+
  facet_wrap(~V3,nrow=2,scales="free")+theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color="black"),
        axis.title = element_text(color="black"),
        strip.background = element_blank())

#======Fig.S1b
load("./_Rdata/sparcc.lsp.pearson.cor.Rdata")
sparcc.lsp.pearson.cor3 <- sparcc.lsp.pearson.cor %>% filter(.id != "AMD" & .id != "Topsoil")
sparcc.lsp.pearson.cor3$.id <-as.character(sparcc.lsp.pearson.cor3$.id)
sparcc.lsp.pearson.cor3 <- cbind(sparcc.lsp.pearson.cor3,ldply(strsplit(sparcc.lsp.pearson.cor3$.id,split = "_")))
Fig.S1b <- sparcc.lsp.pearson.cor3 %>% select(.id,net_dis,Cor.Z.phylo,Cor.Z.habitat,V1,V2) %>%
  melt.data.frame(measure.vars = c("Cor.Z.phylo","Cor.Z.habitat")) %>% 
  within({variable <- factor(variable,levels=c("Cor.Z.habitat","Cor.Z.phylo"));
  .id <- factor(.id,levels = Levels.f);
  V1 <- factor(V1,levels = c("SRF","DCM","MES"));
  variable <- factor(variable,levels=c("Cor.Z.phylo","Cor.Z.habitat"))}) %>%
  ggplot(aes(x=net_dis,y=value)) + 
  geom_bar(aes(group=variable,fill=variable),color="black",stat="identity",position="dodge")+
  labs(x="Network path distance",y="Pearson' r")+
  facet_grid(V2~V1)+
  theme_bw()+
  scale_x_continuous(breaks = seq(1,6))+
  theme(axis.text = element_text(color="black"),
        axis.title = element_text(color="black"),
        strip.background = element_blank(),
        strip.text.x = element_text(color = "black"),
        strip.text.y = element_text(color="black"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")+
  geom_smooth(aes(x=net_dis,y=value,colour=variable))+
  guides(colour=FALSE,fill=guide_legend(title="",byrow=T))+
  scale_colour_manual(values=c("#FFB300","#4DFF00"))+
  scale_fill_manual(values=c("#FFB300","#4DFF00"),labels=c("Phylo","Niche"))

#======
pdf("./_figure/Fig.S1.pdf",width = 5,height = 6)
plot_grid(Fig.S1a,Fig.S1b,nrow=2,labels="AUTO")
dev.off()

tiff("./_figure/Fig.S1.tif",units = "in",res=300,width = 5,height = 6)
plot_grid(Fig.S1a,Fig.S1b,nrow=2,labels="AUTO")
dev.off()

#==========Fig.4
library(igraph)
m<- matrix(0,50,50)
pdf("./_figure/Fig.4B.1.pdf")
plot(graph_from_adjacency_matrix(m),vertex.label="",size=20)
dev.off()

library(reshape)
Fig.4C.1 <- melt(adj_mat) %>% 
  ggplot(aes(x=X1,y=X2)) + 
  geom_tile(fill="white",color="grey")+
  geom_text(aes(label=value),size=2)+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
pdf("./_figure/Fig.4C.1.pdf")
print(Fig.4C.1)
dev.off()

net=network(res$adj_mat[[83]],directed = F)
fig.4C.xx <- ggnet2(net ,node.color = "#FFCC33",
                    edge.color="black",edge.size=1,
                    node.size = 10)
pdf("./_figure/Fig.4C.83.post_steady.pdf")
print(fig.4C.xx)
dev.off()


Fig.4D <- melt.data.frame(link.dynamic,measure.vars = c("TL","LL","GL")) %>% 
  ggplot(aes(x=Time,y=value,colour=variable))+
  geom_line(size=0.8)+theme_classic()+
  labs(x="Time (t)",y="Links")+
  guides(colour=guide_legend(title=""))+
  scale_colour_discrete(labels=c("Total","Loss","Gain"))+
  geom_vline(xintercept = 69,linetype="dashed")+
  theme(axis.text = element_text(color="black",size=12),
        axis.title = element_text(color="black",size=15),
        axis.ticks = element_line(color="black"),
        legend.position = c(0.9,0.55))
pdf("./_figure/Fig.4D.pdf",width = 4,height = 3)
print(Fig.4D)
dev.off()

#========Fig. 5===========#
#===Fig.5A====
load("./_Rdata/data.0.Rdata")
data.0$variable <- as.character(data.0$variable)
data.0.2 <- data.0 %>% filter(variable != "M2/M1")
data.0.2$variable[data.0.2$variable=="G"] <- "Gain Links"
data.0.2$variable[data.0.2$variable=="L"] <- "Loss Links"
Fig.5A <- data.0.2 %>% 
  ggplot(aes(x=l,y=value*100,group=Nspecies,colour=Nspecies))+geom_smooth(se=T)+theme_bw()+
  labs(x="The distance of long-lasting effects (l)",y="Percentage %",color="")+
  theme(axis.text = element_text(color = "black",size=text_size),
        axis.title = element_text(color = "black",size=text_size),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black",size=text_size))+
  scale_x_discrete(expand = c(0.1,0.1,0.1,0),labels=seq(0,5)) +
  facet_wrap(~variable,scales = "free_y")
#====Fig.5B=====
#save(phylo_sig_gain_loss.df,file="./_Rdata/phylo_sig_gain_loss.df.Rdata")
load("./_Rdata/phylo_sig_gain_loss.df.Rdata")
Fig.5B <- phylo_sig_gain_loss.df %>% melt.data.frame(id.vars = c("phylo_isg"),measure.vars = c("Gain","Loss")) %>%
  ggplot(aes(x=phylo_isg,y=value*100))+geom_smooth()+facet_wrap(~variable,scale="free_y")+
  theme_bw()+labs(x="Phylogenetic signal in optimal niche",y="Percentage %")+
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(color="black",size=text_size),
        axis.title = element_text(color="black",size=text_size),
        axis.ticks = element_line(color="black"))
#===Fig.5
Fig.5 <- arrangeGrob(Fig.5A,Fig.5B,ncol=1,nrow = 2,layout_matrix = cbind(c(1,2)))
Fig.5 <- as_ggplot(Fig.5)+draw_plot_label(label=c("A","B"),x=c(0,0),y=c(1,0.57),size=10)
pdf("./_figure/Fig.5.s.pdf",width = 6.6,height = 4)
print(Fig.5)
dev.off()

