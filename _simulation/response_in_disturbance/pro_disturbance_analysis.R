setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./functions.R")
library(stringr)
library(plyr)
library(ggplot2)
library(dplyr)
pro_distubance.var_spec <- list()
rdata <- dir(getwd(),pattern=".Rdata")
N_var <- str_extract(rdata,"N[0-9]*")
for(i in 1:length(rdata))
{
  load(rdata[[i]])
  pro_distubance.var_spec[[i]] <- variable_Nspecies
  print("Sleeping for 3s")
  Sys.sleep(3)
}
names(pro_distubance.var_spec) <- N_var
pro_distubance.var_spec2 <- list()
for(i in 1:length(rdata))
{
  tmp <- pro_distubance.var_spec[[i]]
  names(tmp) <- seq(1,6)
  tmp <- lapply(tmp,pro_df)
  tmp.df <- ldply(tmp)
  colnames(tmp.df)[1] <- "l"
  pro_distubance.var_spec2[[i]] <- tmp.df
}

names(pro_distubance.var_spec2) <- N_var
pro_distubance.var_spec2.df <- ldply(pro_distubance.var_spec2)
colnames(pro_distubance.var_spec2.df) <- c("Nspecies","l","G","L","M2/M1","MD")
pro_distubance.var_spec2.df= pro_distubance.var_spec2.df %>% filter(Nspecies != "N10")
data.0 <- pro_distubance.var_spec2.df %>% within(Nspecies <- factor(Nspecies,levels = c("N30","N50","N70","N100"))) %>% 
 melt.data.frame(measure.vars = c("G","L","M2/M1","MD"))
save(data.0,file="../../_Rdata/data.0.Rdata")
load("../../_Rdata/data.0.Rdata")

p <- data.0 %>% 
  ggplot(aes(x=l,y=value,group=Nspecies,colour=Nspecies))+geom_smooth()+theme_bw()+
  labs(x="l",y="Network persistence",color="")+theme(axis.text = element_text(color = "black" ),axis.title.x = element_text(face="italic"))+
  scale_x_discrete(expand = c(0.1,0.1,0.1,0),labels=seq(0,5)) +
  facet_wrap(~variable,scales = "free_y")


pdf("../../_figure/Simulation.p1.Network_persitence.pdf",width = 8,height = 3)
print(p)
dev.off()

F.value <- data.0 %>%group_by(variable,l) %>% do(data.frame(F=my_aov(.$value,.$Nspecies,1),P=my_aov(.$value,.$Nspecies,2))) %>% 
  recast(l~variable,measure.var = "F")
P.value <- data.0 %>%group_by(variable,l) %>% do(data.frame(F=my_aov(.$value,.$Nspecies,1),P=my_aov(.$value,.$Nspecies,2))) %>% 
  recast(l~variable,measure.var = "P")
p.sig <- matrix("",6,3)
for(i in 1:6)
{
  for(j in 1:3)
  {
    if(P.value[i,j+1] < 0.001) {p.sig[i,j] <- "***";next}
    if(P.value[i,j+1] < 0.01) {p.sig[i,j] <- "**";next}
    if(P.value[i,j+1] < 0.05) {p.sig[i,j] <- "*";next}
  }
}


F.value <- format(F.value[,2:4],digits=3)
anova.res <- matrix(paste(as.matrix(F.value),p.sig,sep=""),6)
colnames(anova.res) <- c("G","L","M2/M1")
rownames(anova.res) <- seq(1:6)-1
write.table(data.frame("l"=rownames(anova.res),anova.res),quote = F,sep="\t",row.names = FALSE,file="../../_table/anova.network_persistence.txt")

sink('../../_table/ANOVA.network_persistence~lxNspecies.txt')
cat("Gained links:\n")
cat("------------\n")
tmp <- data.0 %>% filter(variable=="G");summary(aov(tmp$value~tmp$Nspecies*tmp$l))
cat("===============\n")
cat("Loss links:\n")
cat("------------\n")
tmp <- data.0 %>% filter(variable=="L");summary(aov(tmp$value~tmp$Nspecies*tmp$l))
cat("===============\n")
cat("M2/M1:\n")
cat("------------\n")
tmp <- data.0 %>% filter(variable=="M2/M1");summary(aov(tmp$value~tmp$Nspecies))
cat("===============\n")
sink()

#=====Phylo_sig======
phylo_sig_disturbance <- load("./Pro_disturbance.N50.phy_sig.Rdata")
phylo_sig_gain_loss <- list()
phylo_isg <- seq(0,0.2,0.02)
for(i in 1:length(variable_Nspecies))
{
  tmp <- pro_df(variable_Nspecies[[i]])
  tmp$phylo_isg <- phylo_isg[i]
  phylo_sig_gain_loss[[i]] <- tmp
}

phylo_sig_gain_loss.df <- ldply(phylo_sig_gain_loss)
save(phylo_sig_gain_loss.df,file="./_Rdata/phylo_sig_gain_loss.df.Rdata")
phylo_sig_gain_loss.df %>% melt.data.frame(id.vars = c("phylo_isg"),measure.vars = c("Gain","Loss")) %>%
  ggplot(aes(x=phylo_isg,y=value*100))+geom_smooth()+facet_wrap(~variable,scale="free_y")+
  theme_bw()+labs(x="Phlo Sig",y="Percentage %")+
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
