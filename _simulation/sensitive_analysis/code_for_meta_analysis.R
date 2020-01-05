library(randomForest)
library(reshape)
library(igraph)
library(dplyr)
library(plyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../simulation_function.R")
#setwd("../sensitive_analysis/")
#load("sensitive_analysis_1_30.Rdata")
#data <- variable_Nspecies
#load("sensitive_analysis_31_40.Rdata")
#data[31:40] <- variable_Nspecies
#load("sensitive_analysis_41_50.Rdata")
#data[41:50] <- variable_Nspecies
#load("sensitive_analysis_51_60.Rdata")
#data[51:60] <- variable_Nspecies
#load("sensitive_analysis_61_70.Rdata")
#data[61:70] <- variable_Nspecies
#load("sensitive_analysis_71_80.Rdata")
#data[71:80] <- variable_Nspecies
#load("sensitive_analysis_81_90.Rdata")
#data[81:90] <- variable_Nspecies
#load("sensitive_analysis_91_100.Rdata")
#data[91:100] <- variable_Nspecies
#load("sensitive_anlysis_meta.Rdata")

#res.tmp <-lapply(data,rf.func3)
#names(res.tmp) <- seq(1,100)
#mydata <- ldply(res.tmp)
#mydata$net_dis <- as.numeric(mydata$net_dis)
#vars<- cbind(vars,seq(1,100))
#colnames(vars) <- c("N","l","connection","disconnection","a","th","phy_sig",".id")
#mydata2 <- merge.data.frame(mydata,vars,by.x=".id",by.y=".id",all=T)

#save(mydata2,file="mydata2.Rdata")
load("mydata2")

mydata2.phy_dis <- mydata2 %>% filter(variable=="phy_dis")
mydata2.phy_dis.max_N <- max(mydata2.phy_dis$net_dis)

mydata2.habitat_dis <- mydata2 %>% filter(variable=="habitat_dis")
mydata2.habitat_dis.max_N <- max(mydata2.phy_dis$net_dis)

glm.res.phy_dis <- NULL
for(i in 1:mydata2.phy_dis.max_N)
{
  net_dis.tmp <- as.character(i)
  data.tmp <- mydata2.phy_dis %>% filter(net_dis == net_dis.tmp)
  glm.res.phy_dis[[i]] <- summary(glm(m~N+l+a+connection+disconnection+th+phy_sig,data=data.tmp))
}

glm.res.habitat_dis <- NULL
for(i in 1:mydata2.habitat_dis.max_N)
{
  net_dis.tmp <- as.character(i)
  data.tmp <- mydata2.habitat_dis %>% filter(net_dis == net_dis.tmp)
  glm.res.habitat_dis[[i]] <- summary(glm(m~N+l+a+connection+disconnection+th+phy_sig,data=data.tmp))
}

sink(file="../../_table/glm.res.habitat.txt")
glm.res.habitat_dis
sink()
sink(file="../../_table/glm.res.phy.txt")    
glm.res.phy_dis    
sink()