# this function is for phylogenetic signal test
# As it have been tested for many data sets so far,
# the default setting can be followed, without any changes.
# jianjun wang, 2012-11-01

# phylogenetic signal plotting ####

# standardize the otu.scores
otu.scores.std = otu.scores; # making new matrix that will have standardized species scores
otu.scores.std = decostand(otu.scores.std, "standardize", MARGIN = 2)  

match.phylo.otu = match.phylo.data_JJWang(phylo, otu.scores.std) ; 

# making sure the names on the phylogeny are ordered the same as the names in otu table
match.phylo.otu[[1]]; # tree
dim(match.phylo.otu[[2]]); match.phylo.otu[[2]][1:5, 1:2] ;

phylo.dist.4.reg = cophenetic(match.phylo.otu[[1]]);   
phylo.dist.4.reg = phylo.dist.4.reg/max(phylo.dist.4.reg) ; 

#*******************************************************************************
# Use "all.env.variables" ####
all.env.variables = c("all.vars", colnames(match.phylo.otu[[2]]))   
# go across all environmental variables
# "all.vars" means to combine all variables

for (ii in 1:1) {   
  # length(all.env.variables); if you want to go across all environmental variables; 
  # Slow processes, pay attention #
  (env.variable = all.env.variables[ii]);     
  
  ##########################
  ## calculate the meadian of ecological differences within phylogenetic bins (bands) ####
  
  phylo.dist.4.reg = as.dist(phylo.dist.4.reg); head(phylo.dist.4.reg); 
  # making the distance matrix work for a regression
  
  if (env.variable != 'all.vars') {
    spp.dist = as.matrix(dist(match.phylo.otu[[2]][,env.variable])); 
    # generate distance matrix based on otu habitat values
    colnames(spp.dist) = match.phylo.otu$phy$tip.label; 
    rownames(spp.dist) = match.phylo.otu$phy$tip.label; 
    #spp.dist[1:5,1:5];
    spp.dist.4.reg = as.dist(spp.dist); head(spp.dist.4.reg); # changing the format
  } else{};
  
  if (env.variable == 'all.vars') {
    spp.dist = as.matrix(dist(match.phylo.otu[[2]][,env.vars])); 
    colnames(spp.dist) = match.phylo.otu$phy$tip.label; 
    rownames(spp.dist) = match.phylo.otu$phy$tip.label; 
    #spp.dist[1:5,1:5];
    spp.dist.4.reg = as.dist(spp.dist); head(spp.dist.4.reg);
  } else{};
  

  ##################################
  ## for mantel correlog        ####
  ##################################
  
  if (env.variable != 'all.vars') {
    spp.dist = as.matrix(dist(match.phylo.otu[[2]][,env.variable])); 
    # generate distance matrix based on otu habitat values
    colnames(spp.dist) = match.phylo.otu$phy$tip.label; 
    rownames(spp.dist) = match.phylo.otu$phy$tip.label; #spp.dist[1:5,1:5];
  } else{};
  
  if (env.variable == 'all.vars') {
    spp.dist = as.matrix(dist(match.phylo.otu[[2]][,env.vars])); 
    colnames(spp.dist) = match.phylo.otu$phy$tip.label; 
    rownames(spp.dist) = match.phylo.otu$phy$tip.label; 
  } else{};
  
  
  # head(phylo.dist.4.correlog); # making the distance matrix work with the mantel.correlog function
  
  spp.dist.4.correlog = as.matrix(spp.dist) ;          rm(spp.dist)  # release memory
  phylo.dist.4.correlog = as.matrix(phylo.dist.4.reg); rm(phylo.dist.4.reg)  # release memory
  
  phylo.sig.correlog = mantel.correlog(spp.dist.4.correlog,
                                       phylo.dist.4.correlog,
                                       nperm=permutations,
                                       cutoff=FALSE,
                                       n.class=50,
                                       mult="bonferroni"); 
  
  write.csv(phylo.sig.correlog$mantel.res,
            paste("Signal/",env.variable,"_Phylo_Correlogram_rarefy", rrarefy.number, ".csv", sep=""), quote=F);
  
  # re-read the data and plot
  phylo.sig = read.csv(paste("Signal/",env.variable,"_Phylo_Correlogram_rarefy", rrarefy.number, ".csv", sep=""),head=T,row.names=1)
  
  pdf(paste("Signal/",env.variable,"_phylo.sig_rarefy", rrarefy.number, ".pdf",sep="")) ;
  plot(phylo.sig[,c(1,3)],
       xlab="Phylogenetic Distance Class",
       ylab="Mantel Test Statistic") ;
       title(main = paste(type.of.score, env.variable, sep="  "), cex=0.9) 
  lines(phylo.sig[,c(1,3)]) ;
  points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4) ;
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) ;
  abline(h=0, lty=2, col="red")   
  dev.off()
  
  print(c(subset(colnames(as.matrix(phylo.dist.4.correlog)),colnames(as.matrix(phylo.dist.4.correlog))!=colnames(as.matrix(spp.dist.4.correlog))))); ## should be empty if all aligned correctly
  
  print(paste(env.variable,"finished, ",date(),sep=" "))
  
}

# the end


