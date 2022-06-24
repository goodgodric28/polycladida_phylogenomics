####################################
# ancestral_state_polyclad.R
# Written by: Jessica A. Goodheart
# Last Updated: 24 June 2022
# Purpose: To perform ancestral state reconstruction on polyclad phylogeny
# Inputs used: Polycladida phylogeny and life history traits file
####################################

#########################################
# Script setup
#########################################

directory <- "~/[PATH_TO]/polyclad_ASR"
setwd(directory)

library(corHMM)
library(ape)

#########################################
# Tree setup 
#########################################

####################
### GENUS LEVEL ####
####################

# Assign correct tree topology to the variable tree
tree <- read.tree("polyclad_37_4469_part_besttree_w_500BS_taxfixed.tre")

# Here I am dropping all outgroup tips to only keep one of each genus
tree <- drop.tip(tree, "Armatoplana_leptalea_B")
tree <- drop.tip(tree, "Prostheceraeus_vittatus")
tree <- drop.tip(tree, "Prosthiostomum_siphunculus")
tree <- drop.tip(tree, "Prorhynchus_sp_I")
tree <- drop.tip(tree, "Notocomplana_lapunda")

# Removing species names to leave only genus
tree$tip.label<-gsub("_.*","",tree$tip.label)

# To view tree
plot.phylo(tree)

# Create an ultrametric tree with correlated model, output of chronos 
# needed to be rewritten using read.tree, then test for ultrametricity
chronotree <- chronos(tree, lambda=1, model = "correlated")
ultratree <- read.tree(text=write.tree(chronotree))
is.ultrametric(ultratree)

ultratree <- rotate(ultratree,36)
ultratree <- rotate(ultratree,49)
ultratree <- rotate(ultratree,38)
ultratree <- rotate(ultratree,39)
ultratree <- rotate(ultratree,40)
ultratree <- rotate(ultratree,41)
ultratree <- rotate(ultratree,44)
ultratree <- rotate(ultratree,45)
ultratree <- rotate(ultratree,50)
ultratree <- rotate(ultratree,52)
ultratree <- rotate(ultratree,53)
ultratree <- rotate(ultratree,54)
ultratree <- rotate(ultratree,55)
ultratree <- rotate(ultratree,56)
ultratree <- rotate(ultratree,57)
ultratree <- rotate(ultratree,58)

########################################
# Ancestral state reconstruction 
########################################

# Pull the data matrix of interest into R
lh_data <- read.csv("life_history_traits_full_polycladida_genus.csv",stringsAsFactors=FALSE)

# Change any spaces into underscores for species names, then generate a vector with 
# the appropriate feeding data where the names are set to the species names
lh_data$genus<-gsub(" ","_",x=lh_data$genus)
lh_data2 <- data.frame(lh_data$genus,as.character(lh_data$devel.code))

lh_data2 <- lh_data2$as.character.lh_data.devel.code.
lh_data2 <- setNames(lh_data2, lh_data$genus)
datamatch <- match(names(lh_data2),ultratree$tip.label)
lh_data2 <- lh_data2[order(datamatch)]
lh_data3 <- data.frame(names(lh_data2),lh_data2)
lh_data3$lh_data2 <- as.character(lh_data3$lh_data2)

# corHMM runs
fitcorER <- rayDISC(tree, lh_data3, model="ER", node.states="marginal",state.recon="subsequently")
fitcorARD <- rayDISC(tree, lh_data3, model="ARD", node.states="marginal",state.recon="subsequently")

# Log likelihood test
log.lik.test.ER <- data.frame(fitcorER$states)
log.lik.test.ER$stat <- abs(log(log.lik.test.ER$X0)-log(log.lik.test.ER$X1))
log.lik.test.ER$test <- log.lik.test.ER$stat >= 2
log.lik.test.ER$test <- gsub(TRUE,"*",log.lik.test.ER$test)
log.lik.test.ER$test <- gsub(FALSE,"",log.lik.test.ER$test)

# Set the colors for each set of labels for plotting
cols <- lh_data$devel.color
cols <- setNames(cols, gsub("_.*","",lh_data$genus))
colors <- c("chocolate1","darkslategray3")

# Reordering the vector with the tip colors to match the taxon order in the tree
colsmatch <- match(names(cols),ultratree$tip.label)
cols2 <- cols[order(colsmatch)]

# Phylogeny plot
setEPS()
postscript("ASR_ER_polyclad_genus.eps")
plot.phylo(ultratree, label.offset=0.04, edge.width=2.5, cex=0.8)
tiplabels(pch=22, bg=cols2,cex=1.5)
nodelabels(pie=fitcorER$states,piecol=colors,cex=0.5)
nodelabels(log.lik.test.ER$test, col="purple", bg=NA,frame="none", adj=c(1.7,-0.1), cex=1.2)
dev.off()

# Scaled likelihood data set
scaled_likelihoods <- data.frame(fitcorER$states)
colnames(scaled_likelihoods) <- c("Indirect","Direct")

# Write the data set to a .csv file to be edited in excel
write.csv(scaled_likelihoods,"scaled_likelihoods_genus.csv")

# Save phylogeny with node labels that correspond to the scaled likelihood matrix
setEPS()
postscript("genus_tree_w_node_labels.eps")
plot.phylo(ultratree, label.offset=0.015)
nodelabels(1:30, frame="circle", cex=0.7)
dev.off()

####################
## SPECIES LEVEL ###
####################
# Assign correct tree topology to the variable tree
tree.s <- read.tree("polyclad_37_4469_part_besttree_w_500BS_taxfixed.tre")

# Here I am dropping all outgroup tips to only keep one of each species
tree.s <- drop.tip(tree.s, "Armatoplana_leptalea_B")

# Modifying species label
tree.s$tip.label<-gsub("Armatoplana_leptalea_A","Armatoplana_leptalea",tree.s$tip.label)

# To view tree
plot.phylo(tree.s)

# Create an ultrametric tree with correlated model, output of chronos 
# needed to be rewritten using read.tree, then test for ultrametricity
chronotree.s <- chronos(tree.s, lambda=1, model = "correlated")
ultratree.s <- read.tree(text=write.tree(chronotree.s))
is.ultrametric(ultratree.s)

ultratree.s <- rotate(ultratree.s,40)
ultratree.s <- rotate(ultratree.s,42)
ultratree.s <- rotate(ultratree.s,43)
ultratree.s <- rotate(ultratree.s,44)
ultratree.s <- rotate(ultratree.s,46)
ultratree.s <- rotate(ultratree.s,49)
ultratree.s <- rotate(ultratree.s,50)
ultratree.s <- rotate(ultratree.s,55)
ultratree.s <- rotate(ultratree.s,56)
ultratree.s <- rotate(ultratree.s,58)
ultratree.s <- rotate(ultratree.s,59)
ultratree.s <- rotate(ultratree.s,60)
ultratree.s <- rotate(ultratree.s,61)
ultratree.s <- rotate(ultratree.s,62)
ultratree.s <- rotate(ultratree.s,63)
ultratree.s <- rotate(ultratree.s,64)
ultratree.s <- rotate(ultratree.s,65)

########################################
# Ancestral state reconstruction 
########################################

lh_data.s <- read.csv("life_history_traits_full_polycladida_species.csv",stringsAsFactors=FALSE)

# Change any spaces into underscores for species names, then generate a vector with 
# the appropriate feeding data where the names are set to the species names
lh_data.s$species<-gsub(" ","_",x=lh_data.s$species)
lh_data2.s <- data.frame(lh_data.s$species,as.character(lh_data.s$devel.code))

lh_data2.s <- lh_data2.s$as.character.lh_data.s.devel.code.
lh_data2.s <- setNames(lh_data2.s, lh_data.s$species)
datamatch.s <- match(names(lh_data2.s),ultratree.s$tip.label)
lh_data2.s <- lh_data2.s[order(datamatch.s)]
lh_data3.s <- data.frame(names(lh_data2.s),lh_data2.s)
lh_data3.s$lh_data2.s <- as.character(lh_data3.s$lh_data2.s)

# corHMM runs
fitcorER.s <- rayDISC(tree.s, lh_data3.s, model="ER", node.states="marginal",state.recon="subsequently")
fitcorARD.s <- rayDISC(tree.s, lh_data3.s, model="ARD", node.states="marginal",state.recon="subsequently")

# Log likelihood test
log.lik.test.ER.s <- data.frame(fitcorER.s$states)
log.lik.test.ER.s$stat <- abs(log(log.lik.test.ER.s$X0)-log(log.lik.test.ER.s$X1))
log.lik.test.ER.s$test <- log.lik.test.ER.s$stat >= 2
log.lik.test.ER.s$test <- gsub(TRUE,"*",log.lik.test.ER.s$test)
log.lik.test.ER.s$test <- gsub(FALSE,"",log.lik.test.ER.s$test)

# Set the colors for each set of labels for plotting
cols.s <- lh_data.s$devel.ind
cols.s <- setNames(cols.s, lh_data.s$species)
colors <- c("chocolate1","darkslategray3")

# Reordering the vector with the tip colors to match the taxon order in the tree
colsmatch.s <- match(names(cols.s),ultratree.s$tip.label)
cols2.s <- cols.s[order(colsmatch.s)]

# Phylogeny plot
setEPS()
postscript("ASR_ER_polyclad_species.eps")
plot.phylo(ultratree.s, label.offset=0.04, edge.width=2.5, cex=0.8)
tiplabels(pch=22, bg=cols2.s,cex=1.5)
nodelabels(pie=fitcorER.s$states,piecol=colors,cex=0.5)
nodelabels(log.lik.test.ER.s$test, col="purple", bg=NA,frame="none", adj=c(1.7,-0.1), cex=1.2)
dev.off()

# Scaled likelihood data set
scaled_likelihoods.s <- data.frame(fitcorER.s$states)
colnames(scaled_likelihoods.s) <- c("Indirect","Direct")

# Write the data set to a .csv file to be edited in excel
write.csv(scaled_likelihoods.s,"scaled_likelihoods_species.csv")

# Save phylogeny with node labels that correspond to the scaled likelihood matrix
setEPS()
postscript("species_tree_w_node_labels.eps")
plot.phylo(ultratree.s, label.offset=0.015)
nodelabels(1:34, frame="circle", cex=0.7)
dev.off()

