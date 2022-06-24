####################################
# ancestral_state_platyhelminth.R
# Written by: Jessica A. Goodheart
# Last Updated: 24 June 2022
# Purpose: To perform ancestral state reconstruction on the Platyhelminthes phylogeny
# Inputs used: Platyhelminthes phylogeny (Laumer et al 2015) and life history traits file
####################################

#########################################
# Script setup
#########################################

directory <- "~/[PATH_TO]/platyhelminth_ASR"
setwd(directory)

library(ape)
library(corHMM)

#########################################
# Pull in tree
#########################################

# Assign correct tree topology to the variable tree
tree <- read.tree("platyhelminths_relabeled.tre")

########################################
# Ancestral state reconstruction 
########################################

# Pull the data matrix of interest into R
lh_data <- read.csv("life_history_traits.csv",stringsAsFactors=FALSE)
lh_data2 <- data.frame(lh_data$taxon,lh_data$devel.code)

# corHMM runs
fitcorER <- rayDISC(tree, lh_data2, model="ER", node.states="marginal",state.recon="subsequently")
fitcorARD <- rayDISC(tree, lh_data2, model="ARD", node.states="marginal",state.recon="subsequently")

# Log likelihood test
log.lik.test.ER <- data.frame(fitcorER$states)
log.lik.test.ER$stat <- abs(log(log.lik.test.ER$X0)-log(log.lik.test.ER$X1))
log.lik.test.ER$test <- log.lik.test.ER$stat >= 2
log.lik.test.ER$test <- gsub(TRUE,"*",log.lik.test.ER$test)
log.lik.test.ER$test <- gsub(FALSE,"",log.lik.test.ER$test)

# Set the colors for each set of labels for plotting
cols <- lh_data$devel.color
cols <- setNames(cols, lh_data$taxon)
cols <- gsub("red","chocolate1",cols)
cols <- gsub("blue","darkslategray3",cols)
cols <- gsub("grey","gray48",cols)
colors <- c("chocolate1","darkslategray3")

# Reordering the vector with the tip colors to match the taxon order in the tree
colsmatch <- match(names(cols),tree$tip.label)
cols2 <- cols[order(colsmatch)]

# Phylogeny plot
setEPS()
postscript("ASR_ER_platyhelminth.eps")
plot.phylo(tree, label.offset=0.02, edge.width=2.5, cex=0.8)
tiplabels(pch=22, bg=cols2,cex=1.5)
nodelabels(pie=fitcorER$states,piecol=colors,cex=0.5)
nodelabels(log.lik.test.ER$test, col="purple", bg=NA,frame="none", adj=c(1.7,-0.1), cex=1.2)
dev.off()

# Save phylogeny with node labels that correspond to the scaled likelihood matrix
setEPS()
postscript("Tree_w_node_labels.eps")
plot.phylo(tree, label.offset=0.015)
nodelabels(frame="circle", cex=0.7)
dev.off()

# Scaled likelihood data set
scaled_likelihoods <- data.frame(fitcorER$states)
colnames(scaled_likelihoods) <- c("Indirect","Direct")

# Write the data set to a .csv file to be edited in excel
write.csv(scaled_likelihoods,"scaled_likelihoods.csv")

