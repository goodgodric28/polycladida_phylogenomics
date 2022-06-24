####################################
# rename_tips.R
# Written by: Jessica A. Goodheart
# Last Updated: 24 June 2022
# Purpose: Renames the tips of the Platyhelminthes phylogeny from Laumer et al 2015
# Inputs used: Platyhelminthes phylogeny (Laumer et al 2015) 
####################################

#########################################
# Script setup
#########################################

directory <- "~/[PATH_TO]/platyhelminth_ASR"
setwd(directory)

require(ape)

#########################################
# Pull in and revise tree
#########################################
# Assign correct tree topology to the variable tree
tree <- read.nexus("eLife_platy_BMGE_LG4F_bsbest.tre")

# This roots the tree at Pleurobranchaea californica
tree <- reroot(tree, node.number=76)

# Here I am dropping all outgroup tips to only keep the ingroup taxa
tree <- drop.tip(tree, c("Helo",
                         "Lott",
                         "Capi",
                         "Gnat",
                         "Lepa",
                         "Bpli",
                         "Lepi"))

plot.phylo(tree)

# Create an ultrametric tree with correlated model, output of chronos 
# needed to be rewritten using read.tree, then test for ultrametricity
chronotree <- chronos(tree, lambda=1, model = "correlated")
ultratree <- read.tree(text=write.tree(chronotree))
is.ultrametric(ultratree)



# Here I am dropping all extraneous tips where there are multiple taxa per group
ultratree <- drop.tip(ultratree, c("Mlin",
                                   "Xeno",
                                   "ProI",
                                   "Palp",
                                   "Lepa",
                                   "Stym",
                                   "Sell",
                                   "Mdal",
                                   "Psph",
                                   "Palp",
                                   "Slpr",
                                   "Lehy",
                                   "Dlac",
                                   "Bcan",
                                   "Nmel",
                                   "Ecmu",
                                   "Hyme",
                                   "Csin",
                                   "Sman",
                                   "Pris",
                                   "Bbal"))

# Renaming of tip labels so they are correct
ultratree$tip.label <- gsub("Sleu","Catenulida",ultratree$tip.label)
ultratree$tip.label <- gsub("Mrue","Macrostomorpha",ultratree$tip.label)
ultratree$tip.label <- gsub("Gapp","Prorhynchida",ultratree$tip.label)
ultratree$tip.label <- gsub("Pvit","Polycladida",ultratree$tip.label)
ultratree$tip.label <- gsub("Gnos","Gnosonosemida",ultratree$tip.label)
ultratree$tip.label <- gsub("Rros","Rhabdocoela",ultratree$tip.label)
ultratree$tip.label <- gsub("Feca","Fecampiida",ultratree$tip.label)
ultratree$tip.label <- gsub("Prot","Prolecithophora",ultratree$tip.label)
ultratree$tip.label <- gsub("Smed","Tricladida",ultratree$tip.label)
ultratree$tip.label <- gsub("Gsal","Monogenea",ultratree$tip.label)
ultratree$tip.label <- gsub("Tsol","Cestoda",ultratree$tip.label)
ultratree$tip.label <- gsub("Fmag","Trematoda",ultratree$tip.label)
ultratree$tip.label <- gsub("Bsem","Bothrioplanida",ultratree$tip.label)
ultratree$tip.label <- gsub("Mfus","Proseriata",ultratree$tip.label)

plot.phylo(ultratree)

# Write tree to file
write.tree(ultratree, file="platyhelminths_relabeled.tre")
