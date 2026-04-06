#! /usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

###############################################################
# FULL R PIPELINE FOR TE CLUSTERING / NRI-NTI / DISTANCE TESTS
# Tree file: Jeff_Lsyl_TAR_modHeads.fasttree
# Tip labels contain species in the 2nd "|" field: TE11|Lsyl|TAR|NONE
###############################################################

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(ape)
  library(picante)
  library(phytools)
  library(ggtree)
  library(ggplot2)
})

# ---- Load tree ----
# tree_file <- "JeffLsyl/Jeff_Lsyl_TAR_modHeads.fasttree"
tree_file <- args[1]
te_family <- strsplit(args[1],"_")[[1]][3]
tree <- read.tree(tree_file)

cat("Loaded tree with", length(tree$tip.label), "tips.\n")

# ---- Parse species names from tip labels ----
# Tip format example: TE11|Lsyl|TAR|NONE
split_labels <- strsplit(tree$tip.label, "\\|")
species <- sapply(split_labels, function(x) x[2])   # second field = species
tipnames <- tree$tip.label

# Create a mapping table
tip_map <- data.frame(tip = tipnames, species = species, stringsAsFactors = FALSE)
species_list <- unique(species)
cat("Species present:", paste(species_list, collapse=", "), "\n")

# ---- Compute patristic distances ----
pat <- cophenetic(tree)

# ---- Compute within-species mean distances ----
within_mean <- sapply(species_list, function(sp) {
  pts <- tip_map$tip[tip_map$species == sp]
  m <- pat[pts, pts]
  mean(m[lower.tri(m)])
})

cat("\nWithin-species mean patristic distances:\n")
print(within_mean)

# ---- Compute between-species mean distance ----
if(length(species_list)==2){
  s1 <- species_list[1]
  s2 <- species_list[2]
  m2 <- pat[
    tip_map$tip[tip_map$species == s1],
    tip_map$tip[tip_map$species == s2]
  ]
  between_mean <- mean(m2)
  cat("\nBetween-species mean patristic distance:", between_mean, "\n")
}

# ---- Build community matrix for NRI/NTI ----
comm <- matrix(
  0,
  nrow = length(species_list),
  ncol = length(tipnames),
  dimnames = list(species_list, tipnames)
)
for(i in seq_len(nrow(tip_map))){
  comm[ tip_map$species[i], tip_map$tip[i] ] <- 1
}

# ---- MPD/MNTD with null model (999 permutations) ----
cat("\nRunning ses.mpd (NRI) and ses.mntd (NTI)...\n")
mpd_out <- ses.mpd(comm, cophenetic(tree), null.model="taxa.labels", runs=999)
mntd_out <- ses.mntd(comm, cophenetic(tree), null.model="taxa.labels", runs=999)

cat("\n=== MPD / NRI results ===\n")
print(mpd_out)

cat("\n=== MNTD / NTI results ===\n")
print(mntd_out)

# Interpretation:
# NRI  = -1 * mpd_out$mpd.obs.z
# NTI  = -1 * mntd_out$mntd.obs.z
# Positive and significant -> clustering (few founder lineages)

# ---- Simple permutation test for difference of within-species distances ----
cat("\nRunning permutation test for difference in within-species distances...\n")

nperm <- 1000
obs_diff <- within_mean[1] - within_mean[2]
perm_diffs <- numeric(nperm)

for(i in 1:nperm){
  perm_sp <- sample(tip_map$species)
  perm_w <- sapply(species_list, function(sp){
    pts <- tip_map$tip[perm_sp == sp]
    m <- pat[pts, pts]
    mean(m[lower.tri(m)])
  })
  perm_diffs[i] <- perm_w[1] - perm_w[2]
}

pval <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (nperm + 1)
cat("Observed difference:", obs_diff, "\n")
cat("Permutation P-value:", pval, "\n")

# ---- Estimate number of founders: cluster by patristic distance ----
cat("\nClustering tree by patristic distances to infer founders...\n")
dmat <- as.dist(pat)
hc <- hclust(dmat, method="average")

# choose cut height based on ~10% of max distance (adjustable)
cut_h <- 0.1 * max(dmat)
clusters <- cutree(hc, h = cut_h)

cat("Clusters (founder groups) using height =", cut_h, "\n")
print(table(clusters, tip_map$species))

# ---- Plot tree with species-colored tips ----
pdf(paste0("Jeff_Lsyl_", te_family, "_tree_colored.pdf"), width=10, height=14)

p <- ggtree(tree, layout="rectangular") %<+% tip_map +
  geom_tippoint(aes(color=species), size=2) +
  scale_color_manual(values=c("Jeff"="gray","Lsyl"="#0000fff5")) +
  theme_tree2() +
  ggtitle(paste0(te_family, " TE phylogeny — Jeff vs Lsyl"))

print(p)
dev.off()

cat("\nSaved colored tree to Jeff_Lsyl_TE_tree_colored.pdf\n")

###############################################################
# END OF SCRIPT
###############################################################
