# TimesVector dengan FULL TRACE
args <- commandArgs()
library(skmeans)

#reproducible
set.seed(42)

fname <- args[6]
k <- strtoi(args[7])
outdir <- args[8]
gene_expr <- read.table(fname, header=TRUE)
rownames(gene_expr) <- gene_expr[,1]
gene_exprmat <- as.matrix(gene_expr[,2:ncol(gene_expr)])

cat("=== DATA INPUT ===\n")
cat("Gene expression data:\n")
print(gene_exprmat)

# AUTO-LOAD CUSTOM CENTROID FILE
centroid_file <- "centroid.txt"
if (file.exists(centroid_file)) {
  cat("\n=== LOADING CENTROIDS ===\n")
  cat("Using custom centroids from:", centroid_file, "\n")
  
  init_cent <- as.matrix(read.table(centroid_file, header=FALSE))
  cat("Raw centroids:\n")
  print(init_cent)
  
  # TRACE: Normalize centroids
  cat("\n=== NORMALIZING CENTROIDS ===\n")
  cent_norms <- sqrt(rowSums(init_cent^2))
  cat("Centroid norms:", cent_norms, "\n")
  
  init_cent_norm <- init_cent / cent_norms
  cat("Normalized centroids:\n")
  print(init_cent_norm)
  
  # TRACE: Normalize gene data
  cat("\n=== NORMALIZING GENE DATA ===\n")
  gene_norms <- sqrt(rowSums(gene_exprmat^2))
  cat("Gene norms:\n")
  print(gene_norms)
  
  gene_exprmat_norm <- gene_exprmat / gene_norms
  cat("Normalized gene data:\n")
  print(gene_exprmat_norm)
  
  # TRACE: Calculate cosine distances to initial centroids
  cat("\n=== INITIAL COSINE DISTANCES ===\n")
  for(i in 1:nrow(gene_exprmat_norm)) {
    dist1 <- 1 - sum(gene_exprmat_norm[i,] * init_cent_norm[1,])
    dist2 <- 1 - sum(gene_exprmat_norm[i,] * init_cent_norm[2,])
    cat(rownames(gene_exprmat)[i], 
        "-> Cent1:", round(dist1, 6), 
        "-> Cent2:", round(dist2, 6), 
        "-> Initial Cluster:", ifelse(dist1 < dist2, 1, 2), "\n")
  }
  
  # Run skmeans with custom initialization
  cat("\n=== RUNNING SKMEANS ===\n")
  cl <- skmeans(gene_exprmat_norm, k, control = list(start = list(init_cent_norm)))
  
} else {
  cat("Custom centroid file not found. Using genetic algorithm.\n")
  cl <- skmeans(gene_exprmat, k, method="genetic", m=1, weights=1)
}

cat("\n=== FINAL RESULTS ===\n")
cat("Cluster assignments:\n")
print(cl$cluster)

cat("Final prototypes (centroids):\n")
print(cl$prototypes)

# Calculate raw centroids from original data
cat("\n=== RAW CENTROIDS CALCULATION ===\n")
raw_centroids <- matrix(NA, nrow=k, ncol=ncol(gene_exprmat))
for (i in 1:k) {
  cluster_genes <- gene_exprmat[cl$cluster == i, , drop=FALSE]
  cat("Cluster", i, "has", nrow(cluster_genes), "genes:\n")
  print(rownames(cluster_genes))
  if (nrow(cluster_genes) > 0) {
    raw_centroids[i, ] <- colMeans(cluster_genes)
    cat("Raw centroid", i, ":", round(raw_centroids[i, ], 6), "\n")
  }
}
colnames(raw_centroids) <- colnames(gene_exprmat)

# Output files
outclust <- paste(outdir, "/K", k, ".cluster", sep="")
outprototype <- paste(outdir, "/K", k, ".prototype", sep="")

write.table(cl$cluster, file=outclust, sep="\t", quote=FALSE)
write.table(cl$prototypes, file=outprototype, sep="\t", quote=FALSE)
write.table(raw_centroids, file=paste0(outdir, "/K", k, ".raw_centroids"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat("\n=== FILES SAVED ===\n")
cat("Cluster assignments:", outclust, "\n")
cat("Prototypes:", outprototype, "\n")
cat("Raw centroids:", paste0(outdir, "/K", k, ".raw_centroids"), "\n")
cat("Clustering completed successfully!\n")
