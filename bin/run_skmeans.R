#TimesVector
#
#Created by Inuk Jung on 2016-08-3.
#Copyright (c) 2016 Inuk Jung. All rights reserved.
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

args<-commandArgs()
library(skmeans)
fname=args[6]
k=strtoi(args[7])
outdir=args[8]
#opsional
init_file <- if (length(args) >= 9) args[9] else NA
#
  
gene_expr<-read.table(fname, header=TRUE)
rownames(gene_expr)=gene_expr[,1]
gene_exprmat=as.matrix(gene_expr[,2:ncol(gene_expr)])

outclust=paste(outdir, "/K", k, ".cluster", sep="")
outprototype=paste(outdir, "/K", k, ".prototype", sep="")

#generate intial centers
init_cent <- gene_exprmat[sample(1:nrow(gene_exprmat), k), ]
cl<-skmeans(gene_exprmat, k, method="genetic", m=1, weights=1)
# pakai inisialisasi ini untuk skmeans
cl <- skmeans(gene_exprmat, k, method="genetic", m=1, weights=1,
              control=list(init=init_cent))

# --- Cek apakah ada file inisialisasi custom ---
if (!is.na(init_file) && file.exists(init_file)) {
  cat("🟢 Using custom initial centroids from:", init_file, "\n")
  init_cent <- as.matrix(read.table(init_file))
  
  if (ncol(init_cent) != ncol(gene_exprmat)) {
    stop("❌ Dimension mismatch between init centroid and expression matrix!")
  }
  
  # Jalankan spherical k-means dengan init custom
  cl <- skmeans(gene_exprmat, k, method = "genetic", m = 1, weights = 1,
                control = list(init = init_cent))
  
  # Simpan inisialisasi yang dipakai
  write.table(init_cent,
              file = file.path(outdir, paste0("K", k, ".used_initial_prototype")),
              sep = "\t", quote = FALSE, col.names = FALSE)
} else {
  cat("🟡 No custom init file provided — using default random/genetic initialization.\n")
  cl <- skmeans(gene_exprmat, k, method = "genetic", m = 1, weights = 1)
}
#

# simpan ke file kalau mau lihat nanti
write.table(init_cent,
            file=paste(outdir, "/K", k, ".initial_prototype", sep=""),
            sep="\t", quote=FALSE)
#

write.table(cl$cluster, file=outclust, sep="\t", quote=FALSE)
write.table(cl$prototypes, file=outprototype, sep="\t", quote=FALSE)

#
cat("✅ Done! Saved to:\n")
cat("   ", outclust, "\n")
cat("   ", outproto, "\n")
if (!is.na(init_file) && file.exists(init_file)) {
  cat("   ", file.path(outdir, paste0("K", k, ".used_initial_prototype")), "\n")
}
#
