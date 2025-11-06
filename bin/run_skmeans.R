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

#reproducible
set.seed(42)
#

fname=args[6]
k=strtoi(args[7])
outdir=args[8]
gene_expr<-read.table(fname, header=TRUE)
rownames(gene_expr)=gene_expr[,1]
gene_exprmat=as.matrix(gene_expr[,2:ncol(gene_expr)])

# generate initial centers dari data
init_cent <- as.matrix(read.table("/content/centroid.txt"))

# normalisasi (opsional tapi disarankan untuk spherical)
init_cent <- init_cent / sqrt(rowSums(init_cent^2))
#

outclust=paste(outdir, "/K", k, ".cluster", sep="")
outprototype=paste(outdir, "/K", k, ".prototype", sep="")
cl<-skmeans(gene_exprmat, k, method="genetic", m=1, weights=1)
# pakai inisialisasi ini untuk skmeans
cl <- skmeans(gene_exprmat, k,
              method = "pclust",           # bukan genetic, biar nggak diacak
              control = list(start = list(init_cent)))
#

# RAW centroid (belum dinormalisasi)
raw_centroids <- matrix(NA, nrow=k, ncol=ncol(gene_exprmat))
for (i in 1:k) {
  cluster_genes <- gene_exprmat[cl$cluster == i, , drop=FALSE]
  if (nrow(cluster_genes) > 0) {
    raw_centroids[i, ] <- colMeans(cluster_genes)
  }
}
colnames(raw_centroids) <- colnames(gene_exprmat)
#

# simpan ke file kalau mau lihat nanti
write.table(init_cent,
            file=paste(outdir, "/K", k, ".initial_prototype", sep=""),
            sep="\t", quote=FALSE)
#
# === Simpan ke file ===
write.table(raw_centroids,
            file=paste0(outdir, "/K", k, ".raw_centroid"),
            sep="\t", quote=FALSE, row.names=FALSE)
#


write.table(cl$cluster, file=outclust, sep="\t", quote=FALSE)
write.table(cl$prototypes, file=outprototype, sep="\t", quote=FALSE)
