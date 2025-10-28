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

args <- commandArgs()
library(skmeans)

# --- Deteksi otomatis posisi argumen ---
if (length(args) >= 8) {
  # kalau dipanggil dari Python (TimesVector)
  fname <- args[6]
  k <- as.integer(args[7])
  outdir <- args[8]
  init_file <- if (length(args) >= 9) args[9] else NA
} else {
  # kalau dipanggil langsung (Colab / manual)
  fname <- args[1]
  k <- as.integer(args[2])
  outdir <- args[3]
  init_file <- if (length(args) >= 4) args[4] else NA
}

cat("📁 Expression file:", fname, "\n")
cat("🔢 K =", k, "\n")
cat("💾 Output dir:", outdir, "\n")

# --- Baca data ekspresi ---
gene_expr <- read.table(fname, header = TRUE)
rownames(gene_expr) <- gene_expr[,1]
gene_exprmat <- as.matrix(gene_expr[, -1])

# --- Nama file output ---
outclust <- file.path(outdir, paste0("K", k, ".cluster"))
outproto <- file.path(outdir, paste0("K", k, ".prototype"))

# --- Jika ada file init custom ---
if (!is.na(init_file) && file.exists(init_file)) {
  cat("🟢 Using custom initial centroids from:", init_file, "\n")
  init_cent <- as.matrix(read.table(init_file))
  if (ncol(init_cent) != ncol(gene_exprmat))
    stop("❌ Dimension mismatch between init centroid and expression matrix!")

  cl <- skmeans(gene_exprmat, k, method="genetic", m=1, weights=1,
                control=list(init=init_cent))
  write.table(init_cent,
              file=file.path(outdir, paste0("K", k, ".used_initial_prototype")),
              sep="\t", quote=FALSE, col.names=FALSE)
} else {
  cat("🟡 No init file provided. Using default initialization.\n")
  cl <- skmeans(gene_exprmat, k, method="genetic", m=1, weights=1)
}

# --- Simpan hasil ---
write.table(cl$cluster, file=outclust, sep="\t", quote=FALSE, col.names=FALSE)
write.table(cl$prototypes, file=outproto, sep="\t", quote=FALSE, col.names=FALSE)

cat("\n✅ Done!\n")
cat("   Saved cluster assignment →", outclust, "\n")
cat("   Saved final centroids     →", outproto, "\n")
if (!is.na(init_file) && file.exists(init_file)) {
  cat("   Saved used init centroids →",
      file.path(outdir, paste0("K", k, ".used_initial_prototype")), "\n")
}
