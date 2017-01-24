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
gene_expr<-read.table(fname, header=TRUE)

rownames(gene_expr)=gene_expr[,1]
gene_exprmat=as.matrix(gene_expr[,2:ncol(gene_expr)])

outclust=paste(outdir, "/K", k, ".cluster", sep="")
outprototype=paste(outdir, "/K", k, ".prototype", sep="")
cl<-skmeans(gene_exprmat, k, method="genetic", m=1, weights=1)
write.table(cl$cluster, file=outclust, sep="\t", quote=FALSE)
write.table(cl$prototypes, file=outprototype, sep="\t", quote=FALSE)
