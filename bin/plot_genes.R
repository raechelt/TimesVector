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

library(ggplot2)

args<-commandArgs()

data_f=args[6]
outdir=strsplit(data_f, "_plot.dat")[[1]]
basedir=args[7]

# reading parameter configuration file
param_fname<-paste(basedir, "/param.conf", sep="")
param_f<-file(param_fname, open="r")
lines<-readLines(param_f)
for (l in lines) {
	tok=strsplit(l, "=")
	if (tok[[1]][1]=="class") {
		class=unlist(strsplit(tok[[1]][2], ","))
	} else if (tok[[1]][1]=="timepoints") {
		tps=unlist(strsplit(tok[[1]][2], ","))
	}
}

# reading data
data<-read.table(data_f, header=T)
clids<-unique(data$Clid)

# reordering factors for class and timepoints
data$Samples<-factor(data$Samples, levels=class)
data$TimePoint<-factor(data$TimePoint, levels=tps)

for (clid in clids) {

	data_clid<-data[data$Clid==clid,]
	cluster_type<-data_clid$Type
	gid_count<-length(unique(data_clid$GeneID))/length(class)
	
	pdf(paste(outdir, "/cluster_", clid, "_genes.pdf", sep=""), bg="transparent", height=5, width=6)
	title=paste("Cluster ", clid, " (",cluster_type, ", ", gid_count, " genes)", sep="")

	p<-ggplot(data=data_clid, aes(x=TimePoint, y=Value, group=GeneID, color=Samples))+
	geom_line(size=0.4, alpha=0.8)+
	scale_shape_discrete(solid=T) +
	geom_point(aes(x=TimePoint, y=Value, shape=Samples), size=3)+
	ggtitle(title)+
	ylab("Expression value")+
	xlab("Time points")+
	theme_bw()+ 
	theme(panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border=element_blank())+
	theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))+
	theme(axis.title.x=element_text(vjust=0))+
	theme(axis.title.y=element_text(vjust=1))+
	scale_x_discrete(expand=c(0.05, 0.05))+
	theme(legend.key = element_blank())

	print(p)

	dev.off()
}





