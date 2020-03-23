library(Seurat)
library(ggplot2)
library(ggpubr)
library(stringr)
library(argparse)

print("***Configure Parameters***")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the output path")


parser$add_argument("--markers",
                    type="character",
                    nargs="+",
                    default=NULL,
                    help="markers to be featureplot")

args<-parser$parse_args()

seurat<-readRDS(args$seurat)
genes<-as.character(args$markers)

pdf(file.path(args$outdir,paste0("violin_",paste(genes,collapse="_"),".pdf")),width=1024,height=1024)
VlnPlot(seurat,features=genes,group.by="orig.Cluster",ncol=1) #+stat_compare_means(aes(label = paste0("p = ", ..p.format..)))
dev.off()


plot.list=list()
for(i in 1:length(genes)){
	exprmat<-GetAssayData(seurat,slot = "data")[genes[i],]
        exprmat<-data.frame("expr"=exprmat,"group"=as.factor(seurat$orig.Cluster))
	p<-ggplot(data=exprmat,aes(x=group,y=expr,fill=group))+
		geom_boxplot()+geom_jitter(width=0.1,alpha=0.2,size=0.5)+
		stat_compare_means(aes(label = paste0("p = ", ..p.format..)))+
		theme_bw()+NoGrid()+ggtitle(paste("(",genes[i],")"))
	plot.list[[i]]=p
}

pdf(file.path(args$outdir,paste0("boxplot_",paste(genes,collapse="_"),".pdf")),width=1024,height=1024)
print(CombinePlots(plot.list))
dev.off()


