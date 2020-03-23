library(Seurat)
library(argparse)
library(monocle3)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

args <- parser$parse_args()
dataset<-args$outdir
model.dir<-file.path(dataset,"model")
plot.dir<-file.path(dataset,"plot")
if(!dir.exists(model.dir)){
        dir.create(model.dir,recursive=TRUE)
}

if(!dir.exists(plot.dir)){
        dir.create(plot.dir,recursive=TRUE)
}



path="/home/ye/Data/Zoc/Cell/20200315_Liu/outs"

counts=Read10X(file.path(path,"filtered_feature_bc_matrix"))
cluster<-read.csv(file.path(path,"analysis/clustering/graphclust/clusters.csv"))
rownames(cluster)<-cluster$Barcode
cluster<-subset(cluster,select=Cluster) 
cell=rownames(cluster)

counts<-counts[,cell]
rownames(cluster)=colnames(counts)
print(paste0("Size of counts  [ ",nrow(counts),",",ncol(counts)," ]"))
print("### Create Seurat object")
object<-CreateSeuratObject(counts= counts,
                       assay = "RNA",
                       project ="scRNA",
                       names.delim="_",
                       min.cells=0, 
                       min.features=0)

object<-AddMetaData(object,metadata=cluster,col.name="orig.Cluster")
cells<-colnames(object)
ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
object$ident<-ident
#status<-ifelse(ident%in%c(1,2,3,4),"Normal",ifelse(ident%in%c(12,9,6,8),"LT","ST"))
#object$status<-status

object[["percent.mt"]] <- PercentageFeatureSet(object,pattern = "^MT-")
object[["percent.rpl"]] <- PercentageFeatureSet(object,pattern = "^RPL")
object[["percent.rps"]] <- PercentageFeatureSet(object,pattern = "^RPS")


pDC<-read.csv("PDC.csv")
pDC.Barcode<-as.character(pDC$Barcode)

cells<-colnames(object)
meta.data<-object@meta.data
meta.data$cells<-cells

new.Orig.Cluster<-with(meta.data,ifelse(cells%in%pDC.Barcode,"pDC",orig.Cluster))
print(table(new.Orig.Cluster))
print(table(meta.data$orig.Cluster))

object$new.Orig.Cluster<-new.Orig.Cluster

mat<-object@meta.data[,"new.Orig.Cluster",drop=FALSE]
mat$barcode=rownames(mat)
mat<-mat[,c(2,1)]
write.table(mat,file.path(args$outdir,"new.Orig.clusters.csv"),sep=",",quote=FALSE,row.names=FALSE)
