library(Seurat)
library(argparse)
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



path="/home/ye/Data/Zoc/Cell/CovID/outs"

counts=Read10X(file.path(path,"filtered_feature_bc_matrix"))
cluster<-read.csv(file.path(path,"analysis/clustering/graphclust/clusters.csv"))
rownames(cluster)<-cluster$Barcode
cluster<-subset(cluster,select=Cluster) 
celltype<-ifelse(cluster$Cluster%in%c(1,8,11,13,15),"MC",ifelse(cluster$Cluster%in%c(2,3,4,5,6,16,18,20),"NKT",ifelse(cluster$Cluster%in%c(7,12,17,19,22,23),"BC","Others")))
cell=rownames(cluster)

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
object$celltype<-celltype
print(table(object$celltype))
cells<-colnames(object)
ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
object$ident<-ident
status<-ifelse(ident%in%c(1,2,3,4),"Normal",ifelse(ident%in%c(12,9,6,8),"LT","ST"))
object$status<-status
object[["percent.mt"]] <- PercentageFeatureSet(object,pattern = "^MT-")
object[["percent.rpl"]] <- PercentageFeatureSet(object,pattern = "^RPL")
object[["percent.rps"]] <- PercentageFeatureSet(object,pattern = "^RPS")
object <- FindVariableFeatures(object, selection.method = "vst",
                            nfeatures = 5000,verbose = FALSE)

genes<-rownames(object)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]

print("### Normalize Data")
object<-NormalizeData(object,normalization.method = "LogNormalize",verbose = FALSE)

print("### Scale Data")
object<-ScaleData(object,features=VariableFeatures(object),model.use = "linear",
               vars.to.regress = c("nFeature_RNA"),verbose =FALSE)

saveRDS(object,file.path(model.dir,"seurat.rds"))
print("### Difference Analysis and Find Markers")
Idents(object)<-object$celltype
markers <- FindAllMarkers(object, only.pos = FALSE,
                          features = keep_genes,
                          test.use = "wilcox",
                          min.pct = 0.2,
                          pseudocount.use = 1 )



print("### Save")
saveRDS(markers,file.path(model.dir,"markers.rds"))




