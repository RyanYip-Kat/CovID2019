library(Seurat)
library(argparse)

print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="seurat path")

parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}
seurat<-readRDS(args$seurat)
Idents(seurat)<-seurat$UMAP_clusters
# kangyan
gene1<-c("SOCS2","ARG2","CTLA4","PDCD1","CD274","CD47","SIRPA",
         "CD24","CD52")
gene2<-c("ENTPD1","NT5E","IL10",
         "IL17A","IL37","IL1F10","SIGLEC10","SIGLEC9","SIGLEC6")
gene3<-c(
         "SIGLEC15","CD22","IL1RN","CD33",
         "NLRP12","CD200","CD200R1","IDO1","CD38")
gene4<-c("CD72","FCGR2B","LILRB1","LILRB2","FCRL4","HAVCR2","LAG3",
         "TMEM173","ISG15")



###### cytof_pannel
gene5<-c("FOXP3","HLA-DRA","CD28","ITGAX","NCAM1","CD14","PDCD1","CCR6","CXCR3")
gene6<-c("CD8A","CD4","CD3E","PTPRC","CCR4","PTPRC","CEACAM8","ITGAM","CD19")
gene7<-c("IL2RA","RORC","XCR1","MRC1","B3GAT1","IL3RA","IL7R","CCR7","CLEC10A")
gene8<-c("CLEC4A","TGFB1","SIGLEC6","AXL","CD1C","SELL","FCGR1A","CD163","SEMA4D","CD34")
gene9<-c("SIRPA","CLEC4C","THBD","CLEC9A")



####### zhipu
cd4<-c("CD3E","CD4","CD8A","CXCR3","CCR6","CCR4","PTPRC")
th1<-c("TBX21","IFNG")
th2<-c("GATA3","IL4")
th17<-c("IL17","RORC")
treg<-c("FOXP3","MKI67","IL2RA","KLRB1","CTLA4","PDCD1","IL10")

cd8<-c("GZMB","GZMK")
b_cell<-c("CD19","MS4A1","MZB1","MME")
nk<-c("CD27","NCAM1","B3GAT1","ITGAM")
dc<-c("ITGAX","CD1C","CLEC9A","HLA-DRA","IL3RA")
mono<-c("CD14","FCGR3A")
ilc<-c("IL7R","CCL5")
other<-c("CXCR4","ITGA2B")

jpeg(file.path(args$outdir,"tSNE-Clusters.jpeg"),width = 1024,height = 1024)
DimPlot(seurat,reduction = "tsne",label = TRUE,label.size = 7.5)
dev.off()

jpeg(file.path(args$outdir,"UMAP-Clusters.jpeg"),width = 1024,height = 1024)
DimPlot(seurat,reduction = "umap",label = TRUE,label.size = 7.5)
dev.off()

jpeg(file.path(args$outdir,"cd4.jpeg"),width = 1024,height = 1024)
FeaturePlot(seurat,reduction ="umap",features = cd4)
dev.off()

jpeg(file.path(args$outdir,"th1.jpeg"),width = 1024,height = 1024)
FeaturePlot(seurat,reduction ="umap",features = th1)
dev.off()

jpeg(file.path(args$outdir,"th2.jpeg"),width = 1024,height = 1024)
FeaturePlot(seurat,reduction ="umap",features = th2)
dev.off()

jpeg(file.path(args$outdir,"treg.jpeg"),width = 1024,height = 1024)
FeaturePlot(seurat,reduction ="umap",features = treg)
dev.off()

jpeg(file.path(args$outdir,"cd8.jpeg"),width = 1024,height = 1024)
FeaturePlot(seurat,reduction ="umap",features = cd8)
dev.off()

jpeg(file.path(args$outdir,"B_cell.jpeg"),width = 1024,height = 1024)
FeaturePlot(seurat,reduction ="umap",features = b_cell)
dev.off()

jpeg(file.path(args$outdir,"nk.jpeg"),width = 1024,height = 1024)
FeaturePlot(seurat,reduction ="umap",features = nk)
dev.off()

jpeg(file.path(args$outdir,"dc.jpeg"),width = 1024,height = 1024)
FeaturePlot(seurat,reduction ="umap",features = dc)
dev.off()

jpeg(file.path(args$outdir,"mono.jpeg"),width = 1024,height = 1024)
FeaturePlot(seurat,reduction ="umap",features = mono)
dev.off()

jpeg(file.path(args$outdir,"ilc.jpeg"),width = 1024,height = 1024)
FeaturePlot(seurat,reduction ="umap",features =ilc)
dev.off()

jpeg(file.path(args$outdir,"other.jpeg"),width = 1024,height = 1024)
FeaturePlot(seurat,reduction ="umap",features = other)
dev.off()












