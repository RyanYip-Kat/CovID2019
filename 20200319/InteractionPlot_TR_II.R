library(Seurat)
library(pheatmap)
#library(NMF)
library(dplyr)
library(viridis)
#library(hrbrthemes)
#library(igraph)
library(RColorBrewer)
library(stringr)
library(argparse)

set.seed(7777)

print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--slot",
                    type="character",
                    default="data",
                    help="which slot data to be used")


parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}
combn.ident<-function(idents){
  df<-expand.grid(idents,idents)
  idx<-apply(df,1,function(z){
    if(z[1]==z[2]){
      return(FALSE)
    }else{
      return(TRUE)
    }
  })
  df<-df[idx,]
  rownames(df)<-1:nrow(df)
  return(t(df))
}
Seuratinfo<-function(seurat){
  value<-seurat
  n<-ncol(seurat) # number of cells
  c<-colnames(seurat) # Cells of seurat
  g<-rownames(seurat) # Genes of seurat
  return(list("v"=value,"n"=n,"c"=c,"g"=g))
}

InteractionScore<-function(sL,sR,Lg,Rg){
  # L : Ligand
  # R : Receptor
  vL<-sL[["v"]]
  nL<-sL[["n"]]
  Lc<-sL[["c"]]
  
  vR<-sR[["v"]]
  nR<-sR[["n"]]
  Rc<-sR[["c"]]
  
  Le<-vL[Lg,Lc]
  Re<-vR[Rg,Rc]
  
  p1<-sum(Le)
  p2<-sum(Re)
  
  I<-p1*p2/(nL*nR)
  return(I)
}

Score<-function(sL,sR,interaction_csv){
  score<-c()
  row<-c()
  pairs<-read.csv(interaction_csv,header=TRUE)
  colnames(pairs)<-c("Ligand","Receptor")
  n<-nrow(pairs)
  ligand<-as.character(pairs$Ligand)
  receptor<-as.character(pairs$Receptor)
  eps<-1e-6
  
  Lgs<-sL[["g"]]
  Rgs<-sR[["g"]]
  for(i in 1:n){
    Lg=ligand[i]
    Rg=receptor[i]
    if(Lg%in%Lgs & Rg%in%Rgs){
      I<-InteractionScore(sL,sR,Lg,Rg)
      score<-c(score,I)
      name<-paste0(Lg,"-",Rg)
      row<-c(row,name)
    }
  }
  names(score)<-row
  #score<-as.data.frame(score)
  #rownames(score)<-row
  #colnames(score)<-col.name
  return(score)
}

GetSeuratInteractionScore<-function(seurat,idents,interaction_csv,slot){
  if(length(idents)<2){
    stop("Input idents must more than 2")
  }
  eps<-1e-4
  scores<-list()
  
  mat<-GetAssayData(seurat,slot)
  idents.combn<-combn.ident(idents)
  
  n.combn<-ncol(idents.combn)
  for(i in 1:n.combn){
    sL<-mat[,WhichCells(seurat,idents=idents.combn[1,i])]
    sR<-mat[,WhichCells(seurat,idents=idents.combn[2,i])]
    sLinfo<-Seuratinfo(sL)
    sRinfo<-Seuratinfo(sR)
    
    I1<-Score(sLinfo,sRinfo,interaction_csv)
    I2<-Score(sRinfo,sLinfo,interaction_csv)
    
    s<-data.frame(log2((I1+eps)/(I2+eps)))
    colnames(s)<-paste0(idents.combn[1,i],"_",idents.combn[2,i])
    scores[[i]]<-s
  }
  return(scores)
}


interaction_csv<-"interaction.csv"
bc<-readRDS("../20200315/output/BC/model/seurat.rds")
mc<-readRDS("../20200315/output/MC/model/seurat.rds")
nkt<-readRDS("../20200315/output/NKT/model/seurat.rds")

bc$celltype<-"BC"  # BC
Idents(mc)<-mc$UMAP_clusters  # MC
mc<-RenameIdents(mc,"9"="DC","12"="DC","13"="DC","16"="DC",
                     "2"="Monocyte","3"="Monocyte","4"="Monocyte","6"="Monocyte","7"="Monocyte","14"="Monocyte",
                     "5"="Monocyte","15"="Monocyte","17"="Monocyte",
                     "1"="Monocyte","8"="Monocyte","10"="Monocyte","11"="Monocyte")

mc<-subset(mc,idents=c("DC"))
mc$celltype<-Idents(mc)
Idents(nkt)<-nkt$UMAP_clusters
nkt<-RenameIdents(nkt,"4"="CD4","5"="CD4","8"="CD4","9"="CD4",
                     "10"="CD4","12"="CD4","14"="CD4","20"="CD4","21"="CD4",
                     "2"="CD8","6"="CD8","7"="CD8","11"="CD8","13"="CD8",
                     "15"="CD8","16"="CD8","17"="CD8",
                     "1"="NK1","3"="NK1","19"="NK2",
                     "18"="Pro-T")

nkt<-subset(nkt,idents=c("CD4","CD8"))
nkt$celltype<-Idents(nkt)
seurat<-merge(x=bc,y=c(mc,nkt),merge.data=TRUE)
print(table(seurat$celltype))
Idents(seurat)<-seurat$status
seurat<-subset(seurat,idents="TR")
Idents(seurat)<-seurat$celltype

idents<-c("CD4","CD8","BC","DC")
scores<-GetSeuratInteractionScore(seurat = seurat,
                                  idents=idents ,
                                  interaction_csv = interaction_csv,
                                  slot="data")

scores<-do.call(cbind,scores)
saveRDS(scores,file.path(args$outdir,"InteractionScores.rds"))
idx<-apply(scores,1,function(score){
  if(all(score==0)){
    return(TRUE)
  }else{
    return(FALSE)
  }
})
scores<-scores[!idx,]
mat=as.matrix(scores)
cell2cell<-c("DC_CD4","DC_BC","DC_CD8","BC_CD4","CD4_BC")
mat<-mat[,cell2cell]
pdf(file.path(args$outdir,"InteractionScore-1.pdf"),width=8,height=16)
pheatmap(mat[1:84,],
         border_color=NA,
         cluster_rows =F,
         cluster_cols = F,
         fontsize_row = 10,
         angle_col=90,
         fontsize =10)
dev.off()

pdf(file.path(args$outdir,"InteractionScore-2.pdf"),width=8,height=16)
pheatmap(mat[85:nrow(mat),],
           border_color=NA,
           cluster_rows =F,
           cluster_cols = F,
           fontsize_row = 10,
           angle_col=90,
           fontsize =10)
dev.off()

pdf(file.path(args$outdir,"InteractionScore-3.pdf"),width=8,height=32)
pheatmap(mat,
             border_color=NA,
             cluster_rows =F,
             cluster_cols = F,
             fontsize_row = 10,
             angle_col=90,
             fontsize_col=12)
dev.off()
