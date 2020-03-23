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
parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="seurat path")

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
  idents.combn<-combn(idents,2)
  
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


interaction_csv<-"/home/ye/Work/R/SingleCell/Su/human/VDJ/sample/interaction.csv"
seurat<-readRDS(args$seurat)
idents<-c("Monocyte","NK","Effector memory CD8",
          "Effector memory CD4","cDC1","Activated memory B")

scores<-GetSeuratInteractionScore(seurat = seurat,
                                  idents=idents ,
                                  interaction_csv = interaction_csv,
                                  slot="data")

scores<-do.call(cbind,scores)
pdf(file.path(args$outdir,"InteractionScore.pdf"))
pheatmap(as.matrix(scores),
         border_color=NA,
         cluster_rows =F, 
         cluster_cols = F,
	 fontsize_row = 3.5,
         angle_col=90,
         fontsize = 7)
dev.off()
