library(stringr)
library(dplyr)
library(argparse)

print("*** Configure Parameters ***")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--marker",
                    type="character",
                    default="",
                    help="the marker  to be used")


parser$add_argument("--outdir",
                      type="character",
                      default="",
                      help="save path")

args<-parser$parse_args()
if(!dir.exists(args$outdir)){
            dir.create(args$outdir,recursive=TRUE)
}

markers<-readRDS(args$marker)
topn<-markers%>%group_by(cluster)%>%top_n(70,wt=avg_logFC)
write.table(topn,file.path(args$outdir,"markers.csv"),sep=",",quote=FALSE)

