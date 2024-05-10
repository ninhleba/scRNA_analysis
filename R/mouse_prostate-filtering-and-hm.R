#!/usr/bin/env Rscript

library(argparser)
library(Seurat)
library(scMCA)
library(ggplot2)
library(magrittr)

### User inputs and parameters ###

parser<-arg_parser(name="scMCA_heatmap_cluster.R",description="scMCA heatmaps for Seurat clusters")

parser<-add_argument(
  parser,
  arg='--output_dir',
  short = '-o',
  type="character",
  default='./',
  help="Enter output directory path. Format: path/to/output/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--file_prefix',
  short = '-f',
  type="character",
  default='',
  help="Enter file prefix to mark files for current run")

parser<-add_argument(
  parser,
  arg='--top_celltypes_plotted',
  short = '-t',
  type="numeric",
  default=3,
  help= "Top cell types for each cell to plot on heatmaps.")

args <- parse_args(parser)

#Print run parameters
cat("List of arguments:",'\n')
args
cat('\n')

setwd(args$output_dir)

# Prepare output directory
if(! dir.exists(paste0(args$output_dir,'scMCA_filtering/')))
{
  dir.create(paste0(args$output_dir,'scMCA_filtering/'))
}

z <- readRDS(paste0(args$file_prefix,"_clustered.rds"))
DefaultAssay(z) <- 'RNA'

### Compartment Dotplot Generating ####

#mouse genes
uro<-'Krt5,Krt7,Krt8,Krt13,Krt14,Krt17,Krt19,Krt20,Upk1a,Upk1b,Upk2,Upk3a,Upk3b,Cldn4'
stromal<-'Dcn,Acta2,Tagln,Col1a1,Flt1,Pecam1,Kit,Vim,Cd34,Pdgfra'
immune<-'Ptprc,Cd4,Cd8a,Cd8b1,Trac,Nkg7,Klrd1,Cd79b,Lst1,Cd14,Tpsab1'#No KLRF1
prolif<-'Top2a,Mki67,Shh,Trp53'

genes<-list(uro,stromal,immune,prolif)
final_genes<-lapply(genes,strsplit,split=',') %>% unlist %>% unlist()

dp <- DotPlot(z, assay="RNA", features=final_genes, col.min = 0)
ggsave(paste0(args$file_prefix,'_compartment_dp.png'), width = 25, height = 13, units = 'in', bg = 'white', dp)

### Generate scMCA heatmaps and filter prostate cells identified in the process ###

clusters<-levels(z)

Prostate_cells <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(Prostate_cells) <- c('cell_id', 'cluster', 'cell_type','cor.val')

for (i in clusters) {
   cluster_name = i
   c <- subset(z, cells = colnames(z)[z@meta.data[, "seurat_clusters"] == cluster_name])
   DefaultAssay(c) <- 'RNA'
   c <- NormalizeData(c)
   exprMat <- c@assays$RNA@data
   dge <- read.table(paste0("diff_genes/",args$file_prefix,"_cluster",cluster_name,"_diff_markers.tsv"), header = TRUE)
   dge <- dge$gene
   exprMat_dge <- exprMat[dge, ]
   mca_result <- scMCA(scdata = exprMat_dge, numbers_plot = args$top_celltypes_plotted)
   p <- plotMCA(mca_result, numbers_plot = args$top_celltypes_plotted, show_col = FALSE)
   ggsave(paste0('scMCA_filtering/',args$file_prefix,'_cluster',cluster_name,'_scmca_hm.png'),dpi='retina',p)
   top3_celltype_results <- mca_result[["scMCA_probility"]]
   for (cell_id in unique(top3_celltype_results$Cell)) {
        celltypes <- subset(top3_celltype_results, Cell == cell_id)
        highest_cor <- celltypes[which.max(celltypes$Score), ]
        if (grepl("(Prostate)",highest_cor$`Cell type`)) {
        Prostate_cells[nrow(Prostate_cells) + 1, ] = c(as.character(highest_cor$Cell), i, as.character(highest_cor$'Cell type'), highest_cor$Score)
        }
      }
  }

write.csv(Prostate_cells, file = paste0("scMCA_filtering/",args$file_prefix,"_Prostate_cells.csv"))

if (nrow(Prostate_cells) == 0) {
  message("No prostate cells were found. Exiting program.")
  q() # exit R
}

#Prostate_cells <- read.csv(paste0("scMCA_filtering/",args$file_prefix,"_Prostate_cells.csv"), header = TRUE, row.names = 1)

Prostate_cells$cluster <- as.factor(Prostate_cells$cluster)
b <- ggplot(data.frame(Prostate_cells$cluster), aes(x=Prostate_cells$cluster)) + geom_bar() 
b <- b + labs(x="Cluster",y="# Cells",title ="Prostate cells across Seurat clusters identified by scMCA", subtitle = paste0("Total prostate cells: ",nrow(Prostate_cells)))
ggsave(paste0("scMCA_filtering/",args$file_prefix,'_prostate_cells_dist.png'),dpi='retina',b)

cluster_size <- read.csv(paste0(args$file_prefix,"_cells_by_cluster_by_sample.csv"), header = TRUE, row.names = 1) 
prostate_proportion <- data.frame(table(Prostate_cells$cluster)/cluster_size[names(table(Prostate_cells$cluster)),'Total'])
colnames(prostate_proportion) <- c('cluster','prostate_proportion')
o <- ggplot(prostate_proportion, aes(x=cluster, y=prostate_proportion)) + geom_col()
o <- o + labs(x = 'Cluster', y = 'Prostate proportion', title = 'Prostate proportion across Seurat clusters', subtitle = paste0("Total prostate proportion: ",nrow(Prostate_cells)/cluster_size['Total','Total']))
ggsave(paste0("scMCA_filtering/",args$file_prefix,'_prostate_proportion_by_cluster.png'),dpi='retina',o)

rm(exprMat)
rm(exprMat_dge)
rm(mca_result)

pc <- Prostate_cells$cell_id
filtered_cells =  colnames(z)[-which(colnames(z) %in% pc)]
filtered_seurat <- subset(z, cells = filtered_cells)
saveRDS(filtered_seurat, file=paste0(args$file_prefix,'_prostatefilteredscMCA.rds'))

rm(filtered_seurat)
rm(z)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
