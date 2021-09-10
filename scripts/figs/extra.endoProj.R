rm(list=ls())
source(here::here("scripts/utils.R"))

library(SingleCellExperiment)
library(scuttle)
library(tricycle)
library(ggrepel)
library(scater)

data(neuroRef)
endo.o <- qread(here::here("data/plotdata/endo.qs"))
endo_full.o <- qread(here::here("data/endo.qs"))

metadata(endo_full.o) <- metadata(endo.o)

endo_rotation.m <- attr(reducedDim(endo.o), "rotation")

int_genes.v <- intersect(rownames(endo_full.o), neuroRef$symbol)

set.seed(500)
endo_full.o <- runPCA(endo_full.o, exprs_values = "log.s",  name = "newPCA", ncomponents = 20, subset_row = int_genes.v)

tmp <- attr(reducedDim(endo_full.o, "newPCA"), "percentVar")
reducedDim(endo_full.o, "newPCA")[, 1] <- - reducedDim(endo_full.o, "newPCA")[, 1]
attr(reducedDim(endo_full.o, "newPCA"), "percentVar") <- tmp

#
pca.celltype.p <- plotEmbScat(sce.o = endo.o, dimred = "go.pca", 
															x_lab = str_c("PC1", sprintf(" (%i%%)", round(attr(reducedDim(endo.o, "go.pca"), "percentVar")[1]))),
															y_lab = str_c("PC2", sprintf(" (%i%%)", round(attr(reducedDim(endo.o, "go.pca"), "percentVar")[2]))),
															color_by = "cell_type", facet_by = NULL, facet_labels = NULL, 
															colors.v = metadata(endo.o)$clusters_colors, color.name = "Cell type",
															labels.v = levels(factor(endo.o$cell_type)), title = "PCA of GO (mPancreas)") 

newPCA.celltype.p <- plotEmbScat(sce.o = endo_full.o, dimred = "newPCA", 
						x_lab = str_c("PC1", sprintf(" (%i%%)", round(attr(reducedDim(endo_full.o, "newPCA"), "percentVar")[1]))),
						y_lab = str_c("PC2", sprintf(" (%i%%)", round(attr(reducedDim(endo_full.o, "newPCA"), "percentVar")[2]))),
						color_by = "cell_type", facet_by = NULL, facet_labels = NULL, 
						colors.v = metadata(endo_full.o)$clusters_colors, color.name = "Cell type",
						labels.v = levels(factor(endo_full.o$cell_type)), title = str_c("PCA of ", length(int_genes.v), " genes (mPancreas)"))

mp <- plot_grid(pca.celltype.p + theme(legend.position = "none"), 
								newPCA.celltype.p + theme(legend.position = "none"), 
								get_legend(pca.celltype.p),
								nrow = 1, ncol = 3, rel_widths = c(1, 1, 0.35), label_size = 10, labels = c("a", "b", ""), align = "hv", axis = "tblr")

save_plot(here::here("figs", "response_figs", "extra.endoProj.pdf"), mp,
					base_height = 2 , base_width = 2*1.4, nrow = 1, ncol = 2, device = cairo_pdf)






