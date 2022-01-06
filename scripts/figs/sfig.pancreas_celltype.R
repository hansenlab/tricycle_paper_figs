rm(list=ls())
source(here::here("scripts/utils.R"))


endo.o <- qread(here::here("data/plotdata/endo.qs"))
metadata(endo.o)$point.size <- 4.01

### umap figure
tmp.df <- data.frame(pc1.s = reducedDim(endo.o, "tricycleEmbedding")[, 1], pc2.s = reducedDim(endo.o, "tricycleEmbedding")[, 2],
										 cell_type = endo.o$cell_type)

emb.scat.lp <- plotEmbScat(endo.o, dimred = "tricycleEmbedding", 
													 color_by = "cell_type", 
													 colors.v = metadata(endo.o)$clusters_colors,
													 labels.v = levels(endo.o$cell_type), color.name = "Cell type",
													 facet_by = "cell_type", facet_labels = levels(endo.o$cell_type))
emb.scat.lp <- lapply(emb.scat.lp, function(x) x + theme(legend.position = "none"))
emb.scat.lp[[1]] <- emb.scat.lp[[1]] +
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
	theme(legend.position = c(0, 1),
				legend.justification = c(0, 1),
				legend.key.size = unit(6, "pt"))
	
mp <- plot_grid(plotlist = emb.scat.lp,
								nrow = 3, ncol = 3, label_size = 10, labels = "auto", align = "v", axis = "lr")	
# plotThetaDen(sce.o = endo.o, color.var = "cell_type", color.name = "Cell type", colors.v = metadata(endo.o)$clusters_colors,
# 						 labels.v = levels(endo.o$cell_type), bw = 10)$linear


save_plot(here::here("figs", "sfigs", "sfig.pancreas_celltype.pdf"), mp,
					base_height = 2, base_width = 2 * 1.2, nrow = 3, ncol = 3, device = cairo_pdf)

