rm(list=ls())
source(here::here("scripts/utils.R"))

endo.o <- qread(here::here("data/plotdata/endo.qs"))
endo_full.o <- qread(here::here("data/endo.qs"))

metadata(endo_full.o) <- metadata(endo.o)
### umap figure
tmp.df <- data.frame(pc1.s = reducedDim(endo.o, "umap")[, 1], pc2.s = reducedDim(endo.o, "umap")[, 2],
										 cell_type = endo.o$cell_type)

umap.scat.p <- ggplot(tmp.df, aes(x = pc1.s, y = pc2.s, color = cell_type)) +
	geom_scattermore(pointsize = metadata(endo.o)$point.size, alpha = metadata(endo.o)$point.alpha) +
	scale_color_manual(values = metadata(endo.o)$clusters_colors, name = "Cell type", labels = levels(factor(endo.o$cell_type))) +
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
	labs( y = "UMAP-2", x = "UMAP-1", title = "UMAP (mPancreas)") +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8),
				legend.title = element_text(size = 6),
				legend.text = element_text(size = 6))


### GO PC
pca.celltype.p <- plotEmbScat(sce.o = endo.o, dimred = "go.pca", 
																	 x_lab = str_c("PC1", sprintf(" (%i%%)", round(attr(reducedDim(endo.o, "go.pca"), "percentVar")[1]))),
																	 y_lab = str_c("PC2", sprintf(" (%i%%)", round(attr(reducedDim(endo.o, "go.pca"), "percentVar")[2]))),
																	 color_by = "cell_type", facet_by = NULL, facet_labels = NULL, 
																	 colors.v = metadata(endo.o)$clusters_colors, color.name = "Cell type",
																	 labels.v = levels(factor(endo.o$cell_type)), title = "PCA of GO (mPancreas)") + 
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8),
				legend.title = element_text(size = 6),
				legend.text = element_text(size = 6))





### ductal
GO_Ductal.o <- getGO(endo_full.o[, endo_full.o$cell_type == "Ductal"], row.id = rownames(endo_full.o), id.type = "SYMBOL", runSeuratBy =  NULL)
plot(reducedDim(GO_Ductal.o, "PCA.s"))
reducedDim(GO_Ductal.o, "PCA.s")[, 2] <- - reducedDim(GO_Ductal.o, "PCA.s")[, 2]
metadata(GO_Ductal.o)$point.size <- 4.5

ductal.pca.p <- plotScatCC(sce.o = GO_Ductal.o, dimred = "PCA.s", 
					 x_lab = str_c("PC1", sprintf(" (%i%%)", round(attr(reducedDim(GO_Ductal.o, "PCA.s"), "percentVar")[1]))),
					 y_lab = str_c("PC2", sprintf(" (%i%%)", round(attr(reducedDim(GO_Ductal.o, "PCA.s"), "percentVar")[2]))),
					 title = str_c("PCA of Ductal (n=", ncol(GO_Ductal.o), ")")) +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8),
				legend.title = element_text(size = 6),
				legend.text = element_text(size = 6))

row1.mp <- plot_grid(umap.scat.p + theme(legend.position = "none"), 
										 get_legend(umap.scat.p),
										 pca.celltype.p + theme(legend.position = "none"), 
										 nrow = 1, ncol = 3, rel_widths = c(1, 0.5, 1), label_size = 10, labels = c("a", "", "b"), align = "hv", axis = "tblr")
row2.mp <- plot_grid(ggplot() + theme_nothing(), 
										 ductal.pca.p + theme(legend.position = "none"), 
										 get_legend(ductal.pca.p),
										 ggplot() + theme_nothing(), 
										 nrow = 1, ncol = 4, rel_widths = c(0.5, 1, 0.5, 0.5), label_size = 10, labels = c("", "c", ""), align = "hv", axis = "tblr")

mp <- plot_grid(row1.mp, row2.mp,
								nrow = 2, ncol = 1, label_size = 10, labels = NULL, align = "hv", axis = "tblr")

save_plot(here::here("figs", "main", "main.pancreas.pdf"), mp,
					base_height = 2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2, ncol = 2.2, device = cairo_pdf)



