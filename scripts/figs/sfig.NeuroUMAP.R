rm(list=ls())
source(here::here("scripts/utils.R"))

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))

### plot UMAP
umap.sam.scat.p <- plotEmbScat(sce.o = neurosphere.o, dimred = "umap", x_lab = "UMAP-1", y_lab = "UMAP-2",
															 color_by = "sample", colors.v = c("#729ECE", "#FF9E4A"),
															 labels.v = c("WT1", "WT2"), color.name = 'Sample')
umap.celltype.scat.p <- plotEmbScat(sce.o = neurosphere.o, dimred = "umap",  x_lab = "UMAP-1", y_lab = "UMAP-2",
															 color_by = "cell_type", colors.v = c("#66C2A5", "#8DA0CB",  "#FC8D62", "#A6D854", "#FFD92F", "#B3B3B3"),
															 labels.v = c("aNSCs", "Astrocytes", "NPCs", "Oligodendrocytes", "qNSCs", "Other"), color.name = 'Cell type')
umap.cc.scat.p <- plotEmbScat(sce.o = neurosphere.o, dimred = "umap",  x_lab = "UMAP-1", y_lab = "UMAP-2",
																		color_by = "CCStage", colors.v = ccColors.v,
																		labels.v = ccLabels.v, color.name = 'CC Stage')
umap.cyclone.scat.p <- plotEmbScat(sce.o = neurosphere.o, dimred = "umap",  x_lab = "UMAP-1", y_lab = "UMAP-2",
															color_by = "cyclone", colors.v = cc3Colors.v,
															labels.v = cc3Labels.v, color.name = 'Cyclone')
umap.seurat.scat.p <- plotEmbScat(sce.o = neurosphere.o, dimred = "umap",  x_lab = "UMAP-1", y_lab = "UMAP-2",
																	 color_by = "SeuratCC", colors.v = cc3Colors.v,
																	 labels.v = cc3Labels.v, color.name = 'SeuratCC')

neurosphere.o$log2umi <- log2(neurosphere.o$TotalUMIs)
umap.umi.scat.p <- plotEmbScatViridis(sce.o = neurosphere.o, dimred = "umap",  x_lab = "UMAP-1", y_lab = "UMAP-2",
																			color_by = "log2umi", 
																			color.name =  bquote(paste('log'['2'],"(TotalUMIs)")))

mp <- plot_grid(umap.sam.scat.p, umap.celltype.scat.p, umap.umi.scat.p, 
								umap.cyclone.scat.p, umap.seurat.scat.p, umap.cc.scat.p, 
								nrow = 2, ncol = 3, label_size = 10, labels = "auto", align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.NeuroUMAP.pdf"), mp,
					base_height = 2, base_width = 2*1.5, nrow = 2, ncol = 3, device = cairo_pdf)
