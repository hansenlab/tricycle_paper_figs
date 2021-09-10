rm(list=ls())
source(here::here("scripts/utils.R"))

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))



umap.cyclone.scat.p <- plotEmbScat(sce.o = neurosphere.o, dimred = "umap",  x_lab = "UMAP-1", y_lab = "UMAP-2",
																	 color_by = "cyclone", colors.v = cc3Colors.v,
																	 labels.v = cc3Labels.v, color.name = 'Cyclone')
umap.seurat.scat.p <- plotEmbScat(sce.o = neurosphere.o, dimred = "umap",  x_lab = "UMAP-1", y_lab = "UMAP-2",
																	color_by = "SeuratCC", colors.v = cc3Colors.v,
																	labels.v = cc3Labels.v, color.name = 'SeuratCC')


umap.tricycle.scat.p <- plotEmbScatCyclic(sce.o = neurosphere.o, dimred = "umap", 
									x_lab = "UMAP-1", y_lab = "UMAP-2")


mp <- plot_grid(
								umap.seurat.scat.p, # umap.cyclone.scat.p, 
								umap.tricycle.scat.p, circle_scale_legend(text.size = 1.8, ymax = 4.5, y.outer = 2.5, y.text = 3.1),
								nrow = 1, ncol = 3, label_size = 10, labels = c("a", "b", ""), align = "hv", axis = "tblr")

save_plot(here::here("figs", "response_figs", "response.NeuroSfig.pdf"), mp,
					base_height = 2, base_width = 2*1.5, nrow = 1, ncol = 3, device = cairo_pdf)





