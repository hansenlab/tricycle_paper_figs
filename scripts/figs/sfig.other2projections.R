rm(list=ls())
source(here::here("scripts/utils.R"))


mHSC.o <- qread(here::here("data/plotdata/mHSC.qs"))
HeLa1.o <- qread(here::here("data/plotdata/HeLa1.qs"))
hca_pancreas.o <- qread(here::here("data/fetal/Pancreas.qs"))

hca_pancreas.o <- estimate_cycle_position(hca_pancreas.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
metadata(hca_pancreas.o)$point.size <- ifelse(ncol(hca_pancreas.o) > 1000000, 1.01, 2.01)
metadata(hca_pancreas.o)$point.alpha <- 0.7
metadata(hca_pancreas.o)$species <- "human"
metadata(hca_pancreas.o)$dataname <- str_c("HCA ", "Pancreas")
hca_pancreas.o$top2a <- assay(hca_pancreas.o, "log.s")["TOP2A", ]
reducedDim(hca_pancreas.o, "umap") <- data.frame(umap1 = hca_pancreas.o$Main_cluster_umap_1, umap2 = hca_pancreas.o$Main_cluster_umap_2)




mHSC.scat.p <- plotCirclePlotProjection2(mHSC.o, r = 20, label.x = 15, label.y = 5)
HeLa1.scat.p <- plotCirclePlotProjection2(HeLa1.o, r = 5, label.x = 5, label.y = 1)
hca_pancreas.scat.p <- plotCirclePlotProjection2(hca_pancreas.o, r = 2, label.x = 1, label.y = 0.3)



mHSC.top2a.p <- plotLoess2(mHSC.o, "top2a", col.outname = "Top2A")
HeLa1.top2a.p <- plotLoess2(HeLa1.o, "top2a")
hca_pancreas.top2a.p <- plotLoess2(hca_pancreas.o, "top2a")

mHSC.cyclic.p <- plotEmbScatCyclic(sce.o = mHSC.o, dimred = "umap", 
																	 x_lab = "UMAP-1", y_lab = "UMAP-2", title = str_c( metadata(mHSC.o)$dataname, " UMAP (n=", ncol(mHSC.o), ")"))

HeLa1.cyclic.p <- plotEmbScatCyclic(sce.o = HeLa1.o, dimred = "umap", 
																		x_lab = "UMAP-1", y_lab = "UMAP-2", title = str_c( metadata(HeLa1.o)$dataname, " UMAP (n=", ncol(HeLa1.o), ")"))

hca_pancreas.cyclic.p <- plotEmbScatCyclic(sce.o = hca_pancreas.o, dimred = "umap", 
																		x_lab = "UMAP-1", y_lab = "UMAP-2", title = str_c( metadata(hca_pancreas.o)$dataname, " UMAP (n=", ncol(hca_pancreas.o), ")"))


mp <- plot_grid(mHSC.scat.p +  theme(legend.position = "none") ,
								HeLa1.scat.p + theme(legend.position = "none"), 
								hca_pancreas.scat.p + theme(legend.position = "none"), 
								ggplot() + theme_nothing(),
								mHSC.top2a.p + theme(legend.position = "none"), 
								HeLa1.top2a.p + theme(legend.position = "none"), 
								hca_pancreas.top2a.p + theme(legend.position = "none"), 
								circle_scale_legend(),
								mHSC.cyclic.p + theme(legend.position = "none"), 
								HeLa1.cyclic.p + theme(legend.position = "none"), 
								hca_pancreas.cyclic.p + theme(legend.position = "none"), 
								ggplot() + theme_nothing(),
								nrow = 3, ncol = 4, rel_widths = c(1, 1, 1, 0.5), labels = c("a", rep("", 3), "b", rep("", 3), "c", rep("", 3)), align = "h", axis = "lr")


save_plot(here::here("figs", "sfigs", "sfig.other2projections.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 3, ncol = 3.5, device = cairo_pdf)






