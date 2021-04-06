rm(list=ls())
source(here::here("scripts/utils.R"))


mHSC.o <- qread(here::here("data/plotdata/mHSC.qs"))
HeLa1.o <- qread(here::here("data/plotdata/HeLa1.qs"))



mHSC.scat.p <- plotCirclePlotProjection(mHSC.o, r = 20, label.x = 15, label.y = 5)
HeLa1.scat.p <- plotCirclePlotProjection(HeLa1.o, r = 5, label.x = 5, label.y = 1)


mHSC.top2a.p <- plotLoess(mHSC.o, "top2a", col.outname = "Top2A")
HeLa1.top2a.p <- plotLoess(HeLa1.o, "top2a")


mHSC.cyclic.p <- plotEmbScatCyclic(sce.o = mHSC.o, dimred = "umap", 
																	 x_lab = "UMAP-1", y_lab = "UMAP-2", title = str_c( metadata(mHSC.o)$dataname, " UMAP (n=", ncol(mHSC.o), ")"))

HeLa1.cyclic.p <- plotEmbScatCyclic(sce.o = HeLa1.o, dimred = "umap", 
																		x_lab = "UMAP-1", y_lab = "UMAP-2", title = str_c( metadata(HeLa1.o)$dataname, " UMAP (n=", ncol(HeLa1.o), ")"))




mp <- plot_grid(mHSC.scat.p +  theme(legend.position = "none") ,
								HeLa1.scat.p + theme(legend.position = "none"), 
								get_legend(mHSC.scat.p),
								mHSC.top2a.p + theme(legend.position = "none"), 
								HeLa1.top2a.p + theme(legend.position = "none"), 
								get_legend(mHSC.top2a.p),
								mHSC.cyclic.p + theme(legend.position = "none"), 
								HeLa1.cyclic.p + theme(legend.position = "none"), 
								circle_scale_legend(),
								nrow = 3, ncol = 3, rel_widths = c(1, 1, 0.5), labels = c("a", rep("", 2), "b", rep("", 2), "c", rep("", 2)), align = "h", axis = "lr")


save_plot(here::here("figs", "sfigs", "sfig.other2projections.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 3, ncol = 2.5, device = cairo_pdf)






