rm(list=ls())
source(here::here("scripts/utils.R"))


endo.o <- qread(here::here("data/plotdata/endo.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))
HeLa2.o <- qread(here::here("data/plotdata/HeLa2.qs"))
mRetina.o <- qread(here::here("data/plotdata/mRetina.qs"))

metadata(endo.o)$point.size <- 4.01

hipp.scat.p <- plotCirclePlotProjection2(hipp.o, r = 3, label.x = 2, label.y = 0.5)
endo.scat.p <- plotCirclePlotProjection2(endo.o, r = 2, label.x = 1.8, label.y = 0.5)
mRetina.scat.p <- plotCirclePlotProjection2(mRetina.o, r = 3, label.x = 2.5, label.y = 0.8)
HeLa2.scat.p <- plotCirclePlotProjection2(HeLa2.o, r = 4, label.x = 3.7, label.y = 1)


hipp.top2a.p <- plotLoess2(hipp.o, "top2a", col.outname = "Top2A", title = str_c( metadata(hipp.o)$dataname, " Top2A"), y_lab = bquote(paste('log'['2'],'(expression)')))
endo.top2a.p <- plotLoess2(endo.o, "top2a", col.outname = "Top2A", title = str_c( metadata(endo.o)$dataname, " Top2A"), y_lab = bquote(paste('log'['2'],'(expression)')))
mRetina.top2a.p <- plotLoess2(mRetina.o, "top2a", col.outname = "Top2A", title = str_c( metadata(mRetina.o)$dataname, " Top2A"), y_lab = bquote(paste('log'['2'],'(expression)')))
HeLa2.top2a.p <- plotLoess2(HeLa2.o, "top2a", title = str_c( metadata(HeLa2.o)$dataname, " TOP2A"), y_lab = bquote(paste('log'['2'],'(expression)')))


hipp.cyclic.p <- plotEmbScatCyclic(sce.o = hipp.o, dimred = "umap", 
																	x_lab = "UMAP-1", y_lab = "UMAP-2", title = str_c( metadata(hipp.o)$dataname))
endo.cyclic.p <- plotEmbScatCyclic(sce.o = endo.o, dimred = "umap",
																	x_lab = "UMAP-1", y_lab = "UMAP-2", title = str_c( metadata(endo.o)$dataname ))
mRetina.cyclic.p <- plotEmbScatCyclic(sce.o = mRetina.o, dimred = "umap",
																		 x_lab = "UMAP-1", y_lab = "UMAP-2", title = str_c( metadata(mRetina.o)$dataname ))
HeLa2.cyclic.p <- plotEmbScatCyclic(sce.o = HeLa2.o, dimred = "umap", 
																	 x_lab = "UMAP-1", y_lab = "UMAP-2", title = str_c( metadata(HeLa2.o)$dataname ))




mp <- plot_grid(hipp.scat.p +  theme(legend.position = "none") ,
								# theme(legend.position = c(1, 1),
								# 												legend.justification = c(1, 1),
								# 											 legend.key = element_blank(),
								# 											 legend.key.size = unit(7, "pt")),
								endo.scat.p + theme(legend.position = "none", plot.title = element_text(face = "plain", size = 7.5, hjust = 0.5)), 
								mRetina.scat.p + theme(legend.position = "none"), 
								HeLa2.scat.p + theme(legend.position = "none"), 
								ggplot() + theme_nothing(),
								hipp.top2a.p + theme(legend.position = "none"), 
								endo.top2a.p + theme(legend.position = "none"), 
								mRetina.top2a.p + theme(legend.position = "none"), 
								HeLa2.top2a.p + theme(legend.position = "none"), 
								circle_scale_legend(text.size = 1.8, y.inner =  0.9, ymax = 4.5, y.outer = 2.1, y.text = 3.3),
								hipp.cyclic.p + theme(legend.position = "none"), 
								endo.cyclic.p + theme(legend.position = "none"), 
								mRetina.cyclic.p + theme(legend.position = "none"), 
								HeLa2.cyclic.p + theme(legend.position = "none"), 
								ggplot() + theme_nothing(),
								nrow = 3, ncol = 5, rel_widths = c(1, 1, 1, 1, 0.45), labels = c("a", rep("", 4), "b", rep("", 4), "c", rep("", 4)))


save_plot(here::here("figs", "main", "main.projections.pdf"), mp,
					base_height = 2 / 1.40, base_width = 2*1.2 / (1.2*1.5), nrow = 3, ncol = 4.5, device = cairo_pdf)


