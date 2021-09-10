rm(list=ls())
source(here::here("scripts/utils.R"))


endo.o <- qread(here::here("data/plotdata/endo.qs"))

ori_tricycleEmbedding.v <- endo.o$tricyclePosition


endo.top2a.p <- plotLoess2(endo.o, "top2a", col.outname = "Top2A", title = str_c( metadata(endo.o)$dataname, " Top2A"), y_lab = bquote(paste('log'['2'],'(expression)')))


endo.o <- estimate_cycle_position(endo.o,  dimred = "go.pca")

endo.top2apca.p <- plotLoess2(endo.o, "top2a", col.outname = "Top2A", title = str_c( metadata(endo.o)$dataname, " Top2A"), y_lab = bquote(paste('log'['2'],'(expression)')),
															x_lab = bquote(paste("PCA \u03B8")))

mp <- plot_grid(endo.top2a.p, endo.top2apca.p, 
								nrow = 1, ncol = 2, label_size = 10, labels = "auto", align = "v", axis = "lr")


save_plot(here::here("figs", "response_figs", "response.endopcatheta.pdf"), mp,
					base_height = 2, base_width = 2 * 1.5, nrow = 1, ncol = 2, device = cairo_pdf)



