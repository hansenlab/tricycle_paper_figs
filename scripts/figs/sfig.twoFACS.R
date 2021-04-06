rm(list=ls())
source(here::here("scripts/utils.R"))

mESC.o <- qread(here::here("data/plotdata/mESC.qs"))
hESC.o <- qread(here::here("data/plotdata/hESC.qs"))



mESC.proj.p <- plotCirclePlotProjection(mESC.o, r = 15, label.x = 17, label.y = 5, 
												 col.name = "stage", color.name = "FACS", 
												 colors.v = cc3Colors.v, labels.v = cc3Labels.v) +
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1),
				legend.key = element_blank(),
				legend.key.size = unit(7, "pt"))

hESC.proj.p <- plotCirclePlotProjection(hESC.o, r = 10, label.x = 8, label.y = 3, 
																				col.name = "stage", color.name = "FACS", 
																				colors.v = cc3Colors.v, labels.v = cc3Labels.v) +
	theme(legend.position = "none")

mESC.top2a.p <- plotLoess(mESC.o, col.outname = "Top2A", col.name = "top2a",
													color.name = "FACS", color.var = "stage",
													colors.v = cc3Colors.v, labels.v = cc3Labels.v) +
	theme(legend.position = "none")

hESC.top2a.p <- plotLoess(hESC.o, col.name = "top2a", color.name = "FACS",
													color.var = "stage", colors.v = cc3Colors.v,
													labels.v = cc3Labels.v) +
	theme(legend.position = "none")


mp <- plot_grid(mESC.proj.p,
								mESC.top2a.p,
								hESC.proj.p,
								hESC.top2a.p,
								nrow = 2, ncol = 2,  label_size = 10, labels = c("auto"), align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.twoFACS.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 2, device = cairo_pdf)





