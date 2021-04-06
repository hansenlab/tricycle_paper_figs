rm(list=ls())
source(here::here("scripts/utils.R"))

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))
endo.o <- qread(here::here("data/plotdata/endo.qs"))
mRetina.o <- qread(here::here("data/plotdata/mRetina.qs"))
mHSC.o <- qread(here::here("data/plotdata/mHSC.qs"))
HeLa1.o <- qread(here::here("data/plotdata/HeLa1.qs"))
HeLa2.o <- qread(here::here("data/plotdata/HeLa2.qs"))
hfIntestineSub.o <- qread(here::here("data/plotdata/hfIntestineSub.qs"))

hU2OS.o <- qread(here::here("data/plotdata/hU2OS.qs"))
hiPSCs.o <- qread(here::here("data/plotdata/hiPSCs.qs"))
mESC.o <- qread(here::here("data/plotdata/mESC.qs"))
hESC.o <- qread(here::here("data/plotdata/hESC.qs"))




seuratcc.pl <- do.call(c, lapply(list(neurosphere.o, hipp.o, endo.o, mRetina.o, mHSC.o, HeLa1.o, HeLa2.o, hfIntestineSub.o, hU2OS.o, hiPSCs.o, mESC.o, hESC.o), function(sce.o) {
	if (metadata(sce.o)$dataname == "mRetina") {
		p.l <- list(plotScatCC(sce.o = sce.o, dimred = "tricycleEmbedding", 
													 x_lab = px_lab, y_lab = py_lab, col.name = "SeuratCC", 
													 color.name = "SeuratCC", colors.v = cc3Colors.v, labels.v = cc3Labels.v) +
									theme(legend.position = "none"),
								#plotThetaDen(sce.o, color.var = "SeuratCC", color.name = "SeuratCC", colors.v = cc3Colors.v, labels.v = cc3Labels.v, addP = TRUE)$circular + theme(legend.position = "none"),
								ggplot() + theme_nothing(),
								plotGeneRidge(sce.o, col.name = "top2a", color.var = "SeuratCC", color.name = "SeuratCC", colors.v = cc3Colors.v, labels.v = cc3Labels.v) + theme(legend.position = "none"))
	} else {
		p.l <- list(plotScatCC(sce.o = sce.o, dimred = "tricycleEmbedding", 
													 x_lab = px_lab, y_lab = py_lab, col.name = "SeuratCC", 
													 color.name = "SeuratCC", colors.v = cc3Colors.v, labels.v = cc3Labels.v) +
									theme(legend.position = "none"),
								#plotThetaDen(sce.o, color.var = "SeuratCC", color.name = "SeuratCC", colors.v = cc3Colors.v, labels.v = cc3Labels.v, addP = TRUE)$circular + theme(legend.position = "none"),
								plotSilhouetteBox(sce.o, color.var = "SeuratCC", color.name = "SeuratCC", colors.v = cc3Colors.v, labels.v = cc3Labels.v) + theme(legend.position = "none"),
								plotGeneRidge(sce.o, col.name = "top2a", color.var = "SeuratCC", color.name = "SeuratCC", colors.v = cc3Colors.v, labels.v = cc3Labels.v) + theme(legend.position = "none"))
	}
	print(metadata(sce.o)$dataname)

	return(p.l)
}))




seuratcc.pl[[1]] <- seuratcc.pl[[1]] + theme(legend.position = c(1,1 ),
																					 legend.justification = c(1, 1),
																					 legend.key = element_blank(),
																					 legend.key.size = unit(6.5, "pt"))
seuratcc.pl[[2]] <- seuratcc.pl[[2]] + theme(legend.position = c(1, 0),
																					 legend.justification = c(1, 0),
																					 legend.key = element_blank(),
																					 legend.key.size = unit(6.5, "pt"))
seuratcc.pl[[3]] <- seuratcc.pl[[3]] + theme(legend.position = c(1,1 ),
																						 legend.justification = c(1, 1),
																						 legend.key = element_blank(),
																						 legend.key.size = unit(6.5, "pt"))


mp <- plot_grid(plotlist = seuratcc.pl,
								nrow = 6, ncol = 6, label_size = 10, labels = as.vector(rbind(letters[1:12], rep("", 12), rep("", 12))), align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.seuratcc.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 6, ncol = 6, device = cairo_pdf)





