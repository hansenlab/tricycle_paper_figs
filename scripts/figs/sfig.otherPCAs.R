rm(list=ls())
source(here::here("scripts/utils.R"))


endo.o <- qread(here::here("data/plotdata/endo.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))
mHSC.o <- qread(here::here("data/plotdata/mHSC.qs"))
HeLa1.o <- qread(here::here("data/plotdata/HeLa1.qs"))
HeLa2.o <- qread(here::here("data/plotdata/HeLa2.qs"))
mRetina.o <- qread(here::here("data/plotdata/mRetina.qs"))


pca.lp <- lapply(list(hipp.o, endo.o, mHSC.o, HeLa1.o, HeLa2.o, mRetina.o), function(sce) {
	plotScatCC(sce.o = sce, dimred = "go.pca", 
						 x_lab = str_c("PC1", sprintf(" (%i%%)", round(attr(reducedDim(sce, "go.pca"), "percentVar")[1]))),
						 y_lab = str_c("PC2", sprintf(" (%i%%)", round(attr(reducedDim(sce, "go.pca"), "percentVar")[2]))),
						 title = str_c( metadata(sce)$dataname, " PCA (n=", ncol(sce), ")"))
}) %>% lapply( function(x) x + theme(legend.position = "none") )
pca.lp[[3]] <- pca.lp[[3]] + theme(legend.position = c(1, 0),
																	 legend.justification = c(1, 0),
																	 legend.key = element_blank(),
																	 legend.key.size = unit(6.5, "pt"))

mp <- plot_grid(plotlist = pca.lp,
								nrow = 2, ncol = 3,  label_size = 10, labels = c("auto"), align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.otherPCAs.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 3, device = cairo_pdf)


