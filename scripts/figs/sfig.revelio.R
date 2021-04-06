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


revelio1.pl <- lapply(list(neurosphere.o, hipp.o, endo.o, mRetina.o, mHSC.o, HeLa1.o, HeLa2.o, hfIntestineSub.o, hU2OS.o, hiPSCs.o), function(sce.o) {
	new.o <- SingleCellExperiment(reducedDims = list(dc =  metadata(sce.o)$revelio$dc), colData = data.frame(CCStage = droplevels(metadata(sce.o)$revelio$cc)))
	metadata(new.o) <- metadata(sce.o)

	p <- plotScatCC(sce.o = new.o, dimred = "dc", 
						 x_lab = "Revelio DC-1",
						 y_lab = "Revelio DC-2",
						 title = str_c( metadata(new.o)$dataname, " Revelio (n=", ncol(new.o), ")")) + theme(legend.position = "none")
	return(p)
})
revelio2.pl <- lapply(list(mESC.o, hESC.o), function(sce.o) {
	new.o <- SingleCellExperiment(reducedDims = list(dc =  metadata(sce.o)$revelio$dc), colData = data.frame(stage = sce.o$stage[as.numeric(rownames(metadata(sce.o)$revelio$dc))]))
	metadata(new.o) <- metadata(sce.o)
	
	p <- plotScatCC(sce.o = new.o, dimred = "dc", 
									x_lab = "Revelio DC-1",
									y_lab = "Revelio DC-2",
									color.name = "FACS", col.name = "stage", colors.v = cc3Colors.v, labels.v = cc3Labels.v,
									title = str_c( metadata(new.o)$dataname, " Revelio (n=", ncol(new.o), ")")) + theme(legend.position = "none")
	return(p)
})

revelio.pl <- c(revelio1.pl, revelio2.pl)


revelio.pl[[1]] <- revelio.pl[[1]] + theme(legend.position = c(1,1 ),
																						 legend.justification = c(1, 1),
																						 legend.key = element_blank(),
																						 legend.key.size = unit(6.5, "pt"))
	
mp <- plot_grid(plotlist = revelio.pl,
								nrow = 3, ncol = 4, label_size = 10, labels = "auto", align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.revelio.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 3, ncol = 4, device = cairo_pdf)



