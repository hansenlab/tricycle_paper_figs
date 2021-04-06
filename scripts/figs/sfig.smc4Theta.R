rm(list=ls())
source(here::here("scripts/utils.R"))


endo.o <- qread(here::here("data/plotdata/endo.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))
mHSC.o <- qread(here::here("data/plotdata/mHSC.qs"))
HeLa1.o <- qread(here::here("data/plotdata/HeLa1.qs"))
HeLa2.o <- qread(here::here("data/plotdata/HeLa2.qs"))
mRetina.o <- qread(here::here("data/plotdata/mRetina.qs"))


hipp.smc4.p <- plotLoess(hipp.o, "smc4")
endo.smc4.p <- plotLoess(endo.o, "smc4")
mHSC.smc4.p <- plotLoess(mHSC.o, "smc4")
HeLa1.smc4.p <- plotLoess(HeLa1.o, "smc4")
HeLa2.smc4.p <- plotLoess(HeLa2.o, "smc4")
mRetina.smc4.p <- plotLoess(mRetina.o, "smc4")


mp <- plot_grid(hipp.smc4.p + theme(legend.position = c(1, 1),
																		legend.justification = c(1, 1),
																		legend.key = element_blank(),
																		legend.key.size = unit(6.5, "pt")),
								endo.smc4.p + theme(legend.position = "none"),
								mRetina.smc4.p + theme(legend.position = "none"),
								HeLa2.smc4.p + theme(legend.position = "none"),
								mHSC.smc4.p + theme(legend.position = "none"),
								HeLa1.smc4.p + theme(legend.position = "none"),
								nrow = 2, ncol = 3,  label_size = 10, labels = c("auto"), align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.smc4Theta.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 3, device = cairo_pdf)

