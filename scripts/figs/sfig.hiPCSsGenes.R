rm(list=ls())
source(here::here("scripts/utils.R"))



hiPSCs.o <- qread(here::here("data/plotdata/hiPSCs.qs"))
hiPSCs_full.o <- qread(here::here("data/hiPSCs.qs"))


na.omit(match(c("CDK1", "UBE2C",   "SMC2", "SMC4", "UBE2S", "CKAP2", "MCM6", "HELLS"), rowData(hiPSCs_full.o)$Gene))


SelGenes.lp <- do.call(c, lapply(c("CDK1", "UBE2C",   "SMC2", "SMC4", "UBE2S", "CKAP2", "MCM6", "HELLS"), function(g) {
	colData(hiPSCs.o)[, g] <- assay(hiPSCs_full.o, "log.s")[match(g, rowData(hiPSCs_full.o)$Gene), ]
	p1 <- plotLoess3(sce.o = hiPSCs.o, col.name = g, x_val = "tricyclePosition", col.outname = g, addR2 = TRUE, r2size = 2.5) + theme(legend.position = "none")
	p2 <- plotLoess3(sce.o = hiPSCs.o, col.name = g,  x_val = "Hsiao.theta", x_lab = "FUCCI pseudotime", col.outname = g, addR2 = TRUE, r2size = 2.5)  + theme(legend.position = "none")
	return(list(p1, p2))
}))

SelGenes.mp <- plot_grid(plotlist = c(list(SelGenes.lp[[1]] + theme(legend.position = c(1, 0),
																																		legend.justification = c(1,0))), 
																			SelGenes.lp[2:length(SelGenes.lp)]),
												 nrow = 4, ncol = 4, label_size = 10, labels = as.vector(rbind(letters[1:8], rep("", 8))), align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.hiPCSsGenes.pdf"), SelGenes.mp,
					base_height = 2, base_width = 2*1.2, nrow = 4, ncol = 4, device = cairo_pdf)





