rm(list=ls())
source(here::here("scripts/utils.R"))


### plot top2a dynamics
top2a.lp <- lapply(c("Intestine", "Kidney", "Pancreas", "Stomach"
										), function(o) {
	sce.o <- qread(here::here("data/fetal", str_c(o, ".qs")))
	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
	
	metadata(sce.o)$point.size <- ifelse(ncol(sce.o) > 1000000, 1.01, 2.01)
	metadata(sce.o)$point.alpha <- 0.7
	metadata(sce.o)$species <- "human"
	metadata(sce.o)$dataname <- str_c("HCA ", o)
	
	if (!("TOP2A" %in% rownames(sce.o))) {
		return(NULL)
	} else {
		sce.o$top2a <- assay(sce.o, "log.s")["TOP2A", ]
		return(plotLoess2(sce.o, "top2a"))
	}
	
})

top2a_nuclei.lp <- lapply(c("Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart",
														 "Liver", "Lung", "Muscle", "Placenta", "Spleen", "Thymus"),
													 function(o) {
	sce.o <- qread(here::here("data/fetal/nuclei", str_c(o, ".qs")))
	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
	
	metadata(sce.o)$point.size <- ifelse(ncol(sce.o) > 1000000, 1.01, 2.01)
	metadata(sce.o)$point.alpha <- 0.7
	metadata(sce.o)$species <- "human"
	metadata(sce.o)$dataname <- str_c("HCA ", o)
	
	if (!("TOP2A" %in% rownames(sce.o))) {
		return(NULL)
	} else {
		sce.o$top2a <- assay(sce.o, "log.s")["TOP2A", ]
		return(plotLoess2(sce.o, "top2a"))
	}
	
})
qsave(c(top2a.lp, top2a_nuclei.lp), file = "data/HCA_top2a_p.qs")

mp <- plot_grid(
	plotlist  = c(top2a.lp, top2a_nuclei.lp), 
	nrow = 5, ncol = 3, label_size = 10, labels = c("auto"), align = "hv", axis = "tblr")

save_plot(here::here("figs", "response_figs", "response.HCA_top2a.pdf"), mp,
					base_height = 2, base_width = 2*1.3, nrow = 5, ncol = 3, device = cairo_pdf)


