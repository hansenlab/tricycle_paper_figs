rm(list=ls())
source(here::here("scripts/utils.R"))

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))

neurosphere_full.o <- qread(here::here("data/neurosphere.qs"))
hipp_full.o <- qread(here::here("data/hipp.qs"))


### projection using ksday0
common_genes.v <- intersect(rownames(metadata(neurosphere.o)$rotation), rownames(hipp_full.o))

reducedDim(neurosphere.o, "common") <- project_cycle_space(assay(neurosphere_full.o, "log.s"), ref.m = metadata(neurosphere.o)$rotation[common_genes.v, ])
reducedDim(hipp.o, "common") <- project_cycle_space(assay(hipp_full.o, "log.s"), ref.m = metadata(neurosphere.o)$rotation[common_genes.v, ])


neurosphere.p <- plotEmbScatCyclic(neurosphere.o, dimred = "common", x_lab = px_lab, y_lab = py_lab, title = str_c(metadata(neurosphere.o)$dataname, " projected by mNeurosphere\n(", length(common_genes.v)," common genes)"))
hipp.p <- plotEmbScatCyclic(hipp.o, dimred = "common", x_lab = px_lab, y_lab = py_lab, title = str_c(metadata(hipp.o)$dataname, " projected by mNeurosphere\n(", length(common_genes.v)," common genes)"))

### neuro_projection

mp <- plot_grid(neurosphere.p + theme(legend.position = c(1, 1),
																					 legend.justification = c(1, 1),
																					 legend.key = element_blank(),
																					 legend.key.size = unit(6.5, "pt")), 
								hipp.p + theme(legend.position = "none"),
								nrow = 1, ncol = 2, rel_widths = c(1, 1, 1), labels = c("auto"), align = "h", axis = "lr")
save_plot(here::here("figs", "sfigs", "sfig.commonGeneProj.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 2, device = cairo_pdf)
