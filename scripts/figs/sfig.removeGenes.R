rm(list=ls())
source(here::here("scripts/utils.R"))


neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
neurosphere_full.o <- qread(here::here("data/neurosphere.qs"))

### neurosphere
theta.v <- neurosphere.o$tricyclePosition
ref.m <- metadata(neurosphere.o)$rotation

neurosphere.scat.p <- plotScatCC(sce.o = neurosphere.o, dimred = "tricycleEmbedding", x_lab = px_lab, y_lab = py_lab, title = str_c(metadata(neurosphere.o)$dataname, " (N.genes=", 500, ")")) +
	theme(legend.position = "none")


neurosphere.lp <- lapply(c(400, 300, 200, 100, 50), function(n) {
	set.seed(n )
	rotation.m <- ref.m[sample(seq_len(nrow(ref.m)), size = n, replace = FALSE), ]
	reducedDim(neurosphere.o, "new") <- project_cycle_space(assay(neurosphere_full.o, "log.s"), ref.m = rotation.m)
	scat.p <- plotScatCC(sce.o = neurosphere.o, dimred = "new", x_lab = px_lab, y_lab = py_lab, title = str_c(metadata(neurosphere.o)$dataname, " (N.genes=", n, ")")) + theme(legend.position = "none")

	tmp.df <- data.frame(theta = theta.v,
											 new = as.numeric(circular::coord2rad(reducedDim(neurosphere.o, "new")[, 1:2])),
											 color = neurosphere.o$CCStage)
	tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
	
	theta.p <- ggplot(tmp.df, aes(x = theta, y = new, color = color)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize  = metadata(neurosphere.o)$point.size, alpha = metadata(neurosphere.o)$point.alpha) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize  = metadata(neurosphere.o)$point.size, alpha = metadata(neurosphere.o)$point.alpha) +
		scale_color_manual(values = c(ccColors.v, "grey"), name = "CC Stage", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA")) +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
		labs( 	x = bquote(paste('CC'['ns']," Position \u03B8 (N.genes=500)")), 
					y = bquote(paste('CC'['ns']," Position \u03B8 (N.genes=", .(n), ")")),
					title = str_c(metadata(neurosphere.o)$dataname, " (N.genes=", n, ")")) +
		annotate(geom = "text", x = 2, y = 5, size = 3, hjust = 0.5, vjust = 0.5,
						 label = str_c("Circular correlation \n \u03C1=", sprintf("%.3f", Directional::circ.cor1(tmp.df$theta, tmp.df$new , rads = T)[1])), parse = FALSE) +
		scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi , labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2 )* pi) +
		scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi , labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2)* pi) +
		theme(legend.position = "none")
	return(list(scat.p = scat.p, theta.p = theta.p))
	
})
neurosphere.lp[[1]]$scat.p <- neurosphere.lp[[1]]$scat.p + theme(legend.position = c(0, 1),
																												 legend.justification = c(0, 1))


mp <- plot_grid(plotlist = c(lapply(neurosphere.lp, "[[", 1), lapply(neurosphere.lp, "[[", 2)),
								nrow = 2, ncol = 5, label_size = 10, labels = c("a", rep("", 4), "b"))
save_plot(here::here("figs", "sfigs", "sfig.removeGenes.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 5, device = cairo_pdf)


