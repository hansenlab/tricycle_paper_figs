rm(list=ls())
source(here::here("scripts/utils.R"))

library(ggrepel)

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))

neurosphere.scat.p <- plotCirclePlotProjection(neurosphere.o, r = 7, label.x = 6, label.y = 1.5, label.size = 3)

### compare two theta
tmp.df <- data.frame(theta = neurosphere.o$tricyclePosition, color = neurosphere.o$CCStage,  pca.theta = as.numeric(coord2rad(reducedDim(neurosphere.o, "go.pca")[, 1:2])))
tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)

neurosphere.theta.p <- ggplot(tmp.df, aes(x = pca.theta, y = theta, color = color)) +
	geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize  = metadata(neurosphere.o)$point.size, alpha = metadata(neurosphere.o)$point.alpha) +
	geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize  = metadata(neurosphere.o)$point.size, alpha = metadata(neurosphere.o)$point.alpha) +
	scale_color_manual(values = c(ccColors.v, "grey"), name = "CC Stage", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA")) +
	#scale_shape_manual(name = "Cellcycle", values = c(16, 17, 18)) +
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
	labs( x = "PCA \u03B8",
				y = theta_lab,
				title = "Consistency of two \u03B8s") +
	annotate(geom = "text", x = 2, y = 5, size = 3, hjust = 0.5, vjust = 0.5,
					 label = str_c("Circular correlation \n \u03C1=", sprintf("%.3f", Directional::circ.cor1(tmp.df$theta, tmp.df$pca.theta , rads = T)[1])), parse = FALSE) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi , labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2 )* pi) +
	scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi , labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2)* pi)




### add neurosphere velocity embedding
neuro_vel.o <- qread(here::here("data/neuro_velocity.o.qs"))
sum(is.na(match(colnames(neurosphere.o), colnames(neuro_vel.o))))
neuro_vel.o <- neuro_vel.o[, match(colnames(neurosphere.o), colnames(neuro_vel.o))]

proj_v.m <- projectPCA(assay(neuro_vel.o, "velocity"), metadata(neurosphere.o)$rotation)  

tmp0.df <- data.frame(x = reducedDim(neurosphere.o, "tricycleEmbedding")[, 1], y = reducedDim(neurosphere.o, "tricycleEmbedding")[, 2], 
											color = neurosphere.o$CCStage)
tmp0.df$color <- fct_explicit_na(tmp0.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
grid.s.ldf <- computeVelocityOnGrid( reducedDim(neurosphere.o, "tricycleEmbedding")[, 1:2], proj_v.m[, 1:2],  
																		 x1range = range(reducedDim(neurosphere.o, "tricycleEmbedding")[, 1]), x2range = range(reducedDim(neurosphere.o, "tricycleEmbedding")[, 2]),
																		 adjust = FALSE, density = 0.8, min_mass = 2)
tmp.df <- grid.s.ldf$grid
xlim.v <- range(tmp0.df$x)
ylim.v <- range(tmp0.df$y)

velocity.p <- ggplot(tmp.df, aes(x = x, y = y)) + 
	geom_scattermore(data = tmp0.df %>% dplyr::filter(`color` == "NA"), aes(color= color), pointsize  = metadata(neurosphere.o)$point.size, alpha = 0.4) +
	geom_scattermore(data = tmp0.df %>% dplyr::filter(`color` != "NA"), aes(color= color), pointsize  = metadata(neurosphere.o)$point.size, alpha = 0.4) +
	geom_segment(aes(xend = x + dx * 0.3, yend = y + dy * 0.3), arrow = arrow(length = unit(0.03,"cm")), size = 0.2, color = "gray10") +
	scale_color_manual(values = c(ccColors.v, "grey"), name = "CC Stage", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA")) +
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
	xlim(xlim.v) + ylim(ylim.v) +
	labs( y = py_lab, 
				x =  px_lab,
				title = str_c( metadata(neurosphere.o)$dataname, " (n=", ncol(neurosphere.o), ")")) 








### neuro_projection

mp <- plot_grid(neurosphere.scat.p + theme(legend.position = c(0, 1),
																							legend.justification = c(0, 1),
																							legend.key = element_blank(),
																							legend.key.size = unit(7, "pt")), 
								neurosphere.theta.p + theme(legend.position = "none"),
								velocity.p + theme(legend.position = "none"), 
									nrow = 1, ncol = 3, rel_widths = c(1, 1, 1), labels = c("a", "b", "c"), align = "h", axis = "lr")
save_plot(here::here("figs", "sfigs", "sfig.NeuroSfig.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 3, device = cairo_pdf)


