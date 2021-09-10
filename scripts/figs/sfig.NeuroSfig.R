rm(list=ls())
source(here::here("scripts/utils.R"))

library(ggrepel)
library(ggforce)
neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))

neurosphere.scat.p <- plotCirclePlotProjection2(neurosphere.o, r = 7, label.x = 6, label.y = 1.5, label.size = 3)

### compare two theta
tmp.df <- data.frame(theta = neurosphere.o$tricyclePosition, color = neurosphere.o$CCStage,  pca.theta = as.numeric(coord2rad(reducedDim(neurosphere.o, "go.pca")[, 1:2])))
tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)

neurosphere.theta.p <- ggplot(tmp.df, aes(x = pca.theta, y = theta, color = theta)) +
	geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize  = metadata(neurosphere.o)$point.size, alpha = metadata(neurosphere.o)$point.alpha) +
	geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize  = metadata(neurosphere.o)$point.size, alpha = metadata(neurosphere.o)$point.alpha) +
	scale_color_gradientn(name = NULL, limits = range(0, 2 * pi), 
												breaks = seq(from = 0, to = 2 * pi, length.out = 500) ,
												colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"), 
												guide = FALSE) +
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
											color = neurosphere.o$tricyclePosition)
tmp0.df$color <- fct_explicit_na(tmp0.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
grid.s.ldf <- computeVelocityOnGrid( reducedDim(neurosphere.o, "tricycleEmbedding")[, 1:2], proj_v.m[, 1:2],  
																		 x1range = range(reducedDim(neurosphere.o, "tricycleEmbedding")[, 1]), x2range = range(reducedDim(neurosphere.o, "tricycleEmbedding")[, 2]),
																		 adjust = FALSE, density = 0.8, min_mass = 2)
tmp.df <- grid.s.ldf$grid
xlim.v <- range(tmp0.df$x)
ylim.v <- range(tmp0.df$y)
hue.colors <- c("#2E22EA", "#9E3DFB", "#F86BE2",
								"#FCCE7B", "#C4E416", "#4BBA0F",
								"#447D87", "#2C24E9")

velocity.p <- ggplot(tmp.df, aes(x = x, y = y)) + 
	geom_scattermore(data = tmp0.df, aes(color= color), pointsize  = metadata(neurosphere.o)$point.size, alpha = 0.4) +
	geom_segment(aes(xend = x + dx * 0.3, yend = y + dy * 0.3), arrow = arrow(length = unit(0.03,"cm")), size = 0.2, color = "gray10") +
	#scale_color_manual(values = c(ccColors.v, "grey"), name = "CC Stage", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA")) +
	scale_color_gradientn(name = "tricyclePosition", limits = range(0, 2 * pi), breaks = seq(from = 0, to = 2 * pi, length.out = 500), colors = hue.colors, guide = FALSE) +
	#xlim(xlim.v + c(-2, 2)) + ylim(ylim.v + c(-2, 3)) +
	labs( y = py_lab, 
				x =  px_lab,
				title = str_c( metadata(neurosphere.o)$dataname, " (n=", ncol(neurosphere.o), ")")) +
	geom_arc(aes(x0 = x, y0 = y, r = r,
							 start = .5 * pi, end = -1.3 * pi),
					 data = data.frame(x= -3, y = 0, r = 10),
					 size = 0.4,
					 alpha = 0.5,
					 inherit.aes = FALSE,
					 arrow = arrow(length = unit(0.06, "npc")))


### make legend 

hues.df <- data.frame(theta = seq(from = 0, to = 2 * pi, length.out = 500), colors = colorRampPalette(hue.colors)(500))
hue_text.df <- data.frame(theta = c(0, 0.5 * pi, pi, 1.5 * pi),
													label = c("0/2\u03C0", "0.5\u03C0", "\u03C0", "1.5\u03C0"),
													y.text = rep(4.2, 4),
													hjust = c(0.1, 0.5, 1, 0.5),
													vjust = c(0.5, 0.5, 0.5, 0.5))


legend.p <- ggplot(hues.df) +
	geom_rect(aes(xmin = theta - 0.001, xmax = theta + 0.001, color = colors, fill = colors), ymin = 1.5, ymax = 3, alpha = 0.6) +
	coord_polar(theta = "x", start = -pi / 2, direction = -1, clip = "on") +
	scale_color_identity() +
	scale_fill_identity() +
	guides(fill = FALSE, color = FALSE) +
	theme_void() +
	ylim(c(0, 5)) + 
	geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
							 data = data.frame(x1 = 0, x2 = 1.75 * pi, y1 = 3.5, y2 = 3.5),
							 inherit.aes = FALSE,
							 alpha = 0.5,
							 arrow = arrow(length = unit(0.05, "npc"))) +
	geom_text(data = hue_text.df, aes_string(x = "theta", y = "y.text", label = "label", hjust = "hjust"), size = 3) +
	theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))


### neuro_projection

mp <- plot_grid(neurosphere.scat.p + theme(legend.position = c(0, 1),
																							legend.justification = c(0, 1),
																							legend.key = element_blank(),
																							legend.key.size = unit(7, "pt")), 
								neurosphere.theta.p + theme(legend.position = "none"),
								legend.p, 
								# velocity.p + theme(legend.position = "none"), 
									nrow = 1, ncol = 3, rel_widths = c(1, 1, 0.6), labels = c("a", "b"), align = "h", axis = "lr")
save_plot(here::here("figs", "sfigs", "sfig.NeuroSfig.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 2.5, device = cairo_pdf)


mp <- plot_grid(velocity.p + theme(legend.position = "none"),
								legend.p,
								# velocity.p + theme(legend.position = "none"), 
								nrow = 1, ncol = 2, rel_widths = c(1, 1), labels = NULL, align = "h", axis = "lr")
save_plot(here::here("figs", "sfigs", "sfig.NeuroVelocity.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 2, device = cairo_pdf)

