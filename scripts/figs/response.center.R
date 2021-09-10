rm(list=ls())
source(here::here("scripts/utils.R"))


endo.o <- qread(here::here("data/plotdata/endo.qs"))

ori_tricycleEmbedding.v <- endo.o$tricyclePosition

endo.o <- estimate_cycle_position(endo.o, center.pc1 = -3, center.pc2 = 0)


tmp.df <- data.frame(theta = ori_tricycleEmbedding.v,
										 new = endo.o$tricyclePosition,
										 color = endo.o$CCStage)
tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)

theta.p <- ggplot(tmp.df, aes(x = theta, y = new, color = color)) +
	geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize  = metadata(endo.o)$point.size, alpha = metadata(endo.o)$point.alpha) +
	geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize  = metadata(endo.o)$point.size, alpha = metadata(endo.o)$point.alpha) +
	scale_color_manual(values = c(ccColors.v, "grey"), name = "SchwabeCC", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA")) +
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
	labs( 	x = bquote(paste('CC'['ns']," Position \u03B8 using center (0,0)")), 
				 y = bquote(paste('CC'['ns']," Position \u03B8 using center (-3,0)")),
				 title = str_c(metadata(endo.o)$dataname, " (N=", ncol(endo.o), ")")) +
	annotate(geom = "text", x = 1.5, y = 5, size = 3, hjust = 0.5, vjust = 0.5,
					 label = str_c("Circular correlation \n \u03C1=", sprintf("%.3f", Directional::circ.cor1(tmp.df$theta, tmp.df$new , rads = T)[1])), parse = FALSE) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi , labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2 )* pi) +
	scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi , labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2)* pi) +
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				legend.key = element_blank(),
				legend.key.size = unit(6.5, "pt"))


save_plot(here::here("figs", "response_figs", "response.center.pdf"), theta.p,
					base_height = 2, base_width = 2 * 1.5, nrow = 1, ncol = 1, device = cairo_pdf)

