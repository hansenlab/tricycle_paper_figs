rm(list=ls())
source(here::here("scripts/utils.R"))

mRetina.o <- qread(here::here("data/plotdata/mRetina.qs"))


### facet by both day and ct
tmp.df <- data.frame(x = reducedDim(mRetina.o, "tricycleEmbedding")[, 1], y = reducedDim(mRetina.o, "tricycleEmbedding")[, 2],
										 cc = mRetina.o$CCStage, age = mRetina.o$age, ct = mRetina.o$cell_type)
tmp.df$cc <- fct_explicit_na(tmp.df$cc, na_level = "NA") %>% fct_relevel("NA", after = Inf)
scale_color <- scale_color_manual(values = c(ccColors.v, "grey"), name = "Cellcycle", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA"))

color10.v <- brewer.pal(12, "Set3")[-c(8, 9)]
### day ct
day.ct.lp <- lapply(levels(tmp.df$age), function(a) {
	tmp1.df <- tmp.df %>% dplyr::filter(`age` == a)
	day.p <- ggplot(tmp1.df, aes(x = x, y = y, color = ct)) +
		geom_scattermore(pointsize= 3.5,  alpha = 0.8) +
		scale_color_manual(values = color10.v, name = "Cell type", labels = levels(tmp.df$ct), limits = levels(tmp.df$ct)) +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( y = py_lab, x = px_lab, title =  str_c(a, " (n=", sum(tmp.df$age == a, na.rm = TRUE), ")")) +
		xlim(range(tmp.df$x)) + ylim(range(tmp.df$y))
	
	lp <- lapply(seq_len(nlevels(factor(tmp.df$ct))), function(idx) {
		p <- ggplot(tmp1.df, aes(x = x, y = y, color = ct)) +
			geom_scattermore(data = tmp1.df %>% dplyr::filter(`ct` != levels(factor(tmp.df$ct))[idx]), pointsize = 2.1,  alpha = 0.8, color = "gray90", show.legend = FALSE) +
			geom_scattermore(data = tmp1.df %>% dplyr::filter(`ct` == levels(factor(tmp.df$ct))[idx]), pointsize = 3.5,  alpha = 0.8) +
			scale_color_manual(values = color10.v, name = "Cell type", labels = levels(tmp.df$ct), limits = levels(tmp.df$ct)) +
			guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
			labs( y = py_lab, x = px_lab, title = str_c(levels(factor(tmp.df$ct))[idx], " (n=", sum(tmp1.df$ct == levels(factor(tmp.df$ct))[idx], na.rm = TRUE), ")")) +
			xlim(range(tmp.df$x)) + ylim(range(tmp.df$y)) + theme(legend.position = "none")
		return(p)
	})
	
	legend.p <- get_legend(day.p)
	
	p <- plot_grid(plotlist = c(list(day.p + theme(legend.position = "none"), legend.p), lp),
								 nrow = 1, ncol = 12, rel_widths = c(1, 0.8, rep(1, 10)),label_size = 10, labels = NULL)
	return(p)
})


mp <- plot_grid(plotlist = day.ct.lp,
							 nrow = length(day.ct.lp), ncol = 1, label_size = 10, labels = NULL)

save_plot(here::here("figs", "sfigs", "sfig.retinadayct.pdf"),  mp,
					base_height = 2, base_width = 2*1.2, nrow = length(day.ct.lp), ncol = 11.8, device = cairo_pdf)
save_plot(here::here("figs", "sfigs", "sfig.retinadayct.jpg"),  mp,
					base_height = 2, base_width = 2*1.2, nrow = length(day.ct.lp), ncol = 11.8, type = "cairo", dpi = 130)


