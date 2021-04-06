rm(list=ls())
source(here::here("scripts/utils.R"))

mRetina.o <- qread(here::here("data/plotdata/mRetina.qs"))


projection.day.lp <- plotEmbScat(mRetina.o, dimred = "tricycleEmbedding", color_by = "CCStage", color.name = "CC Stage", labels.v = ccLabels.v,colors.v  = ccColors.v, facet_by = "age",  x_lab = px_lab, y_lab = py_lab)


projection.day.lp <- lapply(projection.day.lp, function(x) x + theme(legend.position = "none"))
projection.day.lp[[1]] <- projection.day.lp[[1]] + theme(legend.position = c(0, 1),
																											 legend.justification = c(0, 1))

mp <- plot_grid(plotlist = projection.day.lp,
							 nrow = 3, ncol = 4, label_size = 10, labels = "auto")
save_plot(here::here("figs", "sfigs", "sfig.retinaday.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 3, ncol = 4, device = cairo_pdf)

