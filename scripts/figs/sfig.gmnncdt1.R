rm(list=ls())
source(here::here("scripts/utils.R"))

hU2OS.o <- qread(here::here("data/plotdata/hU2OS.qs"))
hU2OS_full.o <- qread(here::here("data/hU2OS.qs"))

hU2OS.o$GMNN <- assay(hU2OS_full.o, "log.s")[which(rowData(hU2OS_full.o)$Gene == "GMNN"), ]
hU2OS.o$CDT1 <- assay(hU2OS_full.o, "log.s")[which(rowData(hU2OS_full.o)$Gene == "CDT1"), ]

hU2OS.o$fucci_time.theta <- hU2OS.o$fucci_time * 2 * pi

hU2OS.GMNN.p <- ggplot(data.frame(x = hU2OS.o$fucci_time, y = hU2OS.o$GMNN), aes(x = x, y = y)) +
	geom_point(shape = 1, size = 0.8, alpha = 0.6) +
	labs( x = str_c("FUCCI pseudotime"),
				y = bquote(paste('log'['2'],"(expression of GMNN)")),
				title = str_c(metadata(hU2OS.o)$dataname, " GMNN expression"))


loess.l <- fit_periodic_loess(theta.v = hU2OS.o$fucci_time.theta, y = hU2OS.o$GMNN)
hU2OS.GMNN.p <- plotLoess(sce.o = hU2OS.o, col.name = "GMNN", title =  bquote(paste('hU2OS log'['2'],"(GMNN)")), x_val = "fucci_time.theta") +
	annotate(geom = "text", x = pi, y = 2.5, size = 3, hjust = 0.5, vjust = 0.5,
					 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = str_c(seq(0, 1, 0.25)), limits = c(0, 2) * pi, name = "FUCCI pseudotime") 

loess.l <- fit_periodic_loess(theta.v = hU2OS.o$fucci_time.theta, y = hU2OS.o$CDT1)
hU2OS.CDT1.p <- plotLoess(sce.o = hU2OS.o, col.name = "CDT1", title =  bquote(paste('hU2OS log'['2'],"(CDT1)")), x_val = "fucci_time.theta") +
	annotate(geom = "text", x = pi, y = 1, size = 3, hjust = 0.5, vjust = 0.5,
					 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = str_c(seq(0, 1, 0.25)), limits = c(0, 2) * pi, name = "FUCCI pseudotime") 


loess.l <- fit_periodic_loess(theta.v = hU2OS.o$tricyclePosition, y = hU2OS.o$GMNN)
hU2OS.GMNN.theta.p <- plotLoess(sce.o = hU2OS.o, col.name = "GMNN", title =  bquote(paste('hU2OS log'['2'],"(GMNN)"))) +
	annotate(geom = "text", x = pi, y = 2.5, size = 3, hjust = 0.5, vjust = 0.5,
					 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE)

loess.l <- fit_periodic_loess(theta.v = hU2OS.o$tricyclePosition, y = hU2OS.o$CDT1)
hU2OS.CDT1.theta.p <- plotLoess(sce.o = hU2OS.o, col.name = "CDT1", title =  bquote(paste('hU2OS log'['2'],"(CDT1)"))) +
	annotate(geom = "text", x = pi, y = 1, size = 3, hjust = 0.5, vjust = 0.5,
					 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE)



mp <- plot_grid(hU2OS.GMNN.p + theme(legend.position = c(1, 0),
																		 legend.justification = c(1, 0),
																		 legend.key = element_blank(),
																		 legend.key.size = unit(7, "pt")),
								hU2OS.CDT1.p+ theme(legend.position = "none"),
								hU2OS.GMNN.theta.p + theme(legend.position = "none"),
								hU2OS.CDT1.theta.p + theme(legend.position = "none"),
												 nrow = 2, ncol = 2, label_size = 10, labels = "auto", align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.gmnncdti.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 2, device = cairo_pdf)
