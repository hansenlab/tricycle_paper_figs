rm(list=ls())
source(here::here("scripts/utils.R"))

library(ggrepel)

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))



plotTrendfilter <- function(sce1.o, sce2.o, col.name, col.outname, pointsize = 1.1, log2.trans = FALSE, title = NULL , y_lab = NULL) {
	if (is.null(title)) title <- col.outname
	if (is.null(y_lab)) y_lab <- bquote(paste('log'['2'],'(expression)'))
	theta1.v <- as.numeric(coord2rad(reducedDim(sce1.o, "go.pca")[, 1:2]))
	theta2.v <- as.numeric(coord2rad(reducedDim(sce2.o, "go.pca")[, 1:2]))
	y1 <- colData(sce1.o)[, col.name]
	y2 <- colData(sce2.o)[, col.name]
	if (log2.trans) {
		y1 <- log2(y1)
		y2 <- log2(y2)
	}
	tmp.df <- rbind(data.frame(theta = theta1.v, cc = sce1.o$CCStage, source = "data1", exp = y1),
									data.frame(theta = theta2.v, cc = sce2.o$CCStage, source = "data2", exp = y2))
	tmp.df$cc <- fct_explicit_na(tmp.df$cc, "NA")
	
	loess1.l <- fit_trendfilter(theta.v = theta1.v, y = y1)
	loess2.l <- fit_trendfilter(theta.v = theta2.v, y = y2)
	pred.df <- rbind(loess1.l$pred.df %>% add_column(source = "data1"), loess2.l$pred.df %>% add_column(source = "data2"))
	
	p <- ggplot(tmp.df, aes(x = theta, y = exp, color = source, group = source)) +
		geom_scattermore(pointsize = pointsize, alpha = 0.6) +
		geom_path(data = pred.df, aes(x = x , y = y, linetype = source), color = "black", size = 0.6, alpha = 0.6, inherit.aes = FALSE) +
		scale_color_manual(values = c("#F5C710",  "#2297E6"), name = "Dataset", label =c(metadata(sce1.o)$dataname, metadata(sce2.o)$dataname)) +
		scale_linetype_manual(name = "Dataset-loess", values=c("dotdash", "solid"), label = c(metadata(sce1.o)$dataname, metadata(sce2.o)$dataname)) +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( y = y_lab, x = "PCA \u03B8") +
		ggtitle(title) +
		ylim(range(tmp.df$exp)) +
		scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi) 
	# annotate(geom = "text", x = 0, y = .percent_range(tmp.df$exp, 1), size = 2, hjust = 0, vjust = 1,
	# 				 label = as.character(as.expression(substitute(dataname~" "~~italic(R)^2~"="~rsquared, list(dataname = metadata(sce2.o)$dataname,rsquared = format(loess2.l$rsquared, digits = 3))))), parse = TRUE) + 
	# annotate(geom = "text", x = 0, y = .percent_range(tmp.df$exp, 0.9), size = 2, hjust = 0, vjust = 1,
	# 				 label = as.character(as.expression(substitute(dataname~" "~~italic(R)^2~"="~rsquared, list(dataname = metadata(sce1.o)$dataname,rsquared = format(loess1.l$rsquared, digits = 3))))), parse = TRUE) 
	# 
	p
}

fit_trendfilter <- function (theta.v, y, length.out = 200) 
{
	stopifnot(`theta.v need to be between 0 - 2pi.` = (min(theta.v) >= 
																										 	0) & (max(theta.v) <= 2 * pi), `The length of theta.v and y should be the same.` = length(theta.v) == 
							length(y))
	
	ss.total <- sum(scale(y, scale = FALSE)^2)
	
	od <- order(theta.v)
	x <- c(theta.v - 2 * pi, theta.v, theta.v + 2 * pi)
	
	yy <- y[od]
	yy.rep <- rep(yy, 3)
	include <- rep(c(FALSE, TRUE, FALSE), each = length(yy))
	fit.trend <- trendfilter(yy.rep, ord = 2, approx = FALSE, 
													 maxsteps = 2000)
	cv.trend <- cv.trendfilter(fit.trend)
	which.lambda <- cv.trend$i.1se
	yy.trend.pred <- predict(fit.trend, lambda = cv.trend$lambda.1se, 
													 df = fit.trend$df[which.lambda])$fit
	trend.yy <- yy.trend.pred[, 1]
	fitted.v <- trend.yy[include]
	residual.v <- fitted.v - yy
	rsquared <- 1 - var(yy - yy.trend.pred)/var(yy)
	
	fun_g <- approxfun(x = as.numeric(theta.v[od]), 
										 y = as.numeric(fitted.v), rule = 2)
	
	pred.x <- seq(0, 2 * pi, length.out = length.out)
	pred.y <- fun_g(pred.x)
	pred.df <- data.frame(x = pred.x, y = pred.y)

	return(list(fitted = fitted.v, residual = residual.v, pred.df = pred.df, 
							fit.o = fit.trend, rsquared = rsquared))
}

top2a_t.p <- plotTrendfilter(hipp.o, neurosphere.o, col.name = "top2a", col.outname = "Top2A", title = bquote(italic("Top2A")))
smc4_t.p <- plotTrendfilter(hipp.o, neurosphere.o, col.name = "smc4", col.outname = "Smc4", title = bquote(italic("Smc4")))
umis_t.p <- plotTrendfilter(hipp.o, neurosphere.o, col.name = "TotalUMIs", col.outname = "TotalUMIs", log2.trans = TRUE, y_lab = bquote(paste('log'['2'],'(TotalUMIs)')))


plotLoess2 <- function(sce1.o, sce2.o, col.name, col.outname, pointsize = 1.1, log2.trans = FALSE, title = NULL , y_lab = NULL) {
	if (is.null(title)) title <- col.outname
	if (is.null(y_lab)) y_lab <- bquote(paste('log'['2'],'(expression)'))
	theta1.v <- as.numeric(coord2rad(reducedDim(sce1.o, "go.pca")[, 1:2]))
	theta2.v <- as.numeric(coord2rad(reducedDim(sce2.o, "go.pca")[, 1:2]))
	y1 <- colData(sce1.o)[, col.name]
	y2 <- colData(sce2.o)[, col.name]
	if (log2.trans) {
		y1 <- log2(y1)
		y2 <- log2(y2)
	}
	tmp.df <- rbind(data.frame(theta = theta1.v, cc = sce1.o$CCStage, source = "data1", exp = y1),
									data.frame(theta = theta2.v, cc = sce2.o$CCStage, source = "data2", exp = y2))
	tmp.df$cc <- fct_explicit_na(tmp.df$cc, "NA")
	
	loess1.l <- fit_periodic_loess(theta.v = theta1.v, y = y1)
	loess2.l <- fit_periodic_loess(theta.v = theta2.v, y = y2)
	pred.df <- rbind(loess1.l$pred.df %>% add_column(source = "data1"), loess2.l$pred.df %>% add_column(source = "data2"))
	
	p <- ggplot(tmp.df, aes(x = theta, y = exp, color = source, group = source)) +
		geom_scattermore(pointsize = pointsize, alpha = 0.6) +
		geom_path(data = pred.df, aes(x = x , y = y, linetype = source), color = "black", size = 0.6, alpha = 0.6, inherit.aes = FALSE) +
		scale_color_manual(values = c("#F5C710",  "#2297E6"), name = "Dataset", label =c(metadata(sce1.o)$dataname, metadata(sce2.o)$dataname)) +
		scale_linetype_manual(name = "Dataset-loess", values=c("dotdash", "solid"), label = c(metadata(sce1.o)$dataname, metadata(sce2.o)$dataname)) +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( y = y_lab, x = "PCA \u03B8") +
		ggtitle(title) +
		ylim(range(tmp.df$exp)) +
		scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi) 
	# annotate(geom = "text", x = 0, y = .percent_range(tmp.df$exp, 1), size = 2, hjust = 0, vjust = 1,
	# 				 label = as.character(as.expression(substitute(dataname~" "~~italic(R)^2~"="~rsquared, list(dataname = metadata(sce2.o)$dataname,rsquared = format(loess2.l$rsquared, digits = 3))))), parse = TRUE) + 
	# annotate(geom = "text", x = 0, y = .percent_range(tmp.df$exp, 0.9), size = 2, hjust = 0, vjust = 1,
	# 				 label = as.character(as.expression(substitute(dataname~" "~~italic(R)^2~"="~rsquared, list(dataname = metadata(sce1.o)$dataname,rsquared = format(loess1.l$rsquared, digits = 3))))), parse = TRUE) 
	# 
	p
}


top2a.p <- plotLoess2(hipp.o, neurosphere.o, col.name = "top2a", col.outname = "Top2A", title = bquote(italic("Top2A")))
smc4.p <- plotLoess2(hipp.o, neurosphere.o, col.name = "smc4", col.outname = "Smc4", title = bquote(italic("Smc4")))
umis.p <- plotLoess2(hipp.o, neurosphere.o, col.name = "TotalUMIs", col.outname = "TotalUMIs", log2.trans = TRUE, y_lab = bquote(paste('log'['2'],'(TotalUMIs)')))


mp <- plot_grid(top2a.p + theme(legend.position = c(1, 1),
																legend.justification = c(1, 1),
																legend.key = element_blank(),
																legend.key.size = unit(6.5, "pt"),
																legend.spacing.y = unit(-0.5, "pt"),
																axis.title.y = element_text(size = 8),
																axis.title.x = element_text(size = 8)),
								smc4.p + theme(legend.position = "none", axis.title.y = element_text(size = 8),
															 axis.title.x = element_text(size = 8)),
								umis.p + theme(legend.position = "none", axis.title.y = element_text(size = 8),
															 axis.title.x = element_text(size = 8)),
								top2a_t.p + theme(legend.position = "none", axis.title.y = element_text(size = 8),
																	axis.title.x = element_text(size = 8)),
								smc4_t.p + theme(legend.position = "none", axis.title.y = element_text(size = 8),
															 axis.title.x = element_text(size = 8)),
								umis_t.p + theme(legend.position = "none", axis.title.y = element_text(size = 8),
															 axis.title.x = element_text(size = 8)),
								nrow = 2, ncol = 3,  label_size = 10, labels = c("auto"), align = "h", axis = "l")


save_plot(here::here("figs", "response_figs", "response.trendfilter.pdf"), mp,
					base_height = 2 / 1.2, base_width = 2*1.2 / 1.2, nrow = 2, ncol = 3, device = cairo_pdf)


