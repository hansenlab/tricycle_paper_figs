rm(list=ls())
source(here::here("scripts/utils.R"))

library(ggrepel)
hU2OS.o <- qread(here::here("data/plotdata/hU2OS.qs"))
hiPSCs.o <- qread(here::here("data/plotdata/hiPSCs.qs"))

hU2OS_full.o <- qread(here::here("data/hU2OS.qs"))
hiPSCs_full.o <- qread(here::here("data/hiPSCs.qs"))

hU2OS.fucci.p <- plotEmbScatCyclic(hU2OS.o, dimred = "fucci",x_lab = "GMNN-Green530", y_lab = "CDT1-Red585",
									title = str_c( metadata(hU2OS.o)$dataname, " FUCCI scores (n=", ncol(hU2OS.o), ")"))


### move fucci time
theta.v <- hU2OS.o$tricyclePosition
r.theta.v <- theta.v + 0.85 * pi
r.theta.v[r.theta.v > 2 * pi] <- r.theta.v[r.theta.v > 2 * pi] - 2 * pi
r.theta.v[(r.theta.v < 0.25 * pi) & (hU2OS.o$fucci_time > 0.5)] <- r.theta.v[(r.theta.v < 0.25 * pi) & (hU2OS.o$fucci_time > 0.5)] + 2 * pi

hU2OS.thetapseudotime.p <- ggplot(data.frame(PC1 =  r.theta.v , PC2 = hU2OS.o$fucci_time, type = theta.v), aes(x = PC1, y = PC2 , color = type)) +
	annotate("rect", xmin = 0, xmax = 0.25 * pi, ymin = 0.5, ymax = 1, alpha = .6, fill = NA, linetype = 'dotted', color = "black") +
	annotate("rect", xmin = 2 * pi, xmax = 2.25 * pi, ymin = 0.5, ymax = 1, alpha = .6, fill = NA, linetype = 'dotted', color = "black") +
	annotate("segment", x = 0.5 * pi, xend = 1.25 * pi, y = 0.6, yend = 0.6, colour = "black", size = 0.5, alpha = 0.6, arrow=arrow(length = unit(0.1, "inches")),  linetype = 'solid') + 
	geom_scattermore(pointsize  =metadata(hU2OS.o)$point.size,  alpha = 0.6) +
	geom_smooth(method = "loess", color = "black", se = FALSE, linetype = "dashed", size = 0.8, alpha = 0.8) +
	scale_color_gradientn(name = "Cell cycle position", limits = range(0, 2 * pi), 
												breaks = seq(from = 0, to = 2 * pi, length.out = 500) ,
												colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"), 
												guide = FALSE) +
	labs( y = str_c("FUCCI pseudotime"),
				x = theta_lab,
				title = str_c(metadata(hU2OS.o)$dataname, " FUCCI pseudotime ~ \u03B8  (n=", ncol(hU2OS.o), ")")) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = str_c( (c(0, 0.5, 1, 1.5, 2) + 0.85 - c(0, 0, 0, 2, 2)), "\u03C0"))


### calculate R2 for hU2OS
neuro_ref <- tricycle::neuroRef
genes <- intersect(rowData(hU2OS_full.o)$Gene, neuro_ref$SYMBOL)
hU2OS_full.o <- hU2OS_full.o[match(genes, rowData(hU2OS_full.o)$Gene), ]

hU2OSR2.m <- t(sapply(seq_len(nrow(hU2OS_full.o)), function(i) {
	theta.r2 <- fit_periodic_loess(theta.v = hU2OS_full.o$tricyclePosition, y = assay(hU2OS_full.o, "log.s")[i, ])$rsquared
	fucci.r2 <- fit_periodic_loess(theta.v = hU2OS_full.o$fucci_time * 2 * pi, y = assay(hU2OS_full.o, "log.s")[i, ])$rsquared
	return(c(theta.r2, fucci.r2))
}))
rownames(hU2OSR2.m) <- genes
idx <- which((hU2OSR2.m[, 1] > 0.6))
labels.df <- data.frame(hU2OSR2.m[idx, ], label = rownames(hU2OSR2.m)[idx])
names(labels.df) <- c("y", "x", "label")

hU2OS.r2.p <- ggplot(data.frame(y = hU2OSR2.m[, 1], x = hU2OSR2.m[, 2]), aes(x = x, y = y)) +
	geom_point(shape = 1, size = 0.6) +
	geom_abline(slope = 1, intercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_text_repel(data = labels.df, aes(x = x, y = y, label = label), size = 1.3, segment.size = 0.2, max.overlaps = 30) +
	labs( x = bquote(paste("FUCCI pseudotime R"^"2")),
				y = bquote(paste('CC'['ns']," Position \u03B8", " R"^"2")),
				title = bquote(paste( "hU2OS R"^"2"))) +
	xlim(range(hU2OSR2.m)) + ylim(range(hU2OSR2.m))











### hiPSCs figs
reducedDim(hiPSCs.o, "fucci") <- reducedDim(hiPSCs.o, "fucci")[, c(2, 1)]
hiPSCs.fucci.p <- plotEmbScatCyclic(hiPSCs.o, dimred = "fucci", y_lab = "CDT1-mCherry", x_lab = "GMNN-EGFP",
																		title = str_c( metadata(hiPSCs.o)$dataname, " FUCCI scores (n=", ncol(hiPSCs.o), ")"))

loess.l <- fit_periodic_loess(theta.v = hiPSCs.o$tricyclePosition, y = hiPSCs.o$top2a)
hiPSCs.top2a.p <- plotLoess(sce.o = hiPSCs.o, col.name = "top2a", title =  bquote(paste('hiPSCs log'['2'],"(TOP2A)"))) +
	annotate(geom = "text", x = pi, y = 2.5, size = 3, hjust = 0.5, vjust = 0.5,
					 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE)

loess.l <- fit_periodic_loess(theta.v = hiPSCs.o$Hsiao.theta, y = hiPSCs.o$top2a)
hiPSCs.top2a2.p <- plotLoess(sce.o = hiPSCs.o, col.name = "top2a", x_val = "Hsiao.theta", x_lab = "FUCCI pseudotime", title = bquote(paste('hiPSCs log'['2'],"(TOP2A)"))) +
	annotate(geom = "text", x = pi, y = 2.5, size = 3, hjust = 0.5, vjust = 0.5,
					 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE)



### calculate R2 for hiPSCs
genes <- intersect(rowData(hiPSCs_full.o)$Gene, neuro_ref$SYMBOL)
hiPSCs_full.o <- hiPSCs_full.o[match(genes, rowData(hiPSCs_full.o)$Gene), ]

hiPSCsR2.m <- t(sapply(seq_len(nrow(hiPSCs_full.o)), function(i) {
	theta.r2 <- fit_periodic_loess(theta.v = hiPSCs_full.o$tricyclePosition, y = assay(hiPSCs_full.o, "log.s")[i, ])$rsquared
	hsiao.r2 <- fit_periodic_loess(theta.v = hiPSCs.o$Hsiao.theta, y = assay(hiPSCs_full.o, "log.s")[i, ])$rsquared
	return(c(theta.r2, hsiao.r2))
}))
rownames(hiPSCsR2.m) <- genes
idx <- which((hiPSCsR2.m[, 1] > 0.2) | (hiPSCsR2.m[, 2] > 0.2))
labels.df <- data.frame(hiPSCsR2.m[idx, ], label = rownames(hiPSCsR2.m)[idx])
names(labels.df) <- c("y", "x", "label")

hiPSCs.r2.p <- ggplot(data.frame(y = hiPSCsR2.m[, 1], x = hiPSCsR2.m[, 2]), aes(x = x, y = y)) +
	geom_point(shape = 1, size = 0.6) +
	geom_abline(slope = 1, intercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_text_repel(data = labels.df, aes(x = x, y = y, label = label), size = 1.3, segment.size = 0.2, max.overlaps = 30) +
	labs( x = bquote(paste("FUCCI pseudotime R"^"2")),
				y = bquote(paste('CC'['ns']," Position \u03B8", " R"^"2")),
				title = bquote(paste( "hiPSCs R"^"2"))) +
	xlim(range(hiPSCsR2.m)) + ylim(range(hiPSCsR2.m))





mp <- plot_grid(hU2OS.fucci.p + theme(legend.position = "none"),  
								
								hU2OS.thetapseudotime.p + theme(legend.position = "none"), 
								hU2OS.r2.p,
								circle_scale_legend(),
								hiPSCs.fucci.p + theme(legend.position = "none"),
								hiPSCs.top2a.p + theme(legend.position = c(1, 0),
																			 legend.justification = c(1, 0),
																			 legend.key = element_blank(),
																			 legend.key.size = unit(6, "pt")),
								hiPSCs.top2a2.p + theme(legend.position = "none"),
								hiPSCs.r2.p,
								# hU2OS.thetaTOP2A.p + theme(legend.position = "none"), 
								
								nrow = 2, ncol = 4, label_size = 10, labels = c(letters[1:3], "", "d", "e", "f", "g"), align = "hv", axis = "tblr")

save_plot(here::here("figs", "main", "main.fucci.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 4, device = cairo_pdf)

