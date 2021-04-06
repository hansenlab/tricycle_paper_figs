rm(list=ls())
source(here::here("scripts/utils.R"))

library(ggrepel)

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))

neurosphere_full.o <- qread(here::here("data/neurosphere.qs"))
hipp_full.o <- qread(here::here("data/hipp.qs"))


getPeakPosition <- function(sce1.o, sce2.o, col.name, theta1.v, theta2.v, log2.trans = FALSE) {
	y1 <- colData(sce1.o)[, col.name]
	y2 <- colData(sce2.o)[, col.name]
	if (log2.trans) {
		y1 <- log2(y1)
		y2 <- log2(y2)
	}
	loess1.l <- fit_periodic_loess(theta.v = theta1.v, y = y1)
	loess2.l <- fit_periodic_loess(theta.v = theta2.v, y = y2)
	pred.df <- rbind(loess1.l$pred.df %>% add_column(source = "data1"), loess2.l$pred.df %>% add_column(source = "data2"))
	peak1 <- loess1.l$pred.df$x[which.max(loess1.l$pred.df$y)]
	peak2 <- loess2.l$pred.df$x[which.max(loess2.l$pred.df$y)]
	return(c(peak1, peak2))
}

### 
rotation_neurosphere.m <- metadata(neurosphere.o)$rotation
rotation_hipp.m <- metadata(hipp.o)$rotation
int_rotaGenes.v <- intersect(rownames(rotation_neurosphere.m), rownames(rotation_hipp.m))

### plot all selgenes
SelGenes.v <- int_rotaGenes.v[which((abs(rotation_neurosphere.m[int_rotaGenes.v, 1]) > 0.1) | (abs(rotation_neurosphere.m[int_rotaGenes.v, 2]) > 0.1) | (abs(rotation_hipp.m[int_rotaGenes.v, 1]) > 0.1) | (abs(rotation_hipp.m[int_rotaGenes.v, 2]) > 0.1))]
names(SelGenes.v) <- rowData(neurosphere.o)$Gene[match(SelGenes.v, rownames(neurosphere.o))]
names(SelGenes.v)[1] <- "Top2A"

pca.peak.m <- t(sapply(seq_len(length(SelGenes.v)), function(j) {
	colData(hipp.o)[, names(SelGenes.v)[j]] <- assay(hipp_full.o, "log.s")[match(SelGenes.v[j], rownames(hipp.o)), ]
	colData(neurosphere.o)[, names(SelGenes.v)[j]] <- assay(neurosphere_full.o, "log.s")[match(SelGenes.v[j], rownames(neurosphere.o)), ]
	
	out <- getPeakPosition(hipp.o, neurosphere.o, col.name =  names(SelGenes.v)[j],
												 theta1.v = as.numeric(coord2rad(reducedDim(hipp.o, "go.pca")[, 1:2])),
												 theta2.v = as.numeric(coord2rad(reducedDim(neurosphere.o, "go.pca")[, 1:2])))
	
	out
	
}))

projection.peak.m <- t(sapply(seq_len(length(SelGenes.v)), function(j) {
	colData(hipp.o)[, names(SelGenes.v)[j]] <- assay(hipp_full.o, "log.s")[match(SelGenes.v[j], rownames(hipp.o)), ]
	colData(neurosphere.o)[, names(SelGenes.v)[j]] <- assay(neurosphere_full.o, "log.s")[match(SelGenes.v[j], rownames(neurosphere.o)), ]
	
	out <- getPeakPosition(hipp.o, neurosphere.o, col.name =  names(SelGenes.v)[j],
												 theta1.v = hipp.o$tricyclePosition,
												 theta2.v = neurosphere.o$tricyclePosition)
	
	out
	
}))



tmp.df <- data.frame(x = pca.peak.m[, 2], y = pca.peak.m[, 1], label = names(SelGenes.v))

pca.peak.p <- ggplot(tmp.df, aes(x = x, y = y, label = label)) +
	# geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	# geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_abline(slope = 1, intercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point(size = 0.8,  alpha = 0.6,  shape = 1) +
	geom_text_repel(size = 1.3, segment.size = 0.2, max.overlaps = 30) +
	labs(title = str_c("Exp. peak positions on PCA \u03B8")) +
	scale_x_continuous(limits = c(0, 1.5 * pi),
										 breaks =  c(0, pi / 2, pi, 3 * pi / 2),
										 labels = c(0, str_c(c(0.5, 1, 1.5), "\u03C0")),
										 name = bquote(atop("Exp. peak positions of", paste(' mNeurosphere data on PCA \u03B8')))) +
	scale_y_continuous(limits = c(0, 1.5 * pi),
										 breaks =  c(0, pi / 2, pi, 3 * pi / 2),
										 labels = c(0, str_c(c(0.5, 1, 1.5), "\u03C0")),
										 name = bquote(atop("Exp. peak positions of", paste(' mHippNPC data on PCA \u03B8'))))





tmp.df <- data.frame(x = projection.peak.m[, 2], y = projection.peak.m[, 1], label = names(SelGenes.v))

projection.peak.p <- ggplot(tmp.df, aes(x = x, y = y, label = label)) +
	# geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	# geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_abline(slope = 1, intercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point(size = 0.8,  alpha = 0.6,  shape = 1) +
	geom_text_repel(size = 1.3, segment.size = 0.2, max.overlaps = 30) +
	labs(title = bquote( paste('Exp. peak positions on CC'['ns']," Position \u03B8")) ) +
	scale_x_continuous(limits = c(0, 1.5 * pi),
										 breaks =  c(0, pi / 2, pi, 3 * pi / 2),
										 labels = c(0, str_c(c(0.5, 1, 1.5), "\u03C0")),
										 name = bquote(atop("Exp. peak positions of", paste(' mNeurosphere data on CC'['ns']," Position \u03B8")))) +
	scale_y_continuous(limits = c(0, 1.5 * pi),
										 breaks =  c(0, pi / 2, pi, 3 * pi / 2),
										 labels = c(0, str_c(c(0.5, 1, 1.5), "\u03C0")),
										 name = bquote(atop("Exp. peak positions of", paste(' mHippNPC data on CC'['ns']," Position \u03B8"))))





tmp.df <- data.frame(x = abs(pca.peak.m[, 2]- pca.peak.m[, 1]), y = abs(projection.peak.m[, 2] - projection.peak.m[, 1]), label = names(SelGenes.v))
tmp.df$color <- factor(ifelse(tmp.df$x - tmp.df$y < 0, 1, 2))

distance.scat.p <- ggplot(tmp.df, aes(x = x, y = y, label = label, color = color)) +
	# geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	# geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_abline(slope = 1, intercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point(size = 0.8,  alpha = 0.6,  shape = 1) +
	geom_text_repel(size = 1.3, segment.size = 0.2, max.overlaps = 25) +
	scale_color_brewer(palette = "Set1", guide =  FALSE) + 
	labs(title = str_c("Absolute distance of peaks between \n mNeurosphere and mHippNPC dataset")) +
	scale_x_continuous(limits = c(0, max(tmp.df$x, tmp.df$y)) ,
										 breaks =  c(0, 0.25, 0.5, 0.75) ,
										 labels = c(0, str_c(format(c(0.25, 0.5, 0.75) / pi, digits = 2), "\u03C0")),
										 name = "Absolute distance of peaks on PCA \u03B8") +
	scale_y_continuous(limits = c(0, max(tmp.df$x, tmp.df$y)) ,
										 breaks =  c(0, 0.25, 0.5, 0.75) ,
										 labels = c(0, str_c(format(c(0.25, 0.5, 0.75) / pi, digits = 2), "\u03C0")),
										 name = bquote(atop("Absolute distance of peaks", paste(' on CC'['ns']," Position \u03B8")))) 
mp <- plot_grid(pca.peak.p, projection.peak.p, distance.scat.p,
												 nrow = 1, ncol = 3, label_size = 10, labels = "auto", align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.peakposition.pdf"), mp,
					base_height = 2, base_width = 2*1.5, nrow = 1*1.4, ncol = 3, device = cairo_pdf)



