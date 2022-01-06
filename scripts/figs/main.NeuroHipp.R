rm(list=ls())
source(here::here("scripts/utils.R"))

library(ggrepel)

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))


getCirclePlot <- function(sce.o, dimred = "go.pca", r, label.x, label.y, label.size = 3) {
	
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	title <- str_c( metadata(sce.o)$dataname, " (n=", ncol(sce.o), ")")
	emb.m <- reducedDim(sce.o, dimred)
	
	tmp.df <- data.frame(pc1 = emb.m[, 1], pc2 = emb.m[, 2]) %>% add_column(cc = sce.o$CCStage)
	
	x_lab <- str_c("PC1", sprintf(" (%i%%)", round(attr(emb.m, "percentVar")[1])))
	y_lab <- str_c("PC2", sprintf(" (%i%%)", round(attr(emb.m, "percentVar")[2])))
	xlim <- c(min(-r, min(tmp.df$pc1)), max(r, max(tmp.df$pc1)))
	ylim <- c(min(-r, min(tmp.df$pc2)), max(r, max(tmp.df$pc2)))
	
	tmp.df$cc <- fct_explicit_na(tmp.df$cc, na_level = "NA") %>% fct_relevel("NA", after = Inf)
	
	cc.p <- ggplot(tmp.df, aes(x = pc1, y = pc2 , color = cc)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`cc` == "NA"), pointsize = point.size,  alpha = point.alpha) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`cc` != "NA"), pointsize = point.size,  alpha = point.alpha) +
		scale_color_manual(values = c(ccColors.v, "grey"), name = "SchwabeCC", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA")) +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		annotate("point", x = 0, 
						 y = 0, shape = 16, size = 0.8, alpha = 0.8) + 
		annotate("segment", x = 0, xend =  r ,
						 y = 0, yend = 0, linetype = "dashed", alpha = 0.8) + 
		annotate("segment", x = 0, xend =  cos(pi/6) * r,
						 y = 0, yend =  sin(pi/6) * r,
						 linetype = "dashed", alpha = 0.8) +
		# annotate("path", linetype = "dashed", alpha = 0.8,
		# 				 x = cos(seq(0, 2*pi, length.out=100)) * r,
		# 				 y = sin(seq(0, 2*pi, length.out=100)) * r) +
		# annotate("path", linetype = "solid", alpha = 0.8,
		# 				 x = cos(seq(0, pi/6, length.out=100))/6 ,
		# 				 y = sin(seq(0, pi/6, length.out=100))/6) +
		annotate("text",  alpha = 0.8, x = label.x, y = label.y, label = "\u03B8", parse = TRUE, size = 3) +
		xlim(xlim) + ylim(ylim) +
		labs( y = y_lab, 
					x = x_lab, 
					title = title) 
	
	return(cc.p)
}

### panel a and b
pca_neurosphere.p <- getCirclePlot(neurosphere.o,  r = 10, label.x = 8, label.y = 2, label.size = 3)
pca_hipp.p <- getCirclePlot(hipp.o, r = 10, label.x = 8, label.y = 2, label.size = 3)


### compare rotation - panel c
rotation_neurosphere.m <- metadata(neurosphere.o)$rotation
rotation_hipp.m <- metadata(hipp.o)$rotation
int_rotaGenes.v <- intersect(rownames(rotation_neurosphere.m), rownames(rotation_hipp.m))
tmp.df <- data.frame(pc1.neu = rotation_neurosphere.m[int_rotaGenes.v, 1],
										 pc2.neu = rotation_neurosphere.m[int_rotaGenes.v, 2],
										 pc1.hipp = rotation_hipp.m[int_rotaGenes.v, 1],
										 pc2.hipp = rotation_hipp.m[int_rotaGenes.v, 2])
tmp.df$Gene <- rowData(neurosphere.o)$Gene[match(int_rotaGenes.v, rownames(neurosphere.o))]
tmp.df %<>% mutate(color1 = (abs(pc1.neu) > 0.1) | (abs(pc1.hipp) > 0.1))
tmp.df %<>% mutate(color2 = (abs(pc2.neu) > 0.1) | (abs(pc2.hipp) > 0.1))
tmp1.df <- tmp.df %>% dplyr::filter((abs(pc1.neu) > 0.12) & (abs(pc1.hipp) > 0.15))
tmp2.df <- tmp.df %>% dplyr::filter((abs(pc2.neu) > 0.1) | (abs(pc2.hipp) > 0.1))


pc1.scat.p <- ggplot(tmp.df, aes(x = pc1.neu, y = pc1.hipp, label = Gene)) +
	geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_abline(slope = 1, intercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point(aes(color = color1), size = 0.5,  alpha = 0.4,  shape = 1) +
	geom_text_repel(data = tmp1.df, size = 2, segment.size = 0.2, max.overlaps = 40) +
	scale_color_manual(values = c("black", "red"), name = " ", guide = FALSE) +
	labs( y = "mHippNPC PC1 weights", x = "mNeurosphere PC1 weights", title = str_c("Overlapped genes(n=", length(int_rotaGenes.v), ")")) +
	annotate(geom = "text", x = -0.2, y = 0.05, hjust = 0, vjust = 1,  label = str_c("PCC == ", format(cor(tmp.df$pc1.neu, tmp.df$pc1.hipp), digits = 2)), parse = TRUE, size = 3)




### fit two genes
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


mp <- plot_grid(pca_neurosphere.p + ylim(c(-12.5, max(reducedDim(neurosphere.o, "go.pca")[, 2]))) +
									guides(color = guide_legend(override.aes = list(alpha = 1, size = 0.5), nrow= 2,byrow=TRUE)) +
									theme(legend.position = c(1, 0),
																								legend.justification = c(1, 0),
																								legend.key = element_blank(),
																								legend.key.size = unit(6, "pt"),
																					legend.title = element_text(size = 6),
																					legend.text = element_text(size = 6),
																					axis.title.y = element_text(size = 8),
																					axis.title.x = element_text(size = 8)) ,
								pca_hipp.p + theme(legend.position = "none", axis.title.y = element_text(size = 8),
																	 axis.title.x = element_text(size = 8)), 
								pc1.scat.p + theme(axis.title.y = element_text(size = 8),
																	 axis.title.x = element_text(size = 8)),
								top2a.p + theme(legend.position = c(1, 1),
																legend.justification = c(1, 1),
																legend.key = element_blank(),
																legend.key.size = unit(6.5, "pt"),
																legend.spacing.y = unit(-0.5, "pt"),
																axis.title.y = element_text(size = 8),
																axis.title.x = element_text(size = 8),
																legend.title = element_text(size = 6),
																legend.text = element_text(size = 6)),
								smc4.p + theme(legend.position = "none", axis.title.y = element_text(size = 8),
															 axis.title.x = element_text(size = 8)),
								umis.p + theme(legend.position = "none", axis.title.y = element_text(size = 8),
															 axis.title.x = element_text(size = 8)),
										 nrow = 2, ncol = 3,  label_size = 10, labels = c("auto"), align = "h", axis = "l")


save_plot(here::here("figs", "main", "main.NeuroHipp.pdf"), mp,
					base_height = 2 / 1.2, base_width = 2*1.2 / 1.2, nrow = 2, ncol = 3, device = cairo_pdf)




