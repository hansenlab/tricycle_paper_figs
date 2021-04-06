rm(list=ls())
source(here::here("scripts/utils.R"))

library(ggrepel)

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))

neurosphere_full.o <- qread(here::here("data/neurosphere.qs"))
hipp_full.o <- qread(here::here("data/hipp.qs"))



### fit two genes -  this function is not the same as in main.NeuroHipp.R
plotLoess2 <- function(sce1.o, sce2.o, col.name, col.outname, pointsize = 1.1, log2.trans = FALSE, title = NULL ) {
	if (is.null(title)) title <- col.outname
	theta1.v <- sce1.o$tricyclePosition
	theta2.v <- sce2.o$tricyclePosition
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
		labs( y = bquote(paste('log'['2'],'(expression)')), x = theta_lab) +
		ggtitle(title) +
		ylim(range(tmp.df$exp)) +
		scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2) * pi) +
		annotate(geom = "text", x = 0, y = .percent_range(tmp.df$exp, 1), size = 2, hjust = 0, vjust = 1,
						 label = as.character(as.expression(substitute(dataname~" "~~italic(R)^2~"="~rsquared, list(dataname = metadata(sce2.o)$dataname,rsquared = format(loess2.l$rsquared, digits = 3))))), parse = TRUE) +
		annotate(geom = "text", x = 0, y = .percent_range(tmp.df$exp, 0.9), size = 2, hjust = 0, vjust = 1,
						 label = as.character(as.expression(substitute(dataname~" "~~italic(R)^2~"="~rsquared, list(dataname = metadata(sce1.o)$dataname,rsquared = format(loess1.l$rsquared, digits = 3))))), parse = TRUE)
	
	p
}

### 
rotation_neurosphere.m <- metadata(neurosphere.o)$rotation
rotation_hipp.m <- metadata(hipp.o)$rotation
int_rotaGenes.v <- intersect(rownames(rotation_neurosphere.m), rownames(rotation_hipp.m))

### plot all selgenes
SelGenes.v <- int_rotaGenes.v[which((abs(rotation_neurosphere.m[int_rotaGenes.v, 1]) > 0.1) | (abs(rotation_neurosphere.m[int_rotaGenes.v, 2]) > 0.1) | (abs(rotation_hipp.m[int_rotaGenes.v, 1]) > 0.1) | (abs(rotation_hipp.m[int_rotaGenes.v, 2]) > 0.1))]
names(SelGenes.v) <- rowData(neurosphere.o)$Gene[match(SelGenes.v, rownames(neurosphere.o))]
names(SelGenes.v)[1] <- "Top2A"

SelGenes.lp <- lapply(seq_len(length(SelGenes.v)), function(j) {
	colData(hipp.o)[, names(SelGenes.v)[j]] <- assay(hipp_full.o, "log.s")[match(SelGenes.v[j], rownames(hipp.o)), ]
	colData(neurosphere.o)[, names(SelGenes.v)[j]] <- assay(neurosphere_full.o, "log.s")[match(SelGenes.v[j], rownames(neurosphere.o)), ]
	p <- plotLoess2(hipp.o, neurosphere.o, col.name =  names(SelGenes.v)[j], 
									col.outname =  names(SelGenes.v)[j], title = bquote(italic( .(names(SelGenes.v)[j]))))
	p
	
})

SelGenes.mp <- plot_grid(plotlist = c(list(SelGenes.lp[[1]] + theme(legend.position = c(1, 1),
																																		legend.justification = c(1, 1))), 
																			lapply(SelGenes.lp[2:length(SelGenes.lp)], function(x) x + theme(legend.position = "none"))),
												 nrow = 7, ncol = 5, label_size = 10, labels = NULL, align = "hv", axis = "tblr")

# save_plot(here::here("figs", "sfigs", "sfig.NeuroHippGenesTheta.pdf"), SelGenes.mp,
# 					base_height = 2, base_width = 2*1.2, nrow = 7, ncol = 5, device = cairo_pdf)
save_plot(here::here("figs", "sfigs", "sfig.NeuroHippGenesProjTheta.jpg"), SelGenes.mp,
					base_height = 2, base_width = 2*1.2, nrow = 7, ncol = 5, type = "cairo")
