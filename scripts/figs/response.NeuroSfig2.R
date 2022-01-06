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
x_lab <- str_c("PC1", sprintf(" (%i%%)", round(attr(reducedDim(neurosphere.o, "go.pca"), "percentVar")[1])))
y_lab <- str_c("PC2", sprintf(" (%i%%)", round(attr(reducedDim(neurosphere.o, "go.pca"), "percentVar")[2])))

a.p <- plotEmbScat(sce.o = neurosphere.o, dimred = "go.pca",  x_lab = x_lab, y_lab = y_lab,
						color_by = "CCStage", colors.v = ccColors.v,
						labels.v = ccLabels.v, color.name = 'SchwabeCC') +
	theme(legend.position = c(1, 0),
					legend.justification = c(1, 0))


b.p <- plotEmbScat(sce.o = neurosphere.o, dimred = "go.pca",  x_lab = x_lab, y_lab = y_lab,
						color_by = "SeuratCC", colors.v = cc3Colors.v,
						labels.v = cc3Labels.v, color.name = 'SeuratCC') +
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0))


mp <- plot_grid(
	a.p, b.p,
	nrow = 1, ncol = 2, label_size = 10, labels = c("a", "b"), align = "hv", axis = "tblr")

save_plot(here::here("figs", "response_figs", "response.NeuroSfig2.pdf"), mp,
					base_height = 2, base_width = 2*1.5, nrow = 1.2, ncol = 2, device = cairo_pdf)

