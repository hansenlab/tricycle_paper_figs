rm(list=ls())
source(here::here("scripts/utils.R"))

library(randomcoloR)
colData.df <- qread(here::here("data", "plotdata", str_c("fetal_colData.qs"))) %>% as.data.frame()


colData.df$color <- colData.df$Main_cluster_name
sel <- names(sort(table(colData.df$Main_cluster_name), decreasing = TRUE)[1:50])
colData.df$color[!(colData.df$color %in% sel)] <- NA
colData.df$color <- factor(colData.df$color, levels = sel)


set.seed(100)
ct.p <- ggplot(colData.df, aes(x = Global_umap_1, y = Global_umap_2, color = color)) +
	geom_scattermore(pointsize = 0.3,  alpha = 0.6) +
	scale_color_manual(values = distinctColorPalette(nlevels(factor(colData.df$color))),
										 name = "Cell Type",
										 labels = str_c(levels(factor(colData.df$color)), " (n=", table(factor(colData.df$color)), ")"), 
										 limits = sel)+
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) + 
	labs( x = "UMAP-1", y = "UMAP-2", title = str_c( "Human fetal tissue atalas (n=", nrow(colData.df), ")")) +
	theme(		legend.key.size = unit(8, "pt"),
					legend.key = element_blank(),
					legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
					legend.title = element_text(size = 10),
					legend.title.align = 0.5,
					legend.text.align = 0.5,
					legend.text = element_text(size = 8))
	





colData.df$color <- colData.df$Organ
sel <- names(sort(table(colData.df$Organ), decreasing = TRUE))
colData.df$color[!(colData.df$color %in% sel)] <- NA
colData.df$color <- factor(colData.df$color, levels = sel)


set.seed(100)
Organ.p <- ggplot(colData.df, aes(x = Global_umap_1, y = Global_umap_2, color = color)) +
	geom_scattermore(pointsize = 0.3,  alpha = 0.6) +
	scale_color_manual(values = distinctColorPalette(nlevels(factor(colData.df$color))),
										 name = "Tissue",
										 labels = str_c(levels(factor(colData.df$color)), " (n=", table(factor(colData.df$color)), ")"), 
										 limits = sel)+
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 2), ncol = 3)) + 
	labs( x = "UMAP-1", y = "UMAP-2", title = str_c( "Human fetal tissue atalas (n=", nrow(colData.df), ")")) +
	theme(		legend.key.size = unit(8, "pt"),
					legend.key = element_blank(),
					legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
					legend.title = element_text(size = 10),
					legend.title.align = 0.5,
					legend.text.align = 0.5,
					legend.text = element_text(size = 8))





mp <- plot_grid(Organ.p + theme(legend.position = "none"),
								get_legend(Organ.p),
								ct.p + theme(legend.position = "none"),
								
								get_legend(ct.p), 
								nrow = 2, ncol = 2, rel_widths = c(1, 0.4), label_size = 10, labels = c("a", "","b"), align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.fetalglobalumap.pdf"), mp,
					base_height = 2, base_width = 2, nrow = 4.5 * , ncol = 4.5, device = cairo_pdf)












