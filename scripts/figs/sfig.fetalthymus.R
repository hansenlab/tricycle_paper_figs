rm(list=ls())
source(here::here("scripts/utils.R"))


o <- "Thymus"

sce.o <- qread(here::here("data/fetal/nuclei", str_c(o, ".qs")))
point.size <- ifelse(ncol(sce.o) > 1000000, 1.01, 2.01)
sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
reducedDim(sce.o, "umap") <- data.frame(umap1 = sce.o$Main_cluster_umap_1, umap2 = sce.o$Main_cluster_umap_2)
p1 <- plot_emb_circle_scale(sce.o, dimred = 2, fig.title = str_c("Human fetal ",o, " (n=", ncol(sce.o), ")"), x_lab = "UMAP-1", y_lab = "UMAP-2")

tmp.df <- data.frame(x = reducedDim(sce.o, "umap")[, 1], y = reducedDim(sce.o, "umap")[, 2], color = factor(sce.o$cc))
tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
scale_color <- scale_color_manual(values = c(ccColors.v, "grey"), name = "SchwabeCC", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA"))

p2 <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
	geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = 2.01,  alpha = 0.6) +
	geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = 2.01,  alpha = 0.6) +
	scale_color +
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
	labs( x = "UMAP-1", y = "UMAP-2", title = str_c("Human fetal ", o, " (n=", ncol(sce.o), ")")) +
	theme(legend.position = c(0, 1),
			legend.justification = c(0, 1),
			legend.key = element_blank(),
			legend.key.size = unit(6.5, "pt"))

### get cyclone and seurat
sce.o$cyclone <- runCyclone(sce.o, gname = rownames(sce.o), assay.type = "counts", species = "human", id.type =  "name")
colnames(sce.o) <- sce.o$sample
sce.o$SeuratCC <- runSeuratCC(sce.o, gname = rownames(sce.o), assay.type = "counts", species = "human", id.type =  "name")

tmp.df$cyclone <- sce.o$cyclone
tmp.df$SeuratCC <- sce.o$SeuratCC

p3 <- ggplot(tmp.df, aes(x = x, y = y, color = cyclone)) +
	geom_scattermore(pointsize = 2.01,  alpha = 0.6) +
	scale_color_manual(values = cc3Colors.v, name = "cyclone", labels =  cc3Labels.v, limits =   cc3Labels.v) +
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
	labs( x = "UMAP-1", y = "UMAP-2", title = str_c("Human fetal ", o, " (n=", ncol(sce.o), ")")) +
	theme(legend.position = c(0, 1),
				legend.justification = c(0, 1),
				legend.key = element_blank(),
				legend.key.size = unit(6.5, "pt"))

p4 <- ggplot(tmp.df, aes(x = x, y = y, color = SeuratCC)) +
	geom_scattermore(pointsize = 2.01,  alpha = 0.6) +
	scale_color_manual(values = cc3Colors.v, name = "SeuratCC", labels =  cc3Labels.v, limits =   cc3Labels.v) +
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
	labs( x = "UMAP-1", y = "UMAP-2", title = str_c("Human fetal ", o, " (n=", ncol(sce.o), ")")) +
	theme(legend.position = c(0, 1),
				legend.justification = c(0, 1),
				legend.key = element_blank(),
				legend.key.size = unit(6.5, "pt"))





mp <- plot_grid(p1, p2, p3, p4, 
								nrow = 2, ncol = 2,label_size = 10, labels = "auto", align = "v", axis = "lr")


save_plot(here::here("figs", "sfigs", "sfig.fetalthymus.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 2, device = cairo_pdf)

