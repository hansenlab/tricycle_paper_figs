rm(list=ls())
source(here::here("scripts/utils.R"))

endo.o <- qread(here::here("data/plotdata/endo.qs"))
endo_full.o <- qread(here::here("data/endo.qs"))

metadata(endo_full.o) <- metadata(endo.o)
### umap figure
tmp.df <- data.frame(pc1.s = reducedDim(endo.o, "umap")[, 1], pc2.s = reducedDim(endo.o, "umap")[, 2],
										 cell_type = endo.o$cell_type)

umap.scat.p <- ggplot(tmp.df, aes(x = pc1.s, y = pc2.s, color = cell_type)) +
	geom_scattermore(pointsize = metadata(endo.o)$point.size, alpha = metadata(endo.o)$point.alpha) +
	scale_color_manual(values = metadata(endo.o)$clusters_colors, name = "Cell type", labels = levels(factor(endo.o$cell_type))) +
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
	labs( y = "UMAP-2", x = "UMAP-1", title = "UMAP (mPancreas)") +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))


### GO PC
pca.celltype.p <- plotEmbScat(sce.o = endo.o, dimred = "go.pca", 
															x_lab = str_c("PC1", sprintf(" (%i%%)", round(attr(reducedDim(endo.o, "go.pca"), "percentVar")[1]))),
															y_lab = str_c("PC2", sprintf(" (%i%%)", round(attr(reducedDim(endo.o, "go.pca"), "percentVar")[2]))),
															color_by = "cell_type", facet_by = NULL, facet_labels = NULL, 
															colors.v = metadata(endo.o)$clusters_colors, color.name = "Cell type",
															labels.v = levels(factor(endo.o$cell_type)), title = "PCA of GO (mPancreas)")




### ductal
GO_Ductal.o <- getGO(endo_full.o[, endo_full.o$cell_type == "Ductal"], row.id = rownames(endo_full.o), id.type = "SYMBOL", runSeuratBy =  NULL)
plot(reducedDim(GO_Ductal.o, "PCA.s"))


full_weight.m <- attr(reducedDim(endo.o, "go.pca"), "rotation")
ductal_weight.m <- attr(reducedDim(GO_Ductal.o, "PCA.s"), "rotation")

genes.v <- intersect(rownames(full_weight.m), rownames(ductal_weight.m))

### plot scatter plot
tmp.df <- data.frame(full_pc1 = full_weight.m[genes.v, 1],
										 full_pc2 = full_weight.m[genes.v, 2],
										 ductal_pc1 = ductal_weight.m[genes.v, 1],
										 ductal_pc2 = ductal_weight.m[genes.v, 2],
										 gene = genes.v)

tmp1.df <- tmp.df %>% dplyr::filter((abs(full_pc1) > 0.15) | (abs(ductal_pc1) > 0.15))
tmp2.df <- tmp.df %>% dplyr::filter((abs(full_pc2) > 0.15) | (abs(ductal_pc2) > 0.15))

pc1.p <- ggplot(tmp.df, aes(x = full_pc1, y = ductal_pc1, label = gene)) +
	geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point( size = 0.5,  alpha = 0.4,  shape = 1) +
	geom_text_repel(data = tmp1.df, size = 1.3, segment.size = 0.2, max.overlaps = 25) +
	labs( y = "Ductal PC1 weights", x = "Full PC1 weights", title = str_c("PC1 weights of overlapped genes(n=", length(genes.v), ")")) 

pc2.p <- ggplot(tmp.df, aes(x = full_pc2, y = ductal_pc2, label = gene)) +
	geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point( size = 0.5,  alpha = 0.4,  shape = 1) +
	geom_text_repel(data = tmp2.df, size = 1.3, segment.size = 0.2, max.overlaps = 25) +
	labs( y = "Ductal PC2 weights", x = "Full PC2 weights", title = str_c("PC2 weights of overlapped genes(n=", length(genes.v), ")")) 







set.seed(500)
endo_full.o <- runPCA(endo_full.o, exprs_values = "log.s",  name = "newPCA", ncomponents = 20, subset_row = rownames(ductal_weight.m))
tmp <- attr(reducedDim(endo_full.o, "newPCA"), "percentVar")
reducedDim(endo_full.o, "newPCA")[, 1] <- - reducedDim(endo_full.o, "newPCA")[, 1]
attr(reducedDim(endo_full.o, "newPCA"), "percentVar") <- tmp

tmp <- attr(reducedDim(endo_full.o, "newPCA"), "rotation")
tmp[, 1] <- -tmp[, 1]
attr(reducedDim(endo_full.o, "newPCA"), "rotation") <- tmp 



newPCA.celltype.p <- plotEmbScat(sce.o = endo_full.o, dimred = "newPCA", 
																 x_lab = str_c("PC1", sprintf(" (%i%%)", round(attr(reducedDim(endo_full.o, "newPCA"), "percentVar")[1]))),
																 y_lab = str_c("PC2", sprintf(" (%i%%)", round(attr(reducedDim(endo_full.o, "newPCA"), "percentVar")[2]))),
																 color_by = "cell_type", facet_by = NULL, facet_labels = NULL, 
																 colors.v = metadata(endo_full.o)$clusters_colors, color.name = "Cell type",
																 labels.v = levels(factor(endo_full.o$cell_type)), title = str_c("PCA of ductal genes (mPancreas)"))

mp <- plot_grid(pca.celltype.p + theme(legend.position = "none"), 
								newPCA.celltype.p + theme(legend.position = "none"), 
								get_legend(pca.celltype.p),
								nrow = 1, ncol = 3, rel_widths = c(1, 1, 0.35), label_size = 10, labels = c("a", "b", ""), align = "hv", axis = "tblr")

save_plot(here::here("figs", "response_figs", "extra.endoProj2.pdf"), mp,
					base_height = 2 , base_width = 2*1.4, nrow = 1, ncol = 2, device = cairo_pdf)





genes.v <- union(rownames(full_weight.m), rownames(ductal_weight.m))
sd_full.v <- rowSds(as(assay(endo_full.o[genes.v, ], "log.s"), "dgCMatrix"))
sd_ductal.v <- rowSds(as(assay(endo_full.o[genes.v, endo_full.o$cell_type == "Ductal"], "log.s"), "dgCMatrix"))

a <- genes.v %in% rownames(full_weight.m)
b <- genes.v %in% rownames(ductal_weight.m)

tmp.df <- data.frame(sd_full = sd_full.v,
										 sd_ductal = sd_ductal.v,
										 type = droplevels(interaction(a, b)),
										 gene = genes.v)

tmp1.df <- tmp.df %>% dplyr::filter((abs(sd_full) > 0.9) | (abs(sd_ductal) > 0.9))


sd.p <- ggplot(tmp.df, aes(x = sd_full, y = sd_ductal, label = gene)) +
	geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_abline(slope = 1, intercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point(aes(color = type), size = 0.5,  alpha = 0.4,  shape = 1) +
	geom_text_repel(data = tmp1.df, size = 1.3, segment.size = 0.2, max.overlaps = 25) +
	scale_color_brewer(palette = "Set2", labels = c("Full genes", "Ductal genes", "Both genes"), name = "Type") + 
	labs( y = "SD of Ductal subset", x = "SD of full dataset", title = str_c("SD of top genes(n=", length(genes.v), ")")) 

mp <- plot_grid(pc1.p, pc2.p, sd.p + theme(legend.position = c(0, 1),
																					 legend.justification = c(0, 1)),
								nrow = 1, ncol = 3, label_size = 10, labels = "auto", align = "hv", axis = "tblr")

save_plot(here::here("figs", "response_figs", "extra.endoWeights2.pdf"), mp,
					base_height = 2, base_width = 2*1.5, nrow = 1, ncol = 3, device = cairo_pdf)







