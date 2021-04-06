rm(list=ls())
source(here::here("scripts/utils.R"))

endo.o <- qread(here::here("data/plotdata/endo.qs"))
endo_full.o <- qread(here::here("data/endo.qs"))


c1 <- reducedDim(endo.o, "tricycleEmbedding")[, 1]
c2 <- reducedDim(endo.o, "tricycleEmbedding")[, 2]
adjust.m <- t(apply(assay(endo_full.o, "log.s"), 1, function(x) {
	adj <- mean(x) + residuals(lm(x ~ c1 + c1))
}))

endo.o$top2a.adj <- adjust.m[which(rownames(endo_full.o) == "Top2a"), ]
plot(endo.o$tricyclePosition, endo.o$top2a.adj)
plot(endo.o$tricyclePosition, endo.o$top2a)



set.seed(1000)
adj.pca <- scater::calculatePCA(adjust.m,  ncomponents = 30)
reducedDim(endo.o, "adj.pca") <- adj.pca
endo.o <- scater::runUMAP(endo.o, dimred = 'adj.pca', exprs_values = "log.s", external_neighbors = FALSE,  pca = 30, min_dist = 0.5,  name = "adj.umap")

tmp.df <- data.frame(pc1.s = reducedDim(endo.o, "adj.umap")[, 1], pc2.s = reducedDim(endo.o, "adj.umap")[, 2],
										 cell_type = endo.o$cell_type)

umap.scat.p <- ggplot(tmp.df, aes(x = pc1.s, y = pc2.s, color = cell_type)) +
	geom_scattermore(pointsize = metadata(endo.o)$point.size, alpha = metadata(endo.o)$point.alpha) +
	scale_color_manual(values = metadata(endo.o)$clusters_colors, name = "Cell type", labels = levels(factor(endo.o$cell_type))) +
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
	labs( y = "UMAP-2", x = "UMAP-1", title = "UMAP (cell cycle adjusted)")


endo.top2a.p <- plotLoess(endo.o, "top2a", col.outname = "Top2A")
endo.top2a.adj.p <- plotLoess(endo.o, "top2a.adj", col.outname = "Top2A adj")




mp <- plot_grid(umap.scat.p,
								endo.top2a.p,
								endo.top2a.adj.p,
								nrow = 1, ncol = 3, label_size = 10, labels = "auto", align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.EndoAdj.pdf"), mp,
					base_height = 2, base_width = 2*1.5, nrow = 1, ncol = 3, device = cairo_pdf)




