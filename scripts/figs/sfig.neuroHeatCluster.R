rm(list=ls())
source(here::here("scripts/utils.R"))

library(pheatmap)

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
neurosphere_full.o <- qread(here::here("data/neurosphere.qs"))


neurosphere.o$pca.theta <- as.numeric(coord2rad(reducedDim(neurosphere.o, "go.pca")[, 1:2]))
## heatmap
## neurosphere
data.m <- assay(neurosphere_full.o, "log.s")[rownames(metadata(neurosphere.o)$rotation), ]
data.m %<>% t() %>% scale(center = TRUE, scale = TRUE) %>% t()

loess.ll <- bplapply(seq_len(nrow(data.m)), function(j) {
	v <- data.m[j, ]
	loess.l <- fit_periodic_loess(neurosphere.o$pca.theta, y = v, length.out = 1000)
	return(loess.l)
}, BPPARAM = MulticoreParam(workers = 10L))
qsave(loess.ll, file = here::here("data/plotdata/neuro.loess.qs"))

rank.v <- unlist(lapply(loess.ll, function(x) which.max(x$pred.df$y)))

order.idx <- order(neurosphere.o$pca.theta)
pheno.df <- data.frame(SchwabeCC = neurosphere.o$CCStage[order.idx],
											 "PCA\u03B8.\u03C0" = as.character(c(0.5, 1, 1.5, 2))[as.numeric(cut(neurosphere.o$pca.theta[order.idx], breaks = seq(0, 2 * pi, length.out = 5), include.lowest = TRUE))])
rownames(pheno.df) <- colnames(data.m)[order.idx]
heat.p <- pheatmap(data.m[order(rank.v), order.idx], cluster_rows = FALSE, cluster_cols = FALSE, #breaks = seq(from = -.85, to = .85, length.out = 101),
				 main = str_c("mNeurosphere (n.genes = 500; n.cells=", ncol(neurosphere.o), ")"), show_rownames = FALSE, show_colnames = FALSE, fontsize = 8,
				 clustering_method = "ward.D2", #clustering_distance_rows = dist.o,  clustering_distance_cols = dist.o,
				 annotation_col = pheno.df, annotation_colors = list(SchwabeCC = setNames(ccColors.v, ccLabels.v),
				 																										"PCA\u03B8.\u03C0" = setNames(brewer.pal(9, "Blues")[c(2, 4, 6, 8)], as.character(c(0.5, 1, 1.5, 2)))),
				 breaks = seq(-2, 2, length.out = 100),
				 silent = TRUE,
				 width = 9, height = 6)




## neuroshpere
SelGenes.v <- rownames(metadata(neurosphere.o)$rotation)
tmp.l <- lapply(seq_along(loess.ll), function(i) loess.ll[[i]]$pred.df %>% add_column(idx = i))
	
### 3 clsuters
range.v <- sapply(tmp.l, function(x) diff(range(x$y)))
max.theta.v <-  sapply(tmp.l, function(x) x$x[which.max(x$y)])
cluster1.idx <- which(range.v < 0.5 )
cluster2.idx <- which((range.v >= 0.5 ) & (max.theta.v < 1 * pi))
cluster3.idx <- which((range.v >= 0.5) & (max.theta.v >= 1 * pi))


tmp.df <- do.call(rbind, tmp.l)
ylim.v <- range(tmp.df$y)
neurosphere.all.z.p <- ggplot(tmp.df, aes(x = x, y = y, group = idx)) +
	geom_path(linetype = "solid", color = "black", size = 0.3, alpha = 0.3, show.legend = FALSE) +
	# geom_rug(sides="bl", alpha = 0.5, size = 0.05) +
	#guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
	labs( y = "z-scores", x = "PCA \u03B8", title = str_c("All genes (N=", length(max.theta.v), ")")) +
	ylim(ylim.v) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2) * pi) +
	theme(legend.position = 'none')


tmp.df <- do.call(rbind, tmp.l[cluster1.idx])
neurosphere.cluster1.z.p <- ggplot(tmp.df, aes(x = x, y = y, group =  idx)) +
	geom_path(linetype = "solid", color = "black", size = 0.3, alpha = 0.3, show.legend = FALSE) +
	# geom_rug(sides="bl", alpha = 0.5, size = 0.05) +
	#guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
	labs( y = "z-scores", x = "PCA \u03B8", title = str_c("Cluster1 (N=", length(cluster1.idx), ")")) +
	ylim(ylim.v) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2) * pi) +
	theme(legend.position = 'none')

tmp.df <- do.call(rbind, tmp.l[cluster2.idx])
neurosphere.cluster2.z.p <- ggplot(tmp.df, aes(x = x, y = y, group =  idx)) +
	geom_path(linetype = "solid", color = "black", size = 0.3, alpha = 0.3, show.legend = FALSE) +
	# geom_rug(sides="bl", alpha = 0.5, size = 0.05) +
	#guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
	labs( y = "z-scores", x = "PCA \u03B8", title = str_c("Cluster2 (N=", length(cluster2.idx), ")")) +
	ylim(ylim.v) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2) * pi) +
	theme(legend.position = 'none')

tmp.df <- do.call(rbind, tmp.l[cluster3.idx])
neurosphere.cluster3.z.p <- ggplot(tmp.df, aes(x = x, y = y, group =  idx)) +
	geom_path(linetype = "solid", color = "black", size = 0.3, alpha = 0.3, show.legend = FALSE) +
	# geom_rug(sides="bl", alpha = 0.5, size = 0.05) +
	#guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
	labs( y = "z-scores", x = "PCA \u03B8", title = str_c("Cluster3 (N=", length(cluster3.idx), ")")) +
	ylim(ylim.v) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2) * pi) +
	theme(legend.position = 'none')


mp2 <- plot_grid(neurosphere.all.z.p, 
								neurosphere.cluster1.z.p,
								neurosphere.cluster2.z.p,
								neurosphere.cluster3.z.p,
								nrow = 2, ncol = 2, labels = c("b", "c", "d", "e"), align = "hv", axis = "lrtb")
mp <- plot_grid(heat.p$gtable, mp2,
								nrow = 1, ncol = 2, labels = c("a", " "), align = "none", axis = "none")

save_plot(here::here("figs", "sfigs", "sfig.neuroHeatCluster.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 4, device = cairo_pdf)


### something went wrong for pheatmap jpg output - no background. 
# save_plot(here::here("figs", "sfigs", "sfig.neuroHeatCluster.jpg"), mp,
# 					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 4, type = "cairo")




