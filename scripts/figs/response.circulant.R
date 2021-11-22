rm(list=ls())
source(here::here("scripts/utils.R"))


library(circular)
library(SingleCellExperiment)
library(scuttle)



n_time <- 100
timepoints <- seq(from = 0, to = 2*pi, length.out = n_time)
simul_negbin <- function(locations, amplitudes, offset, bcv = 0.1, libsize = 2000) {
	## FIXME: no dependency on library size
	stopifnot(length(locations) == length(amplitudes))
	nGenes <- length(locations)
	lambda <- t(mapply(function(aa, dd) {
		aa * cos(timepoints - dd)
	}, amplitudes, locations)) + offset
	libsizes <- rep(libsize, nGenes)
	lambda_scaled <- sweep(t(t(lambda) / colSums(lambda)), MARGIN = 1, STATS = libsizes, FUN = "*")
	# lambda_scaled <- lambda
	bcv <- rep(bcv, nGenes)
	lambda_bcv <- rgamma(nGenes * length(timepoints),
											 shape = 1/bcv^2,
											 scale = lambda_scaled * bcv^2)
	counts <- matrix(rpois(nGenes * length(timepoints), lambda = lambda_bcv), nrow = nGenes)
	sce <- SingleCellExperiment(assays = list(counts = counts))
	sce <- logNormCounts(sce)
	sce
}

set.seed(1514)
peak2.m <- assay(simul_negbin(rep(runif(2, min = 0, max = 2*pi), each = 50), rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")
peak3.m <- assay(simul_negbin(rep(runif(3, min = 0, max = 2*pi), each = 34)[1:100], rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")
peak5.m <- assay(simul_negbin(rep(runif(5, min = 0, max = 2*pi), each = 20), rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")
peak10.m <- assay(simul_negbin(rep(runif(10, min = 0, max = 2*pi), each = 10), rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")
peak50.m <- assay(simul_negbin(rep(runif(50, min = 0, max = 2*pi), each = 2), rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")
peak100.m <- assay(simul_negbin(rep(runif(100, min = 0, max = 2*pi), each = 1), rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")



getPCA <- function(xx) {
	pca <- prcomp(t(xx))
	attr(pca, "cov") <- cov(xx)
	return(pca)
}


plotPCA <- function(pca, npeak = "2", point.size = 4.01, gene_xlims, gene_ylims, cell_xlims, cell_ylims) {
	
	imp <- summary(pca)$importance
	
	tmp.df <- data.frame(x = pca$rotation[,1], y = pca$rotation[,2])
	
	p1 <- ggplot(tmp.df, aes(x =x, y = y)) +
		geom_scattermore(shape = 1, pointsize = point.size + 3, alpha = 0.7) +
		labs(title = str_c("Gene PCA (nPeak=", npeak, ")"), x = "PC1", y = "PC2") +
		gene_xlims + gene_ylims
	
	tmp.df <- data.frame(x = pca$x[,1], y = pca$x[,2], timepoints = timepoints)
	p2 <- ggplot(tmp.df, aes(x =x, y = y)) +
		geom_scattermore(shape = 1, pointsize = point.size, alpha = 0.7) +
		labs(title = str_c("Cell PCA (nPeak=", npeak, ")"), x = sprintf("PC1 (%2.1f%%)", imp[2,1]*100), y = sprintf("PC2 (%2.1f%%)", imp[2,2]*100)) +
		annotate("point", x = pca$x[1,1], y = pca$x[1,2], color = "red", shape = 16, size = 1.5) +
		cell_xlims + cell_ylims
	
	p3 <- ggplot(tmp.df, aes(x = timepoints, y = x)) +
		geom_scattermore(shape = 1, pointsize = point.size, alpha = 0.7) +
		labs(title = str_c("Cell PC1 (nPeak=", npeak, ")"), x = "Timepoint \u03B8", y = "PC1")
	
	p4 <- ggplot(tmp.df, aes(x = timepoints, y = y)) +
		geom_scattermore(shape = 1, pointsize = point.size, alpha = 0.7) +
		labs(title = str_c("Cell PC2 (nPeak=", npeak, ")"), x = "Timepoint \u03B8", y = "PC2")
	
	p5 <- pheatmap::pheatmap(attr(pca, "cov"), cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
													 main = str_c("Sample covariance matrix (nPeak=", npeak, ")"), fontsize = 12)[[4]]
		
	mp <- plot_grid(p2 + theme(
		axis.text.y = element_text(size = 10),
		axis.text.x = element_text(size = 10),
		axis.title.y = element_text(size = 11),
		axis.title.x = element_text(size = 11),
		plot.title = element_text(face = "bold", size = 14, hjust = 0.5)),
		p5, nrow = 1, ncol = 2, label_size = 10, labels = NULL, align = "v", axis = "lr")
	return(mp)
}

pca.lo <- list(getPCA(peak2.m),
							 getPCA(peak3.m),
							 getPCA(peak5.m),
							 getPCA(peak10.m),
							 getPCA(peak50.m),
							 getPCA(peak100.m))

tmp1.m <- do.call(rbind, lapply(pca.lo, function(x) x$rotation[,1:2]))
gene_xlims <- xlim(range(tmp1.m[, 1]))
gene_ylims <- ylim(range(tmp1.m[, 2]))

tmp2.m <- do.call(rbind, lapply(pca.lo, function(x) x$x[,1:2]))
cell_xlims <- xlim(range(tmp2.m[, 1]))
cell_ylims <- ylim(range(tmp2.m[, 2]))

mp <- plot_grid(plotlist = lapply(seq_along(c(2, 3, 5, 10, 50, 100)), function(i) {
	plotPCA(pca = pca.lo[[i]], npeak = c(2, 3, 5, 10, 50, 100)[i], gene_xlims = gene_xlims, gene_ylims =  gene_ylims , cell_xlims = cell_xlims, cell_ylims = cell_ylims)
}),
nrow = 3, ncol = 2, label_size = 17, labels = "auto", align = "v", axis = "lr")


save_plot(here::here("figs", "response_figs", "response.circulant.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 3*2, ncol = 4*2, device = cairo_pdf)


