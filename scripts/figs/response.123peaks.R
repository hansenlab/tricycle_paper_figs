rm(list=ls())
source(here::here("scripts/utils.R"))


library(circular)
library(SingleCellExperiment)
library(scuttle)



n_time <- 5000
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
peak1.m <- assay(simul_negbin(rep(runif(1, min = 0, max = 2*pi), each = 100), rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")
peak2.m <- assay(simul_negbin(rep(runif(2, min = 0, max = 2*pi), each = 50), rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")
peak3.m <- assay(simul_negbin(rep(runif(3, min = 0, max = 2*pi), each = 34)[1:100], rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")


getPCA <- function(xx) {
	pca <- prcomp(t(xx))
	imp <- summary(pca)$importance
	return(pca)
}


plotPCA <- function(pca, npeak = "2", point.size = 2.01, gene_xlims, gene_ylims, cell_xlims, cell_ylims) {
	
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
	
	mp <- plot_grid(p1, p2, p3, p4,
									nrow = 1, ncol = 4, label_size = 10, labels = NULL, align = "v", axis = "lr")
	return(mp)
}

pca.lo <- list(getPCA(peak1.m),
							 getPCA(peak2.m),
							 getPCA(peak3.m))

tmp1.m <- do.call(rbind, lapply(pca.lo, function(x) x$rotation[,1:2]))
gene_xlims <- xlim(range(tmp1.m[, 1]))
gene_ylims <- ylim(range(tmp1.m[, 2]))

tmp2.m <- do.call(rbind, lapply(pca.lo, function(x) x$x[,1:2]))
cell_xlims <- xlim(range(tmp2.m[, 1]))
cell_ylims <- ylim(range(tmp2.m[, 2]))

mp <- plot_grid(plotlist = lapply(seq_along(c(1, 2, 3)), function(i) {
	plotPCA(pca = pca.lo[[i]], npeak = c(1, 2, 3)[i], gene_xlims = gene_xlims, gene_ylims =  gene_ylims , cell_xlims = cell_xlims, cell_ylims = cell_ylims)
}),
nrow = 3, ncol = 1, label_size = 10, labels = "auto", align = "v", axis = "lr")


save_plot(here::here("figs", "response_figs", "response.123peaks.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 3, ncol = 4, device = cairo_pdf)


