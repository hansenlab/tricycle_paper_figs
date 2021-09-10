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
	lambda_scaled <- sweep(lambda / colSums(lambda), MARGIN = 1, STATS = libsizes, FUN = "*")
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
simul_negbin2 <- function(locations, amplitudes, offset, bcv = 0.1, libsize = 2000) {
	## FIXME: no dependency on library size
	stopifnot(length(locations) == length(amplitudes))
	nGenes <- length(locations)
	lambda <- t(mapply(function(aa, dd) {
		aa * cos(2 * (timepoints - dd))
	}, amplitudes, locations)) + offset
	libsizes <- rep(libsize, nGenes)
	lambda_scaled <- sweep(lambda / colSums(lambda), MARGIN = 1, STATS = libsizes, FUN = "*")
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
# peak2.m <- assay(simul_negbin(rep(runif(2, min = 0, max = 2*pi), each = 50), rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")
# peak3.m <- assay(simul_negbin(rep(runif(3, min = 0, max = 2*pi), each = 34)[1:100], rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")
# peak5.m <- assay(simul_negbin(rep(runif(5, min = 0, max = 2*pi), each = 20), rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")
# peak10.m <- assay(simul_negbin(rep(runif(10, min = 0, max = 2*pi), each = 10), rep(3, 100), offset = 5, bcv = 0.1, libsize = 2000), "logcounts")
peak50.m <- rbind(assay(simul_negbin(runif(50, min = 0, max = 2*pi),  rep(3, 50), offset = 5, bcv = 0.1, libsize = 2000),"logcounts"), 
									assay(simul_negbin2(runif(50, min = 0, max = 2*pi),  rep(3, 50), offset = 5, bcv = 0.1, libsize = 2000),"logcounts"))


pca.o <- prcomp(t(peak50.m))
imp <- summary(pca.o)$importance

# gene fig
od <- sample(seq_along(timepoints))

tmp.df <- rbind(data.frame(timepoints = timepoints, pt= timepoints[od], exp = peak50.m[1,], g = "g1"),
								data.frame(timepoints = timepoints, pt= timepoints[od], exp = peak50.m[51,], g = "g2"))


a.p <- ggplot(tmp.df, aes(x = timepoints, y = exp, color = g)) +
	geom_scattermore(shape = 1, pointsize = 2.01, alpha = 0.7) +
	scale_color_manual(values = c("#F29360", "#86BBD8"), labels = c("Gene 1", "Gene2"), name = "Gene") +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi)  +
	labs(title = "Unobserved", x = "True timepoint \u03B8", y = "Expression") +
	theme(legend.position = c(0, 1),
				legend.justification = c(0, 1),
				legend.direction="horizontal",
				axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))


### pca fig
theta.hat <- as.numeric(circular::coord2rad(pca.o$x[, 1:2]))
theta.hat.shifted <- ifelse(theta.hat < theta.hat[1], 2*pi - theta.hat[1] + theta.hat, theta.hat - theta.hat[1])

tmp.df <- data.frame(pc1 = pca.o$x[,1], pc2 = pca.o$x[,2], 
										 pc3= pca.o$x[,3], pc4 = pca.o$x[,4],
										 timepoints = timepoints, theta = theta.hat.shifted)
b.p <- ggplot(tmp.df, aes(x = pc1, y = pc2)) +
	geom_scattermore(shape = 1, pointsize = 2.01, alpha = 0.7) +
	labs(title = "PC1 and PC2", x = sprintf("PC1 (%2.1f%%)", imp[2,1]*100), y = sprintf("PC2 (%2.1f%%)", imp[2,2]*100)) +
	annotate("point", x = 0, 
					 y = 0, shape = 16, size = 1, alpha = 0.8) + 
	annotate("segment", x = 0, xend =  4 ,
					 y = 0, yend = 0, linetype = "dashed", alpha = 0.8) + 
	annotate("segment", x = 0, xend =  cos(pi/6) * 3,
					 y = 0, yend =  sin(pi/6) *3,
					 linetype = "dashed", alpha = 0.8) +
	annotate("text",  alpha = 0.8, x = 4, y = 1, label = expression(hat(theta)), parse = FALSE, size = 3) +
	annotate("point", x = pca.o$x[1,1], y = pca.o$x[1,2], color = "red", shape = 16, size = 1.2) +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))

c.p <- ggplot(tmp.df, aes(x = pc3, y = pc4)) +
	geom_scattermore(shape = 1, pointsize = 2.01, alpha = 0.7) +
	labs(title = "PC3 and PC4", x = sprintf("PC3 (%2.1f%%)", imp[2,4]*100), y = sprintf("PC4 (%2.1f%%)", imp[2,4]*100)) +
	# annotate("point", x = 0, 
	# 				 y = 0, shape = 16, size = 1, alpha = 0.8) + 
	# annotate("segment", x = 0, xend =  4 ,
	# 				 y = 0, yend = 0, linetype = "dashed", alpha = 0.8) + 
	# annotate("segment", x = 0, xend =  cos(pi/6) * 3,
	# 				 y = 0, yend =  sin(pi/6) *3,
	# 				 linetype = "dashed", alpha = 0.8) +
	# annotate("text",  alpha = 0.8, x = 3, y = 0.6, label = expression(hat(theta)), parse = FALSE, size = 3) +
	# annotate("point", x = pca$x[1,3], y = pca$x[1,4], color = "red", shape = 16, size = 1.2) +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))

d.p <- ggplot(tmp.df, aes(x = timepoints, y = theta)) +
	geom_scattermore(shape = 1, pointsize = 2.01, alpha = 0.7) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi)  +
	scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi)  +
	labs(title = "Recovered timepoints using PC1-2", x = "True timepoint \u03B8", y = "Recovered timepoint \u03B8") +
	theme(legend.position = c(0, 1),
				legend.justification = c(0, 1),
				legend.direction="horizontal",
				axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))


# tmp.df <- rbind(data.frame(theta = theta.hat.shifted, exp = peak50.m[1,], g = "g1"),
# 								data.frame(theta = theta.hat.shifted, exp = peak50.m[50,], g = "g2"))
# 
# d.p <- ggplot(tmp.df, aes(x = theta, y = exp, color = g)) +
# 	geom_scattermore(shape = 1, pointsize = 2.01, alpha = 0.7) +
# 	scale_color_manual(values = c("#F29360", "#86BBD8"), labels = c("Gene 1", "Gene2"), name = "Gene") +
# 	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi)  +
# 	labs(title = "Recovered", x = expression(paste("Recoverd timepoint ", hat(theta))), y = "Expression") +
# 	theme(legend.position = c(0, 1),
# 				legend.justification = c(0, 1),
# 				legend.direction="horizontal",
# 				axis.title.y = element_text(size = 8),
# 				axis.title.x = element_text(size = 8))



mp <- plot_grid(a.p, b.p, c.p, d.p, 
nrow = 2, ncol = 2, label_size = 10, labels = "auto", align = "v", axis = "lr")


save_plot(here::here("figs", "response_figs", "response.2peaksloc.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 2, device = cairo_pdf)




