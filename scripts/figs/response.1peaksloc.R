rm(list=ls())
source(here::here("scripts/utils.R"))


library(circular)

n_time <- 1000
timepoints <- seq(from = 0, to = 2*pi, length.out = n_time)
simul <- function(locations, amplitudes) {
	t(mapply(function(aa, dd) {
		y <- aa * cos(timepoints - dd)
		y <- y + rnorm(length(timepoints), mean = 0, sd = 0.2)
		return(y)
	}, amplitudes, locations))
}
set.seed(2387546)
xx <- simul( c( rep(0.2, 50), rep(1.2, 50) ), c(rep(0.5, 50), rep(1, 50)))

### figure math
od <- sample(seq_along(timepoints))

tmp.df <- rbind(data.frame(timepoints = timepoints, pt= timepoints[od], exp = -xx[51,], g = "g1"),
								data.frame(timepoints = timepoints, pt= timepoints[od], exp = -xx[1,], g = "g2"))


a.p <- ggplot(tmp.df, aes(x = timepoints, y = exp, color = g)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	scale_color_manual(values = c("#F29360", "#86BBD8"), labels = c("Gene 1", "Gene2"), name = "Gene") +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi)  +
	labs(title = "Two peak locations across genes", x = "Timepoint \u03B8", y = "Expression") +
	theme(legend.position = c(0, 1),
				legend.justification = c(0, 1),
				legend.direction="horizontal",
				axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))

### pca
pca <- prcomp(t(xx))
imp <- summary(pca)$importance


tmp.df <- data.frame(x = pca$x[,1], y = pca$x[,2], t = timepoints)
b.p <- ggplot(tmp.df, aes(x =x, y = y)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	labs(title = "Two peak locations across genes", x = sprintf("PC1 (%2.1f%%)", imp[2,1]*100), y = sprintf("PC2 (%2.1f%%)", imp[2,2]*100)) +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))

c.p <- ggplot(tmp.df, aes(x = t, y = x)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	labs(title = "Two peak locations across genes", y = sprintf("PC1 (%2.1f%%)", imp[2,1]*100), x = "Timepoint \u03B8") +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))

d.p <- ggplot(tmp.df, aes(x = t, y = y)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	labs(title = "Two peak locations across genes", y = sprintf("PC2 (%2.1f%%)", imp[2,2]*100), x = "Timepoint \u03B8") +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))

### simulate all genes peaking at the same location
xx <- simul( rep(1.2, 100), c(rep(0.5, 50), rep(1, 50)))

### figure math
od <- sample(seq_along(timepoints))

tmp.df <- rbind(data.frame(timepoints = timepoints, pt= timepoints[od], exp = -xx[51,], g = "g1"),
								data.frame(timepoints = timepoints, pt= timepoints[od], exp = -xx[1,], g = "g2"))


e.p <- ggplot(tmp.df, aes(x = timepoints, y = exp, color = g)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	scale_color_manual(values = c("#F29360", "#86BBD8"), labels = c("Gene 1", "Gene2"), name = "Gene") +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi)  +
	labs(title = "One peak location across genes", x = "Timepoint \u03B8", y = "Expression") +
	theme(legend.position = c(0, 1),
				legend.justification = c(0, 1),
				legend.direction="horizontal",
				axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))

### pca
pca <- prcomp(t(xx))
imp <- summary(pca)$importance


tmp.df <- data.frame(x = pca$x[,1], y = pca$x[,2], t = timepoints)
f.p <- ggplot(tmp.df, aes(x =x, y = y)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	labs(title = "One peak location across genes", x = sprintf("PC1 (%2.1f%%)", imp[2,1]*100), y = sprintf("PC2 (%2.1f%%)", imp[2,2]*100)) +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))

g.p <- ggplot(tmp.df, aes(x = t, y = x)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	labs(title = "One peak location across genes", y = sprintf("PC1 (%2.1f%%)", imp[2,1]*100), x = "Timepoint \u03B8") +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))

h.p <- ggplot(tmp.df, aes(x = t, y = y)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	labs(title = "One peak location across genes", y = sprintf("PC2 (%2.1f%%)", imp[2,2]*100), x = "Timepoint \u03B8") +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))



mp <- plot_grid(a.p, b.p, c.p, d.p, e.p, f.p, g.p, h.p,
								 nrow = 2, ncol = 4, label_size = 10, labels = "auto", align = "v", axis = "lr")

save_plot(here::here("figs", "response_figs", "response.1peakloc.pdf"), mp,
					base_height = 2/1.1, base_width = 2*1.1 / 1.1, nrow = 2, ncol = 4, device = cairo_pdf)

