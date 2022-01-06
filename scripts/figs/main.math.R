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
xx <- simul( c( rep(0.2, 50), rep(1.2, 50) ) , c(rep(0.5, 50), rep(1, 50)))

### figure math
od <- sample(seq_along(timepoints))

tmp.df <- rbind(data.frame(timepoints = timepoints, pt= timepoints[od], exp = -xx[51,], g = "g1"),
								data.frame(timepoints = timepoints, pt= timepoints[od], exp = -xx[1,], g = "g2"))


a.p <- ggplot(tmp.df, aes(x = timepoints, y = exp, color = g)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	scale_color_manual(values = c("#F29360", "#86BBD8"), labels = c("Gene 1", "Gene2"), name = "Gene") +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi)  +
	labs(title = "Unobserved", x = "Timepoint \u03B8", y = "Expression") +
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				# legend.direction="horizontal",
				axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8),
				legend.title = element_text(size = 6),
				legend.text = element_text(size = 6))

b.p <- ggplot(tmp.df, aes(x = pt, y = exp, color = g)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	scale_color_manual(values = c("#F29360", "#86BBD8"), labels = c("Gene 1", "Gene2"), name = "Gene") +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi)  +
	labs(title = "Observed", x = "Permuted timepoint \u03B8\'", y = "Expression") +
	theme(legend.position = "none",
				legend.justification = c(1, 1),
				legend.direction="horizontal",
				axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))

### pca
pca <- prcomp(t(xx))
imp <- summary(pca)$importance


tmp.df <- data.frame(x = pca$x[,1], y = pca$x[,2])
c.p <- ggplot(tmp.df, aes(x =x, y = y)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	labs(title = "Inference of the time order", x = sprintf("PC1 (%2.1f%%)", imp[2,1]*100), y = sprintf("PC2 (%2.1f%%)", imp[2,2]*100)) +
	annotate("point", x = 0, 
					 y = 0, shape = 16, size = 1, alpha = 0.8) + 
	annotate("segment", x = 0, xend =  4 ,
					 y = 0, yend = 0, linetype = "dashed", alpha = 0.8) + 
	annotate("segment", x = 0, xend =  cos(pi/6) * 3,
					 y = 0, yend =  sin(pi/6) *3,
					 linetype = "dashed", alpha = 0.8) +
	annotate("text",  alpha = 0.8, x = 3, y = 0.6, label = expression(hat(theta)), parse = FALSE, size = 3) +
	annotate("point", x = pca$x[1,1], y = pca$x[1,2], color = "red", shape = 16, size = 1.2) +
	theme(axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))


theta.hat <- as.numeric(circular::coord2rad(pca$x[, 1:2] * rep(c(1, 7.5/3), each = nrow(pca$x))))
theta.hat.shifted <- ifelse(theta.hat < theta.hat[1], 2*pi - theta.hat[1] + theta.hat, theta.hat - theta.hat[1])

tmp.df <- rbind(data.frame(theta = theta.hat.shifted, exp = -xx[51,], g = "g1"),
								data.frame(theta = theta.hat.shifted, exp = -xx[1,], g = "g2"))

d.p <- ggplot(tmp.df, aes(x = theta, y = exp, color = g)) +
	geom_scattermore(shape = 1, pointsize = 4.01, alpha = 0.7) +
	scale_color_manual(values = c("#F29360", "#86BBD8"), labels = c("Gene 1", "Gene2"), name = "Gene") +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi)  +
	labs(title = "Recovered", x = expression(paste("Recoverd timepoint ", hat(theta))), y = "Expression") +
	theme(legend.position = "none",
				legend.justification = c(0, 1),
				legend.direction="horizontal",
				axis.title.y = element_text(size = 8),
				axis.title.x = element_text(size = 8))



mp1 <- plot_grid(a.p, b.p,
								nrow = 2, ncol = 1, label_size = 10, labels = "auto", align = "v", axis = "lr")

mp2 <- plot_grid(ggplot() + theme_nothing(),
								 c.p,
								 ggplot() + theme_nothing(),
								 nrow = 3, ncol = 1, 
								 rel_heights = c(0.5, 1, 0.5),
								 label_size = 10, labels = c("", "c", ""), align = "v", axis = "lr")
mp3 <- plot_grid(ggplot() + theme_nothing(),
								 d.p,
								 ggplot() + theme_nothing(),
								 nrow = 3, ncol = 1, 
								 rel_heights = c(0.5, 1, 0.5),
								 label_size = 10, labels = c("", "d", ""), align = "v", axis = "lr")

mp <- plot_grid(mp1, ggplot() + theme_nothing(), mp2, ggplot() + theme_nothing(), mp3,
								nrow = 1, ncol = 5,  rel_widths = c(1, .3, 1, 0.3, 1),
								label_size = 10, labels = NULL, align = "v", axis = "lr")

mp <- ggdraw(mp) +
	draw_plot(ggplot() + geom_segment(data = data.frame(x = 0, y = 25, xend = 3, yend = 30),
												 aes(x = x, y = y, xend = xend, yend = yend),
												 arrow = arrow(length = unit(0.3, "cm"), type="closed" )) +
							theme_nothing(), 0.28, 0.25, .08, .15, hjust = 0, vjust = 0) +
	draw_plot(ggplot() + geom_segment(data = data.frame(x = 0, y = 30, xend = 3, yend = 25),
																		aes(x = x, y = y, xend = xend, yend = yend),
																		arrow = arrow(length = unit(0.3, "cm"), type="closed")) +
							theme_nothing(), 0.28, 0.65, .08, .15, hjust = 0, vjust = 0) +
	draw_plot(ggplot() + geom_segment(data = data.frame(x = 0, y = 30, xend = 3, yend = 30),
																		aes(x = x, y = y, xend = xend, yend = yend),
																		arrow = arrow(length = unit(0.3, "cm"), type="closed")) +
							theme_nothing(), 0.64, 0.44, .08, .15, hjust = 0, vjust = 0)

save_plot(here::here("figs", "main", "main.math.pdf"), mp,
					base_height = 2/1.1, base_width = 2*1.1 / 1.1, nrow = 2, ncol = 3.2, device = cairo_pdf)

