rm(list=ls())
source(here::here("scripts/utils.R"))


all.df <- qread(here::here("data/plotdata/time.df.qs"))
### convert to mins
all.df[, 1:3] <- (all.df[, 1:3] * 10^(-9)) / 60
names(all.df)[1:3] <- c("seurat", "tricycle", "cyclone")


tmp.df <- all.df %>% gather(key = "method", value = "time", 1:3)
tmp.df$method <- factor(tmp.df$method, levels = c("cyclone", "seurat", "tricycle"))

time.p <- ggplot(tmp.df, aes(x = n, y = time, color = method)) +
	stat_summary(
		fun = median,
		geom = 'line',
		size = 0.5,
		alpha = 0.6,
		aes(group = method, colour = method)
	) +
	geom_jitter(size = 0.8, alpha = 0.6, width = 1000, height = 2, shape = 1, stroke = 0.3) +
	scale_color_brewer(palette = "Set2", name = "Method", labels = c("Cyclone", "Seurat", "TriCycle")) +
	guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
	labs( x = "Number of cells", y = "Time (mins)", title = "Benchmarking running time") +
	scale_x_continuous(name = "Number of cells", breaks = c(5000, 10000, 50000), labels = c("5000", "10000", "50000")) +
	theme(legend.position = c(0, 1),
				legend.justification = c(0, 1),
				legend.key = element_blank(),
				legend.key.size = unit(7, "pt"),
				axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5, angle = 30))


save_plot(here::here("figs", "sfigs", "sfig.runningtime.pdf"), time.p,
					base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 1, device = cairo_pdf)


