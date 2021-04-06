rm(list=ls())
source(here::here("scripts/utils.R"))

library(randomcoloR)
colData.df <- qread(here::here("data", "plotdata", str_c("fetal_colData.qs"))) %>% as.data.frame()

text.df <- data.frame(tissue = levels(factor(colData.df$Organ)),
											x = tapply(colData.df$Global_umap_1, colData.df$Organ, mean),
											y = tapply(colData.df$Global_umap_2, colData.df$Organ, mean))
text.df[text.df$tissue == "Adrenal", 2:3] <- c(13, 15)
text.df[text.df$tissue == "Cerebellum", 2:3] <- c(-10, -8.5)
text.df[text.df$tissue == "Cerebrum", 2:3] <- c(-10, -7.5)
text.df[text.df$tissue == "Eye", 2:3] <- c(14.5, -13)
text.df[text.df$tissue == "Heart", 2:3] <- c(10, -17)
text.df[text.df$tissue == "Intestine", 2:3] <- c(10, 3)
text.df[text.df$tissue == "Kidney", 2:3] <- c(10, -1)
text.df[text.df$tissue == "Liver", 2:3] <- c(18, -1)
text.df[text.df$tissue == "Lung", 2:3] <- c(5, 5)
text.df[text.df$tissue == "Muscle", 2:3] <- c(5, -12)
text.df[text.df$tissue == "Pancreas", 2:3] <- c(9, 2)
text.df[text.df$tissue == "Placenta", 2:3] <- c(7, -3)
text.df[text.df$tissue == "Spleen", 2:3] <- c(12, -5)
text.df[text.df$tissue == "Stomach", 2:3] <- c(3, 3)
text.df[text.df$tissue == "Thymus", 2:3] <- c(13, 0.5)

theta.p <- 	ggplot(colData.df, aes(x = Global_umap_1, y = Global_umap_2, color = tricyclePosition)) +
	geom_scattermore(pointsize = 0.3,  alpha = 0.6) +
	scale_color_gradientn(name = "", limits = range(0, 2 * pi), 
												breaks = seq(from = 0, to = 2 * pi, length.out = 500) ,
												colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"), 
												guide = FALSE) + 
	labs( x = "UMAP-1", y = "UMAP-2", title = str_c( "Human fetal tissues atlas (n=", nrow(colData.df), ")")) +
	theme(plot.margin = unit(c(5, 5, 5, 6), "pt"))  
	#geom_text(data = text.df, aes(x = x, y= y, label = tissue), inherit.aes = FALSE, size = 2, hjust = 0.5, vjust = 0.5)




tmp.df <- colData.df

tmp.df$cycling <- (tmp.df$tricyclePosition >= 0.25 * pi) & (tmp.df$tricyclePosition <= 1.5 * pi)

tmp.df %>% group_by(Organ, cycling) %>% summarise(n = n()) %>%
	mutate(freq = n / sum(n))

frac.df <- tmp.df %>% group_by(Organ, cycling) %>% summarise(n = n()) %>%
	mutate(freq = n / sum(n))

ordered.df <- frac.df %>% dplyr::filter(cycling) %>% arrange(desc(freq))
order.Organ.v <- as.character(ordered.df$Organ)


frac.df$Organ <- factor(frac.df$Organ, levels = order.Organ.v)


pct.p <- ggplot(frac.df, aes(x = Organ, y = freq, color = cycling)) +
	geom_point() +
	scale_color_brewer(palette = "Set1", name = "Actively proliferating") +
	geom_hline(yintercept = 0.5, size = 0.3, linetype = "dashed", alpha = 0.7) +
	scale_x_discrete(labels = str_c(levels(frac.df$Organ), "\n(n=", table(tmp.df$Organ), ")"), name = "") +
	labs(title = "Percentage of actively proliferating cells in human fetal tissue atlas", y = "Percentage") + 
	ylim(c(0, 1)) +
	theme(axis.text.x = element_text(size = 5,  vjust = 0.5, hjust = 0.5, angle = 30),
				legend.position = c(1, 1),
				legend.justification = c(1, 1),
				plot.margin = unit(c(5, 5, 6, 6), "pt"),
				legend.direction="horizontal")





mp <- plot_grid(theta.p, pct.p,
								nrow = 2, ncol = 1, rel_heights = c(1, 0.5), label_size = 10, labels = "auto", align = "v", axis = "lr")

mp <- ggdraw(mp) +
	draw_plot(circle_scale_legend(text.size = 1.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "pt")), 0.99, 0.95, .2, .15, hjust = 1, vjust = 1)

save_plot(here::here("figs", "main", "main.fetalatlas.pdf"), mp,
					base_height = 2, base_width = 2*1.1, nrow = 2.5, ncol = 2, device = cairo_pdf)


