rm(list=ls())
source(here::here("scripts/utils.R"))


### plot the total UMIs
cell.df <- do.call(rbind, lapply(c("Intestine", "Kidney", "Pancreas", "Stomach"), function(o) {
	sce.o <- qread(here::here("data/fetal", str_c(o, ".qs")))
	return(data.frame(tissue = o, totalUMIs = sce.o$TotalUMIs))
}))

nuclei.df <- do.call(rbind, lapply(c("Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart", "Liver", "Lung", "Muscle", "Placenta", "Spleen", "Thymus"), function(o) {
	sce.o <- qread(here::here("data/fetal/nuclei", str_c(o, ".qs"))) 
	return(data.frame(tissue = o, totalUMIs = sce.o$TotalUMIs))
}))


tmp.df <- rbind(cell.df, nuclei.df)
tmp.df$tissue <- factor(tmp.df$tissue, levels = c("Intestine", "Kidney", "Pancreas", "Stomach", "Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart", "Liver", "Lung", "Muscle", "Placenta", "Spleen", "Thymus"))

totalumi.box.p <- ggplot(tmp.df, aes(x = tissue, y = log10(totalUMIs))) + 
	geom_violin(scale = "width", lwd = 0.3, width = 0.8) +
	geom_boxplot(width = .1,  outlier.shape = NA, lwd = 0.3) +
	scale_x_discrete(labels = str_c(levels(tmp.df$tissue), "\n(n=", table(tmp.df$tissue), ")"), name = "") +
	scale_y_continuous(limits = c(log10(200),  max(log10(tmp.df$totalUMIs))), name = "Number of totalUMIs (log10 scaled)", breaks = log10(c(200, 300, 500, 1000, 5000, 10000, 50000, 100000, 200000)),
										 labels = c(200, 300, 500, 1000, 5000, 10000, 50000, 100000, 200000)) +
	labs(title = "Number of totalUMIs of human fetal tissues atlas") +
	geom_vline(xintercept = 4.5, size = 0.3, linetype = "dashed", alpha = 0.7) +
	theme(axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5, angle = 30))

save_plot(here::here("figs", "sfigs", "sfig.fetaltotalumi.pdf"), totalumi.box.p,
					base_height = 2, base_width = 2*1.2, nrow = 1.2, ncol = 2.5, device = cairo_pdf)


tapply(tmp.df$totalUMIs, tmp.df$tissue, median)

#Adrenal Cerebellum   Cerebrum        Eye      Heart  Intestine     Kidney      Liver       Lung     Muscle   Pancreas   Placenta     Spleen    Stomach     Thymus 
#489        796        782        408        409        503        759        452        431        413        892        368        372        429        354

