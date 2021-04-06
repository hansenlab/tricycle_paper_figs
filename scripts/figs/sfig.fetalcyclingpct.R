rm(list=ls())
source(here::here("scripts/utils.R"))


### plot the total UMIs
cell.df <- do.call(rbind, lapply(c("Intestine", "Kidney", "Pancreas", "Stomach"), function(o) {
	sce.o <- qread(here::here("data/fetal", str_c(o, ".qs")))
	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
	return(data.frame(tissue = o, totalUMIs = sce.o$TotalUMIs, tricyclePosition = sce.o$tricyclePosition, cell_type = sce.o$Main_cluster_name, day = sce.o$Development_day))
}))

nuclei.df <- do.call(rbind, lapply(c("Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart", "Liver", "Lung", "Muscle", "Placenta", "Spleen", "Thymus"), function(o) {
	sce.o <- qread(here::here("data/fetal/nuclei", str_c(o, ".qs"))) 
	point.size <- ifelse(ncol(sce.o) > 1000000, 1.01, 2.01)
	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
	
	return(data.frame(tissue = o, totalUMIs = sce.o$TotalUMIs, tricyclePosition = sce.o$tricyclePosition, cell_type = sce.o$Main_cluster_name, day = sce.o$Development_day))
}))


tmp.df <- rbind(cell.df, nuclei.df)
tmp.df$tissue <- factor(tmp.df$tissue, levels = c("Intestine", "Kidney", "Pancreas", "Stomach", "Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart", "Liver", "Lung", "Muscle", "Placenta", "Spleen", "Thymus"))

tmp.df$cycling <- (tmp.df$tricyclePosition >= 0.25 * pi) & (tmp.df$tricyclePosition <= 1.5 * pi)

tmp.df %>% group_by(tissue, cycling) %>% summarise(n = n()) %>%
	mutate(freq = n / sum(n))

frac.df <- tmp.df %>% group_by(tissue, cycling) %>% summarise(n = n()) %>%
	mutate(freq = n / sum(n))

ordered.df <- frac.df %>% dplyr::filter(cycling) %>% arrange(desc(freq))
order.tissues.v <- as.character(ordered.df$tissue)


frac.df$tissue <- factor(frac.df$tissue, levels = order.tissues.v)


pct.p <- ggplot(frac.df, aes(x = tissue, y = freq, color = cycling)) +
	geom_point() +
	scale_color_brewer(palette = "Set1", name = "Actively proliferating") +
	geom_hline(yintercept = 0.5, size = 0.3, linetype = "dashed", alpha = 0.7) +
	scale_x_discrete(labels = str_c(levels(frac.df$tissue), "\n(n=", table(tmp.df$tissue), ")"), name = "") +
	labs(title = "Percentage of actively proliferating cells in human fetal tissues atlas", y = "Percentage") + 
	ylim(c(0, 1)) +
	theme(axis.text.x = element_text(size = 5,  vjust = 0.5, hjust = 0.5, angle = 30))

save_plot(here::here("figs", "sfigs", "sfig.fetalcyclingpct.pdf"), pct.p,
					base_height = 2, base_width = 2*1.2, nrow = 1.5, ncol = 2.2, device = cairo_pdf)
















