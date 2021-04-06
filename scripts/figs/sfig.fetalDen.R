rm(list=ls())
source(here::here("scripts/utils.R"))


### plot the total UMIs
cell.df <- do.call(rbind, lapply(c("Intestine", "Kidney", "Pancreas", "Stomach"), function(o) {
	sce.o <- qread(here::here("data/fetal", str_c(o, ".qs")))
	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
	return(data.frame(tissue = o, totalUMIs = sce.o$TotalUMIs, tricyclePosition = sce.o$tricyclePosition))
}))

nuclei.df <- do.call(rbind, lapply(c("Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart", "Liver", "Lung", "Muscle", "Placenta", "Spleen", "Thymus"), function(o) {
	sce.o <- qread(here::here("data/fetal/nuclei", str_c(o, ".qs"))) 
	point.size <- ifelse(ncol(sce.o) > 1000000, 1.01, 2.01)
	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
	
	return(data.frame(tissue = o, totalUMIs = sce.o$TotalUMIs, tricyclePosition = sce.o$tricyclePosition))
}))


tmp.df <- rbind(cell.df, nuclei.df)
tmp.df$tissue <- factor(tmp.df$tissue, levels = c("Intestine", "Kidney", "Pancreas", "Stomach", "Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart", "Liver", "Lung", "Muscle", "Placenta", "Spleen", "Thymus"))


d <- density(circular(tmp.df$tricyclePosition), bw = 30)
all.df <- data.frame(x = as.numeric(d$x), y = d$y)

tissue.df <- do.call(rbind, lapply(seq_len(nlevels(tmp.df$tissue)), function(idx) {
	d <- density(circular(tmp.df$tricyclePosition[which(tmp.df$tissue == levels(tmp.df$tissue)[idx])]), bw = 30)
	return(data.frame(x = as.numeric(d$x), y = d$y , tissue = levels(tmp.df$tissue)[idx]))
}))
tissue.df$tissue <- factor(tissue.df$tissue, levels = levels(tmp.df$tissue))
max.v <- max(all.df$y, tissue.df$y)


l.den.p <- ggplot(tissue.df, aes(x = x  , y = y )) +
	geom_path(aes(color = tissue), size = 0.3, alpha = 1) +
	geom_path(data = all.df, size = 0.3, alpha = 0.5, color = "black", linetype = "dashed") +
	# geom_jitter(data = tmp.df, aes(x = theta1 * pi + pi / 2, y = tmp.l[[1]]$max * 0.5, color = cc), size = 0.1,  alpha = 0.5, shape = 16) + 
	scale_color_manual(values = randomcoloR::distinctColorPalette(k = nlevels(tmp.df$tissue)), name = "", labels =  str_c(levels(tmp.df$tissue), " (n=", table(tmp.df$tissue), ")")) +
	scale_x_continuous(limits = c(0, 2 * pi),
										 breaks =  c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
										 labels = str_c(c(0, 0.5, 1, 1.5, 2), "\u03C0"),
										 name = theta_lab) +
	labs(title = str_c("Human tissue atlas (n=", nrow(tmp.df), ")"), y = "Density") 


save_plot(here::here("figs", "sfigs", "sfig.fetalDen.pdf"), l.den.p,
					base_height = 2, base_width = 2*1.2, nrow = 1.2, ncol = 1.8 , device = cairo_pdf)






#Adrenal Cerebellum   Cerebrum        Eye      Heart  Intestine     Kidney      Liver       Lung     Muscle   Pancreas   Placenta     Spleen    Stomach     Thymus 
#489        796        782        408        409        503        759        452        431        413        892        368        372        429        354

