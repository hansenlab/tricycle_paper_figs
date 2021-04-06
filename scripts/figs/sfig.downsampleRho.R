rm(list=ls())
source(here::here("scripts/utils.R"))

library(scuttle)
neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))
hipp_full.o <- qread(here::here("data/hipp.qs"))
raw.m <- qread(here::here("data/hippRaw.m.qs"))

ref.m <- metadata(neurosphere.o)$rotation
theta.v <- hipp.o$tricyclePosition

sum(!rownames(hipp_full.o) %in% rownames(raw.m))
sum(!colnames(hipp_full.o) %in% colnames(raw.m))

raw.m <- raw.m[match(rownames(hipp_full.o), rownames(raw.m)), match(colnames(hipp_full.o), colnames(raw.m))]

tmp.df <- do.call(rbind, lapply(seq(from = 0.9, to = 0.1, by = -0.1), function(prop) {
	out.df <- do.call(rbind, bplapply(seq_len(30), function(i) {
		set.seed(i)
		counts.m <- downsampleMatrix(raw.m, prop = prop)
		median.v <- median(colSums(counts.m))
		log.m <- normalizeCounts(counts.m, log=TRUE)
		proj.m <- project_cycle_space(log.m, ref.m = ref.m)
    rho <- Directional::circ.cor1(theta.v, coord2rad(proj.m[, 1:2]), rads = T)[1]
    
		return(data.frame(median = median.v, rho = rho))
	}, BPPARAM = MulticoreParam(workers = 10)))
	out.df$prop <- prop
	
	return(out.df)
}))

rho.box.p <- ggplot(tmp.df, aes(x = factor(prop, levels = seq(from = 0.9, to = 0.1, by = -0.1)), y = rho)) +
	geom_boxplot(fill = "white", outlier.shape=NA, outlier.size = 0.2, size = 0.2, width = 0.8) +
	geom_quasirandom(width = 0.2, size = 0.2, shape = 16, alpha = 0.5) +
	labs( y = "Circular correlation \u03C1", 
				x = "Approx. libray size median", 
				title = str_c("Stability with libray size downsampled (mHippNPC)")) +
	scale_x_discrete(labels = c("9000", "8000", "7000", "6000", "5000", "4000", "3000", "2000", "1000")) +
	scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), limits = c(min(tmp.df$rho, 0), max(tmp.df$rho, 1)))

save_plot(here::here("figs", "sfigs", "sfig.downsampleRho.pdf"), rho.box.p,
					base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 1.2, device = cairo_pdf)












