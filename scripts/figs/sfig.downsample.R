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

p.l <- lapply(c(1, 0.8, 0.6, 0.4, 0.2, 0.1), function(prop) {
	set.seed(prop)
	

	counts.m <- downsampleMatrix(raw.m, prop = prop)
	median.v <- median(colSums(counts.m))
	log.m <- normalizeCounts(counts.m, log=TRUE)
	reducedDim(hipp.o, "new") <- project_cycle_space(log.m, ref.m = ref.m)
	if (prop == 1) {
		scat.p <- plotScatCC(sce.o = hipp.o, dimred = "tricycleEmbedding", x_lab = px_lab, y_lab = py_lab, title = str_c(metadata(hipp.o)$dataname, "\n(original; lib.size median:",format(median.v, digits = 2), ")")) + theme(legend.position = "none")
		
		
		tmp.df <- data.frame(x = reducedDim(hipp.o, "new")[, 1],
												 y = reducedDim(hipp.o, "new")[, 2],
												 theta = theta.v,
												 color = hipp.o$CCStage,
												 totalUMIs = log2(colSums(counts.m)))
		tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
		
		umi.p <- ggplot(tmp.df, aes(x = x, y = y, color = totalUMIs)) +
			geom_scattermore(pointsize = metadata(hipp.o)$point.size,  alpha = metadata(hipp.o)$point.alpha) +
			scale_color_viridis(name = "log2(lib.size)") +  
			labs( 	x = px_lab, 
						 y = py_lab,
						 title = str_c(metadata(hipp.o)$dataname, "\n(original; lib.size median:",format(median.v, digits = 2), ")")) +
			theme(legend.position = c(1, 1),
						legend.justification = c(1, 1))
		
		return(list(scat.p = scat.p, theta.p = NULL, umi.p = umi.p, proj = reducedDim(hipp.o, "tricycleEmbedding")[, 1:2], median.v = median.v))
	}
	
	scat.p <- plotScatCC(sce.o = hipp.o, dimred = "new", x_lab = px_lab, y_lab = py_lab, title = str_c(metadata(hipp.o)$dataname, "\n(downsampled; lib.size median:", format(median.v, digits = 2), ")")) + theme(legend.position = "none")
	
	tmp.df <- data.frame(x = reducedDim(hipp.o, "new")[, 1],
											 y = reducedDim(hipp.o, "new")[, 2],
											 theta = theta.v,
											 new = as.numeric(circular::coord2rad(reducedDim(hipp.o, "new")[, 1:2])),
											 color = hipp.o$CCStage,
											 totalUMIs = log2(colSums(counts.m)))
	tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
	
	
	theta.p <- ggplot(tmp.df, aes(x = theta, y = new, color = color)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize  = metadata(hipp.o)$point.size, alpha = metadata(hipp.o)$point.alpha) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize  = metadata(hipp.o)$point.size, alpha = metadata(hipp.o)$point.alpha) +
		scale_color_manual(values = c(ccColors.v, "grey"), name = "CC Stage", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA")) +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
		labs( 	x = bquote(paste('CC'['ns']," Position \u03B8 (original)")), 
					 y = bquote(paste('CC'['ns']," Position \u03B8 (downsampled)")),
					 title = str_c(metadata(hipp.o)$dataname, "\n(downsampled; lib.size median:", format(median.v, digits = 2), ")")) +
		annotate(geom = "text", x = 2, y = 5, size = 3, hjust = 0.5, vjust = 0.5,
						 label = str_c("Circular correlation \n \u03C1=", sprintf("%.3f", Directional::circ.cor1(tmp.df$theta, tmp.df$new , rads = T)[1])), parse = FALSE) +
		scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi , labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2 )* pi) +
		scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi , labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2)* pi) +
		theme(legend.position = "none")
	
	umi.p <- ggplot(tmp.df, aes(x = x, y = y, color = totalUMIs)) +
		geom_scattermore(pointsize = metadata(hipp.o)$point.size,  alpha = metadata(hipp.o)$point.alpha) +
		scale_color_viridis(name = "log2(lib.size)") +  
		labs( 	x = px_lab, 
					 y = py_lab,
					 title = str_c(metadata(hipp.o)$dataname, "\n(downsampled; lib.size median:", format(median.v, digits = 2), ")")) +
		theme(legend.position = c(1, 1),
						legend.justification = c(1, 1))
	
	
	return(list(scat.p = scat.p, theta.p = theta.p, umi.p = umi.p, proj = reducedDim(hipp.o, "new")[, 1:2], median.v = median.v))
	
})
p.l[[1]]$scat.p <- p.l[[1]]$scat.p + theme(legend.position = c(1, 1),
																																 legend.justification = c(1, 1))

### putting together
tmp.df <- do.call(rbind, lapply(seq_along(c(1, 0.8, 0.6, 0.4, 0.2, 0.1)), function(i) {
	df <- data.frame(p.l[[i]]$proj, s = c(1, 0.8, 0.6, 0.4, 0.2, 0.1)[i])
	df
}))
tmp.df$s <- factor(tmp.df$s, levels = c(1, 0.8, 0.6, 0.4, 0.2, 0.1))

p.l[[1]]$theta.p <- ggplot(tmp.df, aes(x = PC1, y = PC2, color = s)) +
	geom_scattermore(pointsize = metadata(hipp.o)$point.size,  alpha = metadata(hipp.o)$point.alpha) +
	scale_color_brewer(palette = "Set2", name = "lib.s median", labels = format(sapply(p.l, "[[", 5), digits = 2)) +
	labs( y = py_lab, x = px_lab, title = str_c(metadata(hipp.o)$dataname, "\n(multiple lib.size)")) +
	theme(legend.position = c(1, 1),
					legend.justification = c(1, 1))


mp <- plot_grid(plotlist = c(lapply(p.l, "[[", 1), lapply(p.l, "[[", 2), lapply(p.l, "[[", 3)),
								nrow = 3, ncol = 6, label_size = 10, labels = c("a", "b", rep("", 4), "c", "d", rep("", 4), "e", "f"))

save_plot(here::here("figs", "sfigs", "sfig.downsample.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 3, ncol = 6, device = cairo_pdf)

