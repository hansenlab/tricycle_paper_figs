rm(list=ls())
source(here::here("scripts/utils.R"))


cell.lp  <- lapply(c("Intestine", "Kidney", "Pancreas", "Stomach"), function(o) {
	sce.o <- qread(here::here("data/fetal", str_c(o, ".qs")))
	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")

	p <- plot_emb_circle_scale(sce.o, dimred = 1, fig.title = str_c(o, " (n=", ncol(sce.o), ")"), x_lab = "UMAP-1", y_lab = "UMAP-2") +
		annotate(geom = "point", x = 0, y = 0, size = 1)
	
	return(p)
})

mp <- plot_grid(plotlist = cell.lp,
								nrow = 2, ncol = 2,  label_size = 10, labels = NULL, align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.test.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 2, device = cairo_pdf)








# getDenPlots <- function(theta, color_var.v, color_name, title, bw = 30, colors.v = NULL, scale = FALSE, line.alpha = 1) {
# 	d <- density(circular(theta), bw = bw)
# 	all.df <- data.frame(x = as.numeric(d$x), y = d$y)
# 	color_var.v[!(color_var.v %in% names(sort(table(color_var.v), decreasing = TRUE))[1:9])] <- NA
# 	color_var.v <- factor(color_var.v)
# 	if (scale) {
# 		alpha <- table(color_var.v) / sum(!is.na(color_var.v), na.rm = TRUE)
# 	} else {
# 		alpha <- rep(1, length(color_var.v))
# 	}
# 	
# 	
# 	stage.df <- do.call(rbind, lapply(seq_len(nlevels(color_var.v)), function(idx) {
# 		d <- density(circular(theta[which(color_var.v == levels(color_var.v)[idx])]), bw = bw)
# 		return(data.frame(x = as.numeric(d$x), y = d$y * alpha[idx] , color = levels(color_var.v)[idx]))
# 	}))
# 	stage.df$color <- factor(stage.df$color, levels = levels(color_var.v))
# 	max.v <- max(all.df$y, stage.df$y)
# 	
# 	if (is.null(colors.v)) colors.v <- brewer.pal(nlevels(color_var.v), "Set1")
# 	
# 	l.den.p <- ggplot(stage.df, aes(x = x  , y = y )) +
# 		geom_path(aes(color = color), size = 0.5, alpha = line.alpha) +
# 		geom_path(data = all.df, size = 0.5, alpha = line.alpha, color = "black", linetype = "dashed") +
# 		# geom_jitter(data = tmp.df, aes(x = theta1 * pi + pi / 2, y = tmp.l[[1]]$max * 0.5, color = cc), size = 0.1,  alpha = 0.5, shape = 16) + 
# 		scale_color_manual(values = colors.v, name = color_name, labels =  str_c(levels(color_var.v), "\nn=", table(color_var.v))) +
# 		scale_x_continuous(limits = c(0, 2 * pi),
# 											 breaks =  c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
# 											 labels = str_c(c(0, 0.5, 1, 1.5, 2), "\u03C0"),
# 											 name = "\u03B8") +
# 		labs(title = title, y = "Density") 
# 	
# 	cir.den.p <- ggplot(stage.df, aes(x = x  , y = y + max.v )) +
# 		geom_path(aes(color = color), size = 0.5, alpha = line.alpha) +
# 		geom_path(data = all.df, size = 0.5, alpha = line.alpha, color = "black", linetype = "dashed") +
# 		# geom_jitter(data = tmp.df, aes(x = theta1 * pi + pi / 2, y = tmp.l[[1]]$max * 0.5, color = cc), size = 0.1,  alpha = 0.5, shape = 16) + 
# 		scale_color_manual(values = colors.v, name = color_name, labels =  str_c(levels(color_var.v), "\nn=", table(color_var.v))) +
# 		coord_polar(theta = "x", start = - pi / 2, direction = -1, clip = "on") +
# 		scale_x_continuous(limits = c(0, 2 * pi),
# 											 breaks =  c(0, pi / 2, pi, 3 * pi / 2),
# 											 labels = str_c(c(0, 0.5, 1, 1.5), "\u03C0"),
# 											 name = "") +
# 		scale_y_continuous(limits = c(0, max.v * 2),
# 											 breaks =  c(max.v, max.v * 1.5, max.v * 2),
# 											 labels = c("0", format(max.v * c(0.5, 1), digits = 3)),
# 											 name = "Density") +
# 		labs(title = title) 
# 	
# 	return(list(linear = l.den.p, circular = cir.den.p))
# 	
# }


cell.lp  <- lapply(c("Intestine", "Kidney", "Pancreas", "Stomach"), function(o) {
	sce.o <- qread(here::here("data/fetal", str_c(o, ".qs")))
	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
	reducedDim(sce.o, "umap") <- data.frame(umap1 = sce.o$Main_cluster_umap_1, umap2 = sce.o$Main_cluster_umap_2)
	p1 <- plot_emb_circle_scale(sce.o, dimred = 2, fig.title = str_c(o, " (n=", ncol(sce.o), ")"), x_lab = "UMAP-1", y_lab = "UMAP-2")
	
	tmp.df <- data.frame(x = reducedDim(sce.o, "tricycleEmbedding")[, 1], y = reducedDim(sce.o, "tricycleEmbedding")[, 2], color = factor(sce.o$cc))
	tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
	scale_color <- scale_color_manual(values = c(ccColors.v, "grey"), name = "CC Stage", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA"))
	
	p2 <- plot_emb_circle_scale(sce.o, dimred = 1, fig.title = str_c(o, " (n=", ncol(sce.o), ")"), x_lab = px_lab, y_lab =  py_lab) +
		annotate(geom = "point", x = 0, y = 0, size = 1)
		# ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		# geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = 2.01,  alpha = 0.6) +
		# geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = 2.01,  alpha = 0.6) +
		# scale_color +
		# guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		# labs( x = px_lab, y = py_lab, title = str_c(o, " (n=", ncol(sce.o), ")")) +
		# theme(legend.position = c(0, 1),
		# 			legend.justification = c(0, 1),
		# 			legend.key = element_blank(),
		# 			legend.key.size = unit(6.5, "pt"))
	
	if (nlevels(factor(sce.o$Development_day)) > 9) {
		sel.ct.v <- names(sort(table(sce.o$Development_day), decreasing = TRUE))[1:9]
		sce.o$Development_day[!(sce.o$Development_day %in% sel.ct.v)] <- NA
	}
	tmp.df <- data.frame(x = sce.o$Main_cluster_umap_1, y = sce.o$Main_cluster_umap_2, color = factor(sce.o$Development_day), tricyclePosition = sce.o$tricyclePosition)
	tmp.df$cycling <- (tmp.df$tricyclePosition >= 0.25 * pi) & (tmp.df$tricyclePosition <= 1.5 * pi)
	p3 <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(pointsize = 2.01,  alpha = 0.6) +
		scale_color_brewer(palette = "Set1", name = "Day", labels = str_c(levels(factor(sce.o$Development_day)), " (n=", table(factor(sce.o$Development_day)), ")"))+
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( x = "UMAP-1", y = "UMAP-2", title = str_c(o, " (n=", ncol(sce.o), ")")) 

	frac.df <- tmp.df %>% group_by(color, cycling) %>% summarise(n = n()) %>%
		mutate(freq = n / sum(n))
	ordered.df <- frac.df %>% dplyr::filter(cycling) %>% arrange(desc(freq))
	order.color.v <- as.character(ordered.df$color)
	
	frac.df$color <- factor(frac.df$color, levels = as.character(sort(as.numeric(levels(frac.df$color)))))
	
	pct1.p <- ggplot(frac.df, aes(x = color, y = freq, color = cycling)) +
		geom_point(shape = 1, size = 1) +
		scale_color_brewer(palette = "Set1", name = "Actively proliferating") +
		geom_hline(yintercept = 0.5, size = 0.3, linetype = "dashed", alpha = 0.7) +
		scale_x_discrete(labels = str_c(levels(frac.df$color), "\n(n=", table(tmp.df$color), ")"), name = "Development day") +
		labs(title = str_c("Percentage of actively \nproliferating cells in ", tolower(o), " (n=", ncol(sce.o),")" ),  y = "Percentage") + 
		ylim(c(0, 1)) +
		theme(axis.text.x = element_text(size = 4.5,  vjust = 0.5, hjust = 0.5, angle = 30),
					legend.position = c(0, 1),
					legend.justification = c(0, 1),
					legend.direction="horizontal")
	
	
	
	
	if (nlevels(factor(sce.o$Main_cluster_name)) > 9) {
		sel.ct.v <- names(sort(table(sce.o$Main_cluster_name), decreasing = TRUE))[1:9]
		sce.o$Main_cluster_name[!(sce.o$Main_cluster_name %in% sel.ct.v)] <- NA
	}
	tmp.df <- data.frame(x = sce.o$Main_cluster_umap_1, y = sce.o$Main_cluster_umap_2, color = factor(sce.o$Main_cluster_name), tricyclePosition = sce.o$tricyclePosition)
	tmp.df$cycling <- (tmp.df$tricyclePosition >= 0.25 * pi) & (tmp.df$tricyclePosition <= 1.5 * pi)
	p4 <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(pointsize = 2.01,  alpha = 0.6) +
		scale_color_brewer(palette = "Set1", name = "Cell Type", labels = str_c(levels(factor(sce.o$Main_cluster_name)), " (n=", table(factor(sce.o$Main_cluster_name)), ")"))+
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( x = "UMAP-1", y = "UMAP-2", title = str_c(o, " (n=", ncol(sce.o), ")")) 
	
	frac.df <- tmp.df %>% group_by(color, cycling) %>% summarise(n = n()) %>%
		mutate(freq = n / sum(n))
	ordered.df <- frac.df %>% dplyr::filter(cycling) %>% arrange(desc(freq))
	order.color.v <- as.character(ordered.df$color)
	
	frac.df$color <- factor(frac.df$color, levels = order.color.v)
	
	pct2.p <- ggplot(frac.df %>% dplyr::filter(!is.na(color)), aes(x = color, y = freq, color = cycling)) +
		geom_point(shape = 1, size = 1) +
		scale_color_brewer(palette = "Set1", name = "Actively proliferating") +
		geom_hline(yintercept = 0.5, size = 0.3, linetype = "dashed", alpha = 0.7) +
		scale_x_discrete(labels = str_c(levels(frac.df$color), "\n(n=", table(tmp.df$color), ")"), name = "Cell type") +
		labs(title = str_c("Percentage of actively \nproliferating cells in ", tolower(o), " (n=", ncol(sce.o),")" ),  y = "Percentage") + 
		ylim(c(0, 1)) +
		theme(axis.text.x = element_text(size = 4.5,  vjust = 0.5, hjust = 0.5, angle = 30),
					axis.title.x = element_blank(),
					legend.position = c(0, 1),
					legend.justification = c(0, 1),
					legend.direction="horizontal")

	
	mp <- plot_grid(p2, p1, p4 + theme(legend.position = "none"), get_legend(p4), pct2.p,
									p3 + theme(legend.position = "none"), get_legend(p3), pct1.p,
									nrow = 1, ncol = 8,  rel_widths = c(1, 1, 1, 0.8, 1, 1, 0.5, 1),label_size = 10, labels = NULL, align = "v", axis = "lr")
	
	return(mp)
})
# save_plot(here::here("figs", "sfigs", "sfig.pdf"), mp,
# 					base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 7.5, device = cairo_pdf)

mp <- plot_grid(plotlist = cell.lp,
								nrow = 4, ncol = 1,  label_size = 10, labels = "auto", align = "hv", axis = "tblr")
save_plot(here::here("figs", "sfigs", "sfig.fetalsinglecell.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 4, ncol = 7.5, device = cairo_pdf)



nuclei.lp <- lapply(c("Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart", "Liver", "Lung", "Muscle", "Placenta", "Spleen", "Thymus"), function(o) {
	sce.o <- qread(here::here("data/fetal/nuclei", str_c(o, ".qs")))
	point.size <- ifelse(ncol(sce.o) > 1000000, 1.01, 2.01)
	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
	reducedDim(sce.o, "umap") <- data.frame(umap1 = sce.o$Main_cluster_umap_1, umap2 = sce.o$Main_cluster_umap_2)
	p1 <- plot_emb_circle_scale(sce.o, dimred = 2, fig.title = str_c(o, " (n=", ncol(sce.o), ")"), x_lab = "UMAP-1", y_lab = "UMAP-2")
	
	tmp.df <- data.frame(x = reducedDim(sce.o, "tricycleEmbedding")[, 1], y = reducedDim(sce.o, "tricycleEmbedding")[, 2], color = factor(sce.o$cc))
	tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
	scale_color <- scale_color_manual(values = c(ccColors.v, "grey"), name = "CC Stage", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA"))
	
	p2 <- plot_emb_circle_scale(sce.o, dimred = 1, fig.title = str_c(o, " (n=", ncol(sce.o), ")"), x_lab = px_lab, y_lab =  py_lab) +
		annotate(geom = "point", x = 0, y = 0, size = 1)
	# ggplot(tmp.df, aes(x = x, y = y, color = color)) +
	# geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = 2.01,  alpha = 0.6) +
	# geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = 2.01,  alpha = 0.6) +
	# scale_color +
	# guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
	# labs( x = px_lab, y = py_lab, title = str_c(o, " (n=", ncol(sce.o), ")")) +
	# theme(legend.position = c(0, 1),
	# 			legend.justification = c(0, 1),
	# 			legend.key = element_blank(),
	# 			legend.key.size = unit(6.5, "pt"))
	
	if (nlevels(factor(sce.o$Development_day)) > 9) {
		sel.ct.v <- names(sort(table(sce.o$Development_day), decreasing = TRUE))[1:9]
		sce.o$Development_day[!(sce.o$Development_day %in% sel.ct.v)] <- NA
	}
	tmp.df <- data.frame(x = sce.o$Main_cluster_umap_1, y = sce.o$Main_cluster_umap_2, color = factor(sce.o$Development_day), tricyclePosition = sce.o$tricyclePosition)
	tmp.df$cycling <- (tmp.df$tricyclePosition >= 0.25 * pi) & (tmp.df$tricyclePosition <= 1.5 * pi)
	p3 <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(pointsize = 2.01,  alpha = 0.6) +
		scale_color_brewer(palette = "Set1", name = "Day", labels = str_c(levels(factor(sce.o$Development_day)), " (n=", table(factor(sce.o$Development_day)), ")"))+
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( x = "UMAP-1", y = "UMAP-2", title = str_c(o, " (n=", ncol(sce.o), ")")) 
	
	frac.df <- tmp.df %>% group_by(color, cycling) %>% summarise(n = n()) %>%
		mutate(freq = n / sum(n))
	ordered.df <- frac.df %>% dplyr::filter(cycling) %>% arrange(desc(freq))
	order.color.v <- as.character(ordered.df$color)
	
	frac.df$color <- factor(frac.df$color, levels = as.character(sort(as.numeric(levels(frac.df$color)))))
	
	pct1.p <- ggplot(frac.df, aes(x = color, y = freq, color = cycling)) +
		geom_point(shape = 1, size = 1) +
		scale_color_brewer(palette = "Set1", name = "Actively proliferating") +
		geom_hline(yintercept = 0.5, size = 0.3, linetype = "dashed", alpha = 0.7) +
		scale_x_discrete(labels = str_c(levels(frac.df$color), "\n(n=", table(tmp.df$color), ")"), name = "Development day") +
		labs(title = str_c("Percentage of actively \nproliferating cells in ", tolower(o), " (n=", ncol(sce.o),")" ),  y = "Percentage") + 
		ylim(c(0, 1)) +
		theme(axis.text.x = element_text(size = 4.5,  vjust = 0.5, hjust = 0.5, angle = 30),
					legend.position = c(0, 1),
					legend.justification = c(0, 1),
					legend.direction="horizontal")
	
	
	
	
	if (nlevels(factor(sce.o$Main_cluster_name)) > 9) {
		sel.ct.v <- names(sort(table(sce.o$Main_cluster_name), decreasing = TRUE))[1:9]
		sce.o$Main_cluster_name[!(sce.o$Main_cluster_name %in% sel.ct.v)] <- NA
	}
	tmp.df <- data.frame(x = sce.o$Main_cluster_umap_1, y = sce.o$Main_cluster_umap_2, color = factor(sce.o$Main_cluster_name), tricyclePosition = sce.o$tricyclePosition)
	tmp.df$cycling <- (tmp.df$tricyclePosition >= 0.25 * pi) & (tmp.df$tricyclePosition <= 1.5 * pi)
	p4 <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(pointsize = 2.01,  alpha = 0.6) +
		scale_color_brewer(palette = "Set1", name = "Cell Type", labels = str_c(levels(factor(sce.o$Main_cluster_name)), " (n=", table(factor(sce.o$Main_cluster_name)), ")"))+
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( x = "UMAP-1", y = "UMAP-2", title = str_c(o, " (n=", ncol(sce.o), ")")) 
	
	frac.df <- tmp.df %>% group_by(color, cycling) %>% summarise(n = n()) %>%
		mutate(freq = n / sum(n))
	ordered.df <- frac.df %>% dplyr::filter(cycling) %>% arrange(desc(freq))
	order.color.v <- as.character(ordered.df$color)
	
	frac.df$color <- factor(frac.df$color, levels = order.color.v)
	
	pct2.p <- ggplot(frac.df %>% dplyr::filter(!is.na(color)), aes(x = color, y = freq, color = cycling)) +
		geom_point(shape = 1, size = 1) +
		scale_color_brewer(palette = "Set1", name = "Actively proliferating") +
		geom_hline(yintercept = 0.5, size = 0.3, linetype = "dashed", alpha = 0.7) +
		scale_x_discrete(labels = str_c(levels(frac.df$color), "\n(n=", table(tmp.df$color), ")"), name = "Cell type") +
		labs(title = str_c("Percentage of actively \nproliferating cells in ", tolower(o), " (n=", ncol(sce.o),")" ),  y = "Percentage") + 
		ylim(c(0, 1)) +
		theme(axis.text.x = element_text(size = 4.5,  vjust = 0.5, hjust = 0.5, angle = 30),
					axis.title.x = element_blank(),
					legend.position = c(0, 1),
					legend.justification = c(0, 1),
					legend.direction="horizontal")

	

	
	mp <- plot_grid(p2, p1, p4 + theme(legend.position = "none"), get_legend(p4), pct2.p,
									p3 + theme(legend.position = "none"), get_legend(p3), pct1.p,
									nrow = 1, ncol = 8,  rel_widths = c(1, 1, 1, 0.8, 1, 1, 0.5, 1),label_size = 10, labels = NULL, align = "v", axis = "lr")
	
	rm(sce.o)
	gc()
	message(o)
	return(mp)
})



mp1 <- plot_grid(plotlist = nuclei.lp[1:6],
								nrow = 6, ncol = 1,  label_size = 10, labels = "auto", align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.fetalsinglenuclei1.pdf"), mp1,
					base_height = 2, base_width = 2*1.2, nrow = 6, ncol = 7.5, device = cairo_pdf)


mp2 <- plot_grid(plotlist = nuclei.lp[7:12],
								 nrow = 5, ncol = 1,  label_size = 10, labels = letters[7:12], align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.fetalsinglenuclei2.pdf"), mp2,
					base_height = 2, base_width = 2*1.2, nrow = 5, ncol = 7.5, device = cairo_pdf)


# mp <- plot_grid(plotlist = nuclei.lp,
# 								nrow = 6, ncol = 2,  label_size = 10, labels = "auto", align = "hv", axis = "tblr")
# save_plot(here::here("figs", "sfigs", "sfig.fetalsinglenuclei2.pdf"), mp,
# 					base_height = 2, base_width = 2*1.2, nrow = 6, ncol = 5.5 *2 , device = cairo_pdf)


# Cerebellum.lp <- lapply(c("Cerebellum"), function(o) {
# 	sce.o <- qread(here::here("data/fetal/nuclei", str_c(o, ".qs")))
# 	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
# 	reducedDim(sce.o, "umap") <- data.frame(umap1 = sce.o$Main_cluster_umap_1, umap2 = sce.o$Main_cluster_umap_2)
# 	p1 <- plot_emb_circle_scale(sce.o, dimred = 2, fig.title = str_c(o, " (n=", ncol(sce.o), ")"), x_lab = "UMAP-1", y_lab = "UMAP-2", point.size = 1.01)
# 	
# 	tmp.df <- data.frame(x = reducedDim(sce.o, "tricycleEmbedding")[, 1], y = reducedDim(sce.o, "tricycleEmbedding")[, 2], color = factor(sce.o$cc))
# 	tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
# 	scale_color <- scale_color_manual(values = c(ccColors.v, "grey"), name = "CC Stage", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA"))
# 	
# 	p2 <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
# 		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = 1.01,  alpha = 0.6) +
# 		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = 1.01,  alpha = 0.6) +
# 		scale_color +
# 		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
# 		labs( x = px_lab, y = py_lab, title = str_c(o, " (n=", ncol(sce.o), ")")) +
# 		theme(legend.position = c(0, 1),
# 					legend.justification = c(0, 1),
# 					legend.key = element_blank(),
# 					legend.key.size = unit(6.5, "pt"))
# 	
# 	if (nlevels(factor(sce.o$Development_day)) > 9) {
# 		sel.ct.v <- names(sort(table(sce.o$Development_day), decreasing = TRUE))[1:9]
# 		sce.o$Development_day[!(sce.o$Development_day %in% sel.ct.v)] <- NA
# 	}
# 	tmp.df <- data.frame(x = sce.o$Main_cluster_umap_1, y = sce.o$Main_cluster_umap_2, color = factor(sce.o$Development_day))
# 	p3 <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
# 		geom_scattermore(pointsize = 1.01,  alpha = 0.6) +
# 		scale_color_brewer(palette = "Set1", name = "Day", labels = str_c(levels(factor(sce.o$Development_day)), " (n=", table(factor(sce.o$Development_day)), ")"))+
# 		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
# 		labs( x = "UMAP-1", y = "UMAP-2", title = str_c(o, " (n=", ncol(sce.o), ")")) 
# 	
# 	if (nlevels(factor(sce.o$Main_cluster_name)) > 9) {
# 		sel.ct.v <- names(sort(table(sce.o$Main_cluster_name), decreasing = TRUE))[1:9]
# 		sce.o$Main_cluster_name[!(sce.o$Main_cluster_name %in% sel.ct.v)] <- NA
# 	}
# 	tmp.df <- data.frame(x = sce.o$Main_cluster_umap_1, y = sce.o$Main_cluster_umap_2, color = factor(sce.o$Main_cluster_name))
# 	p4 <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
# 		geom_scattermore(pointsize = 1.01,  alpha = 0.6) +
# 		scale_color_brewer(palette = "Set1", name = "Cell Type", labels = str_c(levels(factor(sce.o$Main_cluster_name)), " (n=", table(factor(sce.o$Main_cluster_name)), ")"))+
# 		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
# 		labs( x = "UMAP-1", y = "UMAP-2", title = str_c(o, " (n=", ncol(sce.o), ")")) 
# 	
# 	mp <- plot_grid(p2, p1, p4 + theme(legend.position = "none"), get_legend(p4), p3 + theme(legend.position = "none"), get_legend(p3), 
# 									nrow = 1, ncol = 6,  rel_widths = c(1, 1, 1, 0.8, 1, 0.5),label_size = 10, labels = NULL, align = "hv", axis = "tblr")
# 	
# 	return(mp)
# })
# 
# save_plot(here::here("figs", "sfigs", "sfig.fetalCerebellum.pdf"), Cerebellum.lp[[1]],
# 					base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 5.5, device = cairo_pdf)
# 
# 
# Cerebrum.lp <- lapply(c("Cerebrum"), function(o) {
# 	sce.o <- qread(here::here("data/fetal/nuclei", str_c(o, ".qs")))
# 	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
# 	reducedDim(sce.o, "umap") <- data.frame(umap1 = sce.o$Main_cluster_umap_1, umap2 = sce.o$Main_cluster_umap_2)
# 	p1 <- plot_emb_circle_scale(sce.o, dimred = 2, fig.title = str_c(o, " (n=", ncol(sce.o), ")"), x_lab = "UMAP-1", y_lab = "UMAP-2", point.size = 1.01)
# 	
# 	tmp.df <- data.frame(x = reducedDim(sce.o, "tricycleEmbedding")[, 1], y = reducedDim(sce.o, "tricycleEmbedding")[, 2], color = factor(sce.o$cc))
# 	tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
# 	scale_color <- scale_color_manual(values = c(ccColors.v, "grey"), name = "CC Stage", labels =  c(ccLabels.v, "NA"), limits =   c(ccLabels.v, "NA"))
# 	
# 	p2 <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
# 		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = 1.01,  alpha = 0.6) +
# 		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = 1.01,  alpha = 0.6) +
# 		scale_color +
# 		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
# 		labs( x = px_lab, y = py_lab, title = str_c(o, " (n=", ncol(sce.o), ")")) +
# 		theme(legend.position = c(0, 1),
# 					legend.justification = c(0, 1),
# 					legend.key = element_blank(),
# 					legend.key.size = unit(6.5, "pt"))
# 	
# 	if (nlevels(factor(sce.o$Development_day)) > 9) {
# 		sel.ct.v <- names(sort(table(sce.o$Development_day), decreasing = TRUE))[1:9]
# 		sce.o$Development_day[!(sce.o$Development_day %in% sel.ct.v)] <- NA
# 	}
# 	tmp.df <- data.frame(x = sce.o$Main_cluster_umap_1, y = sce.o$Main_cluster_umap_2, color = factor(sce.o$Development_day))
# 	p3 <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
# 		geom_scattermore(pointsize = 1.01,  alpha = 0.6) +
# 		scale_color_brewer(palette = "Set1", name = "Day", labels = str_c(levels(factor(sce.o$Development_day)), " (n=", table(factor(sce.o$Development_day)), ")"))+
# 		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
# 		labs( x = "UMAP-1", y = "UMAP-2", title = str_c(o, " (n=", ncol(sce.o), ")")) 
# 	
# 	if (nlevels(factor(sce.o$Main_cluster_name)) > 9) {
# 		sel.ct.v <- names(sort(table(sce.o$Main_cluster_name), decreasing = TRUE))[1:9]
# 		sce.o$Main_cluster_name[!(sce.o$Main_cluster_name %in% sel.ct.v)] <- NA
# 	}
# 	tmp.df <- data.frame(x = sce.o$Main_cluster_umap_1, y = sce.o$Main_cluster_umap_2, color = factor(sce.o$Main_cluster_name))
# 	p4 <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
# 		geom_scattermore(pointsize = 1.01,  alpha = 0.6) +
# 		scale_color_brewer(palette = "Set1", name = "Cell Type", labels = str_c(levels(factor(sce.o$Main_cluster_name)), " (n=", table(factor(sce.o$Main_cluster_name)), ")"))+
# 		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
# 		labs( x = "UMAP-1", y = "UMAP-2", title = str_c(o, " (n=", ncol(sce.o), ")")) 
# 	
# 	mp <- plot_grid(p2, p1, p4 + theme(legend.position = "none"), get_legend(p4), p3 + theme(legend.position = "none"), get_legend(p3), 
# 									nrow = 1, ncol = 6,  rel_widths = c(1, 1, 1, 0.8, 1, 0.5),label_size = 10, labels = NULL, align = "hv", axis = "tblr")
# 	
# 	return(mp)
# })
# 
# save_plot(here::here("figs", "sfigs", "sfig.fetalCerebrum.pdf"), Cerebrum.lp[[1]],
# 					base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 5.5, device = cairo_pdf)





