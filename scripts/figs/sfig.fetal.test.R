rm(list=ls())
source(here::here("scripts/utils.R"))

o <- "Lung"

sce.o <- qread(here::here("data/fetal/nuclei", str_c(o, ".qs")))
point.size <- ifelse(ncol(sce.o) > 1000000, 1.01, 2.01)
sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")

r <- sqrt(apply(reducedDim(sce.o, "tricycleEmbedding")[, 1:2], 1, function(x) sum(x ^ 2)))

retain.idx <- which(r > 0.2)

sce.o <- sce.o[, retain.idx]

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

frac.df$color <- factor(frac.df$color, levels = order.color.v)

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


save_plot(here::here("figs", "sfigs", "sfig.test.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 7.5, device = cairo_pdf)
