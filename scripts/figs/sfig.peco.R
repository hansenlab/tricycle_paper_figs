rm(list=ls())
source(here::here("scripts/utils.R"))

library(ggrepel)

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))
endo.o <- qread(here::here("data/plotdata/endo.qs"))
mHSC.o <- qread(here::here("data/plotdata/mHSC.qs"))
HeLa1.o <- qread(here::here("data/plotdata/HeLa1.qs"))
HeLa2.o <- qread(here::here("data/plotdata/HeLa2.qs"))
hfIntestineSub.o <- qread(here::here("data/plotdata/hfIntestineSub.qs"))

hU2OS.o <- qread(here::here("data/plotdata/hU2OS.qs"))
hiPSCs.o <- qread(here::here("data/plotdata/hiPSCs.qs"))
mESC.o <- qread(here::here("data/plotdata/mESC.qs"))
hESC.o <- qread(here::here("data/plotdata/hESC.qs"))

gene_n.v <- sapply(list(neurosphere.o, hipp.o, endo.o, mHSC.o, HeLa1.o, HeLa2.o, hfIntestineSub.o, hU2OS.o, hiPSCs.o, mESC.o, hESC.o), function(x) length(metadata(x)$peco$g.v))
genes.v <- unlist(sapply(list(neurosphere.o, hipp.o, endo.o, mHSC.o, HeLa1.o, HeLa2.o, hfIntestineSub.o, hU2OS.o, hiPSCs.o, mESC.o, hESC.o), function(x) metadata(x)$peco$g.v))
length(genes.v != "UBC")

getR2 <- function(sce.o) {
	peco.l <- metadata(sce.o)$peco
	tmp.df <- data.frame(x = peco.l$fit_cyclic$cellcycle_peco_ordered, 
											 theta = colData(sce.o)[names(peco.l$fit_cyclic$cellcycle_peco_ordered), "tricyclePosition"])
	
	r2.df <- do.call(rbind, lapply(seq_along(peco.l$g.v), function(i) {
		if (metadata(sce.o)$species == "mouse") {
			if (peco.l$g.v[i] == "TOP2A") {
				col.outname <- "Top2A"
			} else if (peco.l$g.v[i] == "UBE2C") {
				col.outname <- "Ube2C"
			} else {
				col.outname <- str_to_title(peco.l$g.v[i])
			}
		} else {
			col.outname <- peco.l$g.v[i]
		} 
		
		tmp.df$y <- peco.l$yy_input[i, names(peco.l$fit_cyclic$cellcycle_peco_ordered)]
		r2theta <- fit_periodic_loess(theta.v = tmp.df$theta, y = tmp.df$y)$rsquared
		tmp.df$peco.fitted <- peco.l$fit_cyclic$cellcycle_function[[i]](tmp.df$theta)
		ss.total <- sum(scale(tmp.df$y, scale = FALSE) ^ 2)
		r2peco <- 1 - sum((tmp.df$y - tmp.df$peco.fitted) ^ 2) / ss.total
		
		
		return(data.frame(gene = col.outname, r2theta = r2theta, r2peco = r2peco))
	}))
	
	r2.df$data_name <- metadata(sce.o)$dataname
	return(r2.df)
}


tmp.df <- do.call(rbind, lapply(list(neurosphere.o, hipp.o, endo.o, mHSC.o, HeLa1.o, HeLa2.o, hfIntestineSub.o, hU2OS.o, hiPSCs.o, mESC.o, hESC.o), function(x) getR2(x)))

tmp.df$data_name <- factor(tmp.df$data_name, levels = unique(tmp.df$data_name))

r2.p <- ggplot(tmp.df, aes(x = r2peco, y = r2theta, label = gene)) +
	geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_abline(slope = 1, intercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point(aes(color = data_name), size = 1,  alpha = 1,  shape = 1) +
	geom_text_repel(data = tmp.df, size = 1.3, segment.size = 0.2, max.overlaps = 25) +
	scale_color_brewer(palette = "Set3", name = "Dataset") +
	labs( y = bquote(paste("tricycle R"^"2")), x = bquote(paste("peco R"^"2")), title = bquote(paste("Comparisons of R"^"2"))) 


plotPeco <- function(sce.o, color.var = "CCStage",  color.name = "SchwabeCC", colors.v = NULL, labels.v = NULL) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (is.null(colors.v)) colors.v <- ccColors.v
	if (is.null(labels.v)) labels.v <- ccLabels.v
	
	peco.l <- metadata(sce.o)$peco
	tmp.df <- data.frame(x = peco.l$fit_cyclic$cellcycle_peco_ordered, color = colData(sce.o)[names(peco.l$fit_cyclic$cellcycle_peco_ordered), color.var],
											 theta = colData(sce.o)[names(peco.l$fit_cyclic$cellcycle_peco_ordered), "tricyclePosition"])
	if (any(is.na(tmp.df$color))) {
		tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
		scale_color <- scale_color_manual(values = c(colors.v, "grey"), name = color.name, labels =  c(labels.v, "NA"), limits =   c(labels.v, "NA"))
	} else {
		scale_color <- scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = labels.v)
	}
	
	p.l <- lapply(which(peco.l$g.v == "TOP2A"), function(i) {
		if (metadata(sce.o)$species == "mouse") {
			if (peco.l$g.v[i] == "TOP2A") {
				col.outname <- "Top2A"
			} else if (peco.l$g.v[i] == "UBE2C") {
				col.outname <- "Ube2C"
			} else {
				col.outname <- str_to_title(peco.l$g.v[i])
			}
		} else {
			col.outname <- peco.l$g.v[i]
		} 
		
		tmp.df$y <- peco.l$yy_input[i, names(peco.l$fit_cyclic$cellcycle_peco_ordered)]
		r2theta <- fit_periodic_loess(theta.v = tmp.df$theta, y = tmp.df$y)$rsquared
		tmp.df$peco.fitted <- peco.l$fit_cyclic$cellcycle_function[[i]](tmp.df$theta)
		ss.total <- sum(scale(tmp.df$y, scale = FALSE) ^ 2)
		r2peco <- 1 - sum((tmp.df$y - tmp.df$peco.fitted) ^ 2) / ss.total
		if (r2theta > r2peco) {
			label_color <- c("red", "black")
		} else {
			label_color <- c("black", "red")
		}
		
		x <- seq(0,2*pi, length.out = 100)
		line.df <- data.frame(x = x, y = peco.l$fit_cyclic$cellcycle_function[[i]](x))
		
		
		sact.p <- ggplot(data = tmp.df, aes(x = x, y = y, color = color)) +
			geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size, alpha = point.alpha) +
			geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size, alpha = point.alpha) +
			geom_path(data = line.df, aes(x = x , y = y), linetype = "dashed", color = "black", size = 0.8, alpha = 0.6, inherit.aes = FALSE) +
			scale_color +
			guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
			labs( y = "peco normalized expression", x = "peco \u03B8", title = str_c(metadata(sce.o)$dataname, " peco (", col.outname, ")")) +
			scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2 * pi)) +
			theme(legend.position = "none") +
			annotate(geom = "text", x = pi, y = .percent_range(tmp.df$y, 0), size = 2, hjust = 0.5, vjust = 0, color = label_color[1],
							 label = as.character(as.expression(substitute('CC'['ns']~" Position \u03B8"~~italic(R)^2~"="~rsquared, list(rsquared = format(r2theta, digits = 3))))), parse = TRUE) +
			annotate(geom = "text", x = pi, y = .percent_range(tmp.df$y, 0.1), size = 2, hjust = 0.5, vjust = 0, color = label_color[2],
							 label = as.character(as.expression(substitute('peco \u03B8'~~italic(R)^2~"="~rsquared, list(rsquared = format(r2peco, digits = 3))))), parse = TRUE)
		
		return(sact.p)
		
	})
	return(p.l)
}

plotPeco2 <- function(sce.o, color.var = "CCStage") {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha

	
	peco.l <- metadata(sce.o)$peco
	tmp.df <- data.frame(x = peco.l$fit_cyclic$cellcycle_peco_ordered,
											 theta = colData(sce.o)[names(peco.l$fit_cyclic$cellcycle_peco_ordered), "tricyclePosition"])
	
	p.l <- lapply(which(peco.l$g.v == "TOP2A"), function(i) {
		if (metadata(sce.o)$species == "mouse") {
			if (peco.l$g.v[i] == "TOP2A") {
				col.outname <- "Top2A"
			} else if (peco.l$g.v[i] == "UBE2C") {
				col.outname <- "Ube2C"
			} else {
				col.outname <- str_to_title(peco.l$g.v[i])
			}
		} else {
			col.outname <- peco.l$g.v[i]
		} 
		
		tmp.df$y <- peco.l$yy_input[i, names(peco.l$fit_cyclic$cellcycle_peco_ordered)]
		r2theta <- fit_periodic_loess(theta.v = tmp.df$theta, y = tmp.df$y)$rsquared
		tmp.df$peco.fitted <- peco.l$fit_cyclic$cellcycle_function[[i]](tmp.df$theta)
		ss.total <- sum(scale(tmp.df$y, scale = FALSE) ^ 2)
		r2peco <- 1 - sum((tmp.df$y - tmp.df$peco.fitted) ^ 2) / ss.total
		if (r2theta > r2peco) {
			label_color <- c("red", "black")
		} else {
			label_color <- c("black", "red")
		}
		
		x <- seq(0,2*pi, length.out = 100)
		line.df <- data.frame(x = x, y = peco.l$fit_cyclic$cellcycle_function[[i]](x))
		
		
		sact.p <- ggplot(data = tmp.df, aes(x = x, y = y)) +
			geom_scattermore(pointsize = point.size, alpha = point.alpha, color = "grey") +
			geom_path(data = line.df, aes(x = x , y = y), linetype = "dashed", color = "black", size = 0.8, alpha = 0.6, inherit.aes = FALSE) +
			labs( y = "peco normalized expression", x = "peco \u03B8", title = str_c(metadata(sce.o)$dataname, " peco (", col.outname, ")")) +
			scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2 * pi)) +
			theme(legend.position = "none") +
			annotate(geom = "text", x = pi, y = .percent_range(tmp.df$y, 0), size = 2, hjust = 0.5, vjust = 0, color = label_color[1],
							 label = as.character(as.expression(substitute('CC'['ns']~" Position \u03B8"~~italic(R)^2~"="~rsquared, list(rsquared = format(r2theta, digits = 3))))), parse = TRUE) +
			annotate(geom = "text", x = pi, y = .percent_range(tmp.df$y, 0.1), size = 2, hjust = 0.5, vjust = 0, color = label_color[2],
							 label = as.character(as.expression(substitute('peco \u03B8'~~italic(R)^2~"="~rsquared, list(rsquared = format(r2peco, digits = 3))))), parse = TRUE)
		
		return(sact.p)
		
	})
	return(p.l)
}

peco1.pl <- do.call(c, lapply(list(neurosphere.o, hipp.o, endo.o, mHSC.o, HeLa1.o, HeLa2.o, hfIntestineSub.o, hU2OS.o, hiPSCs.o), function(x) plotPeco2(x)))
peco2.pl <- do.call(c, lapply(list(mESC.o, hESC.o), function(x) plotPeco2(x)))
# peco2.pl[[1]] <- peco2.pl[[1]]  + theme(legend.position = c(1,1 ),
# 																				legend.justification = c(1, 1),
# 																				legend.key = element_blank(),
# 																				legend.key.size = unit(6.5, "pt"))

peco.pl <- c(list(r2.p), peco1.pl, peco2.pl) 


mp <- plot_grid(plotlist = peco.pl,
								nrow = 4, ncol = 3, label_size = 10, labels = "auto")

save_plot(here::here("figs", "sfigs", "sfig.peco.pdf"), mp,
					base_height = 2, base_width = 2*1.5, nrow = 4, ncol = 3, device = cairo_pdf)


save_plot(here::here("figs", "response_figs", "response.r2.pdf"), r2.p,
					base_height = 2, base_width = 2*1.5, nrow = 1, ncol = 1, device = cairo_pdf)



