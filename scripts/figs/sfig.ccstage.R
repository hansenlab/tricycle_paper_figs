rm(list=ls())
source(here::here("scripts/utils.R"))

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))
endo.o <- qread(here::here("data/plotdata/endo.qs"))
mRetina.o <- qread(here::here("data/plotdata/mRetina.qs"))
HeLa1.o <- qread(here::here("data/plotdata/HeLa1.qs"))
HeLa2.o <- qread(here::here("data/plotdata/HeLa2.qs"))


# plotThetaDen <- function(sce.o, title = NULL, bw = 30, scale = FALSE) {
# 	if (is.null(title)) title <- metadata(sce.o)$dataname
# 	theta <- sce.o$tricyclePosition
# 	cc <- fct_drop(sce.o$CCStage)
# 	d <- density(circular(theta), bw = bw)
# 	all.df <- data.frame(x = as.numeric(d$x), y = d$y)
# 	
# 	
# 	if (scale) {
# 		alpha <- table(cc) / sum(!is.na(cc), na.rm = TRUE)
# 	} else {
# 		alpha <- rep(1, 5)
# 	}
# 	
# 	stage.df <- do.call(rbind, lapply(seq_len(nlevels(cc)), function(idx) {
# 		d <- density(circular(theta[which(cc == levels(cc)[idx])]), bw = bw)
# 		return(data.frame(x = as.numeric(d$x), y = d$y * alpha[idx], cc = levels(cc)[idx]))
# 	}))
# 	stage.df$cc <- factor(stage.df$cc, levels = levels(cc))
# 	max.v <- max(all.df$y, stage.df$y)
# 	
# 	
# 	l.den.p <- ggplot(stage.df, aes(x = x  , y = y )) +
# 		geom_path(aes(color = cc), size = 0.5, alpha = 0.5) +
# 		geom_path(data = all.df, size = 0.5, alpha = 0.5, color = "black", linetype = "dashed") +
# 		# geom_jitter(data = tmp.df, aes(x = theta1 * pi + pi / 2, y = tmp.l[[1]]$max * 0.5, color = cc), size = 0.1,  alpha = 0.5, shape = 16) + 
# 		scale_color_manual(values = ccColors.v, name = "Cellcycle", labels =  ccLabels.v) +
# 		scale_x_continuous(limits = c(0, 2 * pi),
# 											 breaks =  c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
# 											 labels = str_c(c(0, 0.5, 1, 1.5, 2), "\u03C0"),
# 											 name = "\u03B8") +
# 		labs(title = title, y = "Density") 
# 	
# 	cir.den.p <- ggplot(stage.df, aes(x = x  , y = y + max.v )) +
# 		geom_path(aes(color = cc), size = 0.5, alpha = 0.5) +
# 		geom_path(data = all.df, size = 0.5, alpha = 0.5, color = "black", linetype = "dashed") +
# 		# geom_jitter(data = tmp.df, aes(x = theta1 * pi + pi / 2, y = tmp.l[[1]]$max * 0.5, color = cc), size = 0.1,  alpha = 0.5, shape = 16) + 
# 		scale_color_manual(values = ccColors.v, name = "Cellcycle", labels =  ccLabels.v) +
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

plotPct0Box <- function(sce.o, title = NULL) {
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " (n=", ncol(sce.o), ")")
	theta.v <- sce.o$tricyclePosition

	g0g1.v <- (theta.v < 0.25 * pi) | (theta.v > 1.5 * pi)
	
	tmp.df <- data.frame(pct0 = sce.o$pct0, cc = sce.o$CCStage, g0g1 = g0g1.v)
	tmp.df$cc <- fct_explicit_na(tmp.df$cc, "NA")
	
	p <- ggplot(tmp.df, aes(x = cc, y = pct0, color = cc,  dodge = g0g1, fill = g0g1)) +
		geom_boxplot(aes(fill = g0g1), outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.8) +
		scale_color_manual(values = c(ccColors.v, "grey60"), name = "Schwabe", labels = c(ccLabels.v, "NA")) +
		scale_fill_manual(values = c("white", "grey90"), name = "Flag") + 
		# geom_quasirandom(width = 0.08, size = 0.2, shape = 16, alpha = 0.5, dodge.width = 0.8) +
		scale_x_discrete(name = "SchwabeCC", labels = str_c(levels(tmp.df$cc), "\n(n=", table(tmp.df$cc), ")")) +
		labs( y =  "Pct of 0 (projection genes)",  title = title) +
		theme(axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5, angle = 30),
					axis.title.x = element_blank())
	return(p)
}

plotGeneBox <- function(sce.o, col.name, col.outname = NULL, title = NULL) {
	theta.v <- sce.o$tricyclePosition
	if (metadata(sce.o)$species == "mouse") {
		if (col.name == "top2a") {
			col.outname <- "Top2A"
		}  else {
			col.outname <- str_to_title(col.name)
		}
	} else {
		col.outname <- str_to_upper(col.name)
	} 
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " ", col.outname, " (n=", ncol(sce.o), ")")
	gene.v <- colData(sce.o)[[col.name]]
	g0g1.v <- (theta.v < 0.25 * pi) | (theta.v > 1.5 * pi)
	tmp.df <- data.frame(y = gene.v, cc = sce.o$CCStage, g0g1 = g0g1.v)
	tmp.df$cc <- fct_explicit_na(tmp.df$cc, "NA")
	
	p <- ggplot(tmp.df, aes(x = cc, y = y, color = cc,  dodge = g0g1, fill = g0g1)) +
		geom_boxplot(aes(fill = g0g1), outlier.shape=1, outlier.size = 0.2, size = 0.2, width = 0.8) +
		scale_color_manual(values = c(ccColors.v, "grey60"), name = "Schwabe", labels = c(ccLabels.v, "NA")) +
		scale_fill_manual(values = c("white", "grey90"), name = "Flag") + 
		# geom_quasirandom(width = 0.08, size = 0.2, shape = 16, alpha = 0.5, dodge.width = 0.8) +
		scale_x_discrete(name = "SchwabeCC", labels = str_c(levels(tmp.df$cc), "\n(n=", table(tmp.df$cc), ")")) +
		labs( y = bquote(paste('log'['2'],'(expression of ', .(col.outname), ")")),  title = title) +
		theme(axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5, angle = 30),
					axis.title.x = element_blank())
	return(p)
}


ccstage.pl <- do.call(c, lapply(list(neurosphere.o, hipp.o, endo.o, mRetina.o, HeLa1.o, HeLa2.o), function(sce.o) {
	pct.theta.p <- plotLoess(sce.o, col.name = "pct0", col.outname = "Pct0", y_lab = "Pct of 0 (projection genes)") + theme(legend.position = "none")
	den.p <- plotThetaDen(sce.o)$circular + theme(legend.position = "none")
	pct0.box.p <- plotPct0Box(sce.o) 
	top2a.box.p <- plotGeneBox(sce.o, "top2a") + theme(legend.position = "none")
	smc4.box.p <- plotGeneBox(sce.o, "smc4") + theme(legend.position = "none")
	legend.p <- get_legend(pct0.box.p)
	return(list(pct.theta.p, den.p, pct0.box.p + theme(legend.position = "none"), top2a.box.p, smc4.box.p, legend.p))
}))

mp <- plot_grid(plotlist  = ccstage.pl, 
	nrow = 6, ncol = 6, rel_widths = c(1, 1, rep(1, 3), 0.3), labels = do.call(c, lapply(letters[1:6], function(x) c(x, rep("", 5)))), align = "hv", axis = "lrtb")


save_plot(here::here("figs", "sfigs", "sfig.ccstage.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 6, ncol = 4.5, device = cairo_pdf)






#### use revelio cc labels 
rownames(metadata(neurosphere.o)$revelio$dc)
rownames(metadata(hipp.o)$revelio$dc)
rownames(metadata(endo.o)$revelio$dc)
rownames(metadata(HeLa1.o)$revelio$dc)
rownames(metadata(HeLa2.o)$revelio$dc)

useRevelioCC <- function(sce.o) {
	metadata(sce.o)$N <- ncol(sce.o)
	sce.o <- sce.o[, match(rownames(metadata(sce.o)$revelio$dc), colnames(sce.o))]
	sce.o$CCStage <- metadata(sce.o)$revelio$cc
	print(ncol(sce.o))
	sce.o
}

neurosphere.o <- useRevelioCC(neurosphere.o)
hipp.o <- useRevelioCC(hipp.o)
endo.o <- useRevelioCC(endo.o)
HeLa1.o <- useRevelioCC(HeLa1.o)
HeLa2.o <- useRevelioCC(HeLa2.o)

### mRetina.o  special
metadata(mRetina.o)$N <- ncol(mRetina.o)
mRetina.o <- mRetina.o[, as.numeric(str_split_fixed(rownames(metadata(mRetina.o)$revelio$dc), fixed("_"), 2)[, 2])]
mRetina.o$CCStage <- metadata(mRetina.o)$revelio$cc


ccstage2.pl <- do.call(c, lapply(list(neurosphere.o, hipp.o, endo.o, mRetina.o, HeLa1.o, HeLa2.o), function(sce.o) {
	pct.theta.p <- plotLoess(sce.o, col.name = "pct0", col.outname = "Pct0", y_lab = "Pct of 0 (projection genes)") + theme(legend.position = "none")
	den.p <- plotThetaDen(sce.o)$circular + theme(legend.position = "none")
	pct0.box.p <- plotPct0Box(sce.o) 
	top2a.box.p <- plotGeneBox(sce.o, "top2a") + theme(legend.position = "none")
	smc4.box.p <- plotGeneBox(sce.o, "smc4") + theme(legend.position = "none")
	legend.p <- get_legend(pct0.box.p)
	return(list(pct.theta.p, den.p, pct0.box.p + theme(legend.position = "none"), top2a.box.p, smc4.box.p, legend.p))
}))

mp2 <- plot_grid(plotlist  = ccstage2.pl, 
								nrow = 6, ncol = 6, rel_widths = c(1, 1, rep(1, 3), 0.3), labels = do.call(c, lapply(letters[1:6], function(x) c(x, rep("", 5)))), align = "hv", axis = "lrtb")


save_plot(here::here("figs", "sfigs", "sfig.ccstage2.pdf"), mp2,
					base_height = 2, base_width = 2*1.2, nrow = 6, ncol = 4.5, device = cairo_pdf)












