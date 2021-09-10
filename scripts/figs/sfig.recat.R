rm(list=ls())
source(here::here("scripts/utils.R"))


endo.o <- qread(here::here("data/plotdata/endo.qs"))
mHSC.o <- qread(here::here("data/plotdata/mHSC.qs"))
HeLa1.o <- qread(here::here("data/plotdata/HeLa1.qs"))
HeLa2.o <- qread(here::here("data/plotdata/HeLa2.qs"))
hfIntestineSub.o <- qread(here::here("data/plotdata/hfIntestineSub.qs"))

hU2OS.o <- qread(here::here("data/plotdata/hU2OS.qs"))
hiPSCs.o <- qread(here::here("data/plotdata/hiPSCs.qs"))
mESC.o <- qread(here::here("data/plotdata/mESC.qs"))
hESC.o <- qread(here::here("data/plotdata/hESC.qs"))

sapply(list(endo.o, mHSC.o, HeLa1.o, HeLa2.o, hfIntestineSub.o, hU2OS.o, hiPSCs.o, mESC.o, hESC.o), function(x) !is.null(metadata(x)$reCAT))


plot_bayes <- function (sce.o, title = NULL) 
{
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " reCAT (n=", ncol(sce.o), ")")
	
	bayes_score <- metadata(sce.o)$reCAT$scores$bayes_score
	ordIndex <- metadata(sce.o)$reCAT$ordInde
	cell = bayes_score[ordIndex, ]
	data <- data.frame(x = rep(1:length(ordIndex), 3), 
										 y = c(cell$G1.score, cell$S.score, cell$G2M.score), 
										 z = factor(c(rep("G1 score", length(ordIndex)), 
										 						 rep("S score", length(ordIndex)), 
										 						 rep("G2/M score", length(ordIndex))), 
										 					 levels = c("G1 score", "S score", "G2/M score")))
	p <- ggplot(data)+  
		geom_scattermore(aes(x = x, y = y, color = z), pointsize = point.size, alpha = point.alpha) + 
		# geom_line(aes(x = x, y = y, color = z, group = z), size = 0.1, alpha = alpha) +
		scale_color_manual(name = "", values = c(`G1 score` = "#6699FF", `S score` = "#66CC00", `G2/M score` = "#FF3366"), 
											 breaks = c("G1 score", "S score", "G2/M score"), 
											 guide = guide_legend(keywidth = 1, keyheight = 0.6)) + 
		labs(title = title, x = "Time Series (t)", y = "Bayes-scores")
	return(p)
}

plot_mean <- function (sce.o, title= NULL) 
{
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " reCAT (n=", ncol(sce.o), ")")
	
	mean_score <- metadata(sce.o)$reCAT$scores$mean_score
	ordIndex <- metadata(sce.o)$reCAT$ordInde
	
	cell = mean_score[ordIndex, ]
	data2 <- data.frame(x = rep(1:length(ordIndex), 6), 
											y = c(cell$G1Score, cell$G1SScore, cell$SScore, cell$G2Score, cell$G2MScore, cell$MScore), 
											z = factor(c(rep("G1 score", length(ordIndex)), 
																	 rep("G1/S score", length(ordIndex)), 
																	 rep("S score", length(ordIndex)), 
																	 rep("G2 score", length(ordIndex)), 
																	 rep("G2/M score", length(ordIndex)), 
																	 rep("M score", length(ordIndex))), 
																 levels = c("G1 score", "G1/S score", "S score", "G2 score", "G2/M score", "M score")))
	
	p <- ggplot(data2) + 
		geom_scattermore(aes(x = x, y = y, color = z), pointsize = point.size, alpha = point.alpha) + 
		#geom_line(aes(x = x, y = y, color = z, group = z), size = 0.1, alpha = alpha) +
		scale_color_manual(name = "", values = c(`G1 score` = "#6699FF", 
																						 `G1/S score` = "black", `S score` = "#66CC00", `G2 score` = "orange", 
																						 `G2/M score` = "#FF3366", `M score` = "yellow"), 
											 breaks = c("G1 score", "G1/S score", "S score", "G2 score", 
											 					 "G2/M score", "M score"), guide = guide_legend(keywidth = 1, 
											 					 																							 keyheight = 0.6)) + 
		labs(title = title, x = "Time Series (t)", y = "Mean-scores")
	
	return(p)
}

plotGenereCAT <- function(sce.o, col.name = "top2a", col.outname = NULL, title = NULL, color.var = "CCStage", color.name = "SchwabeCC", colors.v = NULL, labels.v = NULL) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	
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
	if (is.null(colors.v)) colors.v <- ccColors.v
	if (is.null(labels.v)) labels.v <- ccLabels.v
	
	ordIndex <- metadata(sce.o)$reCAT$ordInde
	gene.v <- colData(sce.o)[[col.name]]
	tmp.df <- data.frame(x = seq_along(ordIndex), y =  gene.v[ordIndex], color = colData(sce.o)[, color.var][ordIndex])
	if (any(is.na(tmp.df$color))) {
		tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
		scale_color <- scale_color_manual(values = c(colors.v, "grey"), name = color.name, labels =  c(labels.v, "NA"), limits =   c(labels.v, "NA"))
	} else {
		scale_color <- scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = labels.v)
	}
	p <- ggplot(tmp.df, aes(x = x, y = y)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size, alpha = point.alpha, color = "grey") +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size, alpha = point.alpha, color = "grey") +
		geom_smooth(aes(x = x, y = y), inherit.aes = FALSE, method = "loess", se = FALSE, color = "black", linetype = "dashed", size = 0.8, alpha = 0.8) +
		labs( x = "reCAT Time Series (t)", y = bquote(paste('log'['2'],'(expression of ', .(col.outname), ")")), title = title) 
	
	return(p)
}

plotreCATTheta <- function(sce.o, title = NULL, color.var = "CCStage", color.name = "SchwabeCC", colors.v = NULL, labels.v = NULL) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha

	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " \u03B8 (n=", ncol(sce.o), ")")
	if (is.null(colors.v)) colors.v <- ccColors.v
	if (is.null(labels.v)) labels.v <- ccLabels.v
	
	theta.v <- sce.o$tricyclePosition
	ordIndex <- metadata(sce.o)$reCAT$ordInde
	tmp.df <- data.frame(x = seq_along(ordIndex), y =  theta.v[ordIndex], color = colData(sce.o)[, color.var][ordIndex])
	if (any(is.na(tmp.df$color))) {
		tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
		scale_color <- scale_color_manual(values = c(colors.v, "grey"), name = color.name, labels =  c(labels.v, "NA"), limits =   c(labels.v, "NA"))
	} else {
		scale_color <- scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = labels.v)
	}
	
	p <- ggplot(data = tmp.df , aes(x = x, y = y )) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size, alpha = point.alpha, color = "grey") +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size, alpha = point.alpha, color = "grey") +
		labs( x = "reCAT Time Series (t)", y = theta_lab, title = title) +
		scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2)* pi) 
	return(p)
}


recat.lp <- do.call(c, lapply(list(endo.o, mHSC.o, HeLa1.o, HeLa2.o, hfIntestineSub.o, hU2OS.o, hiPSCs.o, mESC.o, hESC.o), function(sce.o) {
	bayes.p <- plot_bayes(sce.o) + theme(legend.position = "none")
	mean.p <- plot_mean(sce.o) + theme(legend.position = "none")
	top2a.p <- plotGenereCAT(sce.o) + theme(legend.position = "none")
	theta.p <- plotreCATTheta(sce.o) + theme(legend.position = "none")
	list(bayes.p, mean.p, top2a.p, theta.p) 
}))
# recat2.lp <- do.call(c, lapply(list(mESC.o, hESC.o), function(sce.o) {
# 	bayes.p <- plot_bayes(sce.o) + theme(legend.position = "none")
# 	mean.p <- plot_mean(sce.o) + theme(legend.position = "none")
# 	top2a.p <- plotGenereCAT(sce.o, color.name = "FACS", color.var = "stage", colors.v = cc3Colors.v, labels.v = cc3Labels.v) + theme(legend.position = "none")
# 	theta.p <- plotreCATTheta(sce.o, color.name = "FACS", color.var = "stage", colors.v = cc3Colors.v, labels.v = cc3Labels.v) + theme(legend.position = "none")
# 	list(bayes.p, mean.p, top2a.p, theta.p) 
# }))

# recat.lp <- c(recat1.lp, recat2.lp)

recat.lp[[1]] <- recat.lp[[1]] + theme(legend.position = c(0, 0),
																								 legend.justification = c(0, 0),
																								 legend.key = element_blank(),
																								 legend.key.size = unit(6, "pt"))
recat.lp[[2]] <- recat.lp[[2]] + theme(legend.position = c(0, 1.05),
																								 legend.justification = c(0, 1),
																								 legend.key = element_blank(),
																								 legend.key.size = unit(6, "pt"))
recat.lp[[3]] <- recat.lp[[3]] + theme(legend.position = c(0, 1),
																			 legend.justification = c(0, 1),
																			 legend.key = element_blank(),
																			 legend.key.size = unit(6, "pt"))
recat.lp[[31]] <- recat.lp[[31]] + theme(legend.position = c(1, 0),
																			 legend.justification = c(1, 0),
																			 legend.key = element_blank(),
																			 legend.key.size = unit(6, "pt"))

mp <- plot_grid(plotlist = recat.lp,
								nrow = 5, ncol = 8,  labels = as.vector(rbind(letters[1:9], rep("", 9), rep("", 9), rep("", 9))), align = "vh", axis = "tblr")
save_plot(here::here("figs", "sfigs", "sfig.recat.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 5, ncol = 8, device = cairo_pdf)






