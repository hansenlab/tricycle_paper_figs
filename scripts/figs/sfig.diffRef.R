rm(list=ls())
source(here::here("scripts/utils.R"))


neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))
mRetina.o <- qread(here::here("data/plotdata/mRetina.qs"))

neurosphere_full.o <- qread(here::here("data/neurosphere.qs"))
hipp_full.o <- qread(here::here("data/hipp.qs"))
mRetina_full.o <- qread(here::here("data/mRetina.qs"))


rotation_hipp.m <- metadata(hipp.o)$rotation
rotation_hipp_acc.m <- rotation_hipp.m
rownames(rotation_hipp_acc.m) <- rowData(hipp.o)$Accession[match(rownames(rotation_hipp.m), rownames(hipp.o))]


reducedDim(neurosphere.o, "hippP") <- project_cycle_space(assay(neurosphere_full.o, "log.s"), ref.m = rotation_hipp.m)
reducedDim(hipp.o, "hippP") <- project_cycle_space(assay(hipp_full.o, "log.s"), ref.m = rotation_hipp.m)
reducedDim(mRetina.o, "hippP") <- project_cycle_space(assay(mRetina_full.o, "log.s"), ref.m = rotation_hipp_acc.m)

neurosphere.o$tricyclePosition2 <- as.numeric(coord2rad(reducedDim(neurosphere.o, "hippP")[, 1:2]))
hipp.o$tricyclePosition2 <- as.numeric(coord2rad(reducedDim(hipp.o, "hippP")[, 1:2]))
mRetina.o$tricyclePosition2 <- as.numeric(coord2rad(reducedDim(mRetina.o, "hippP")[, 1:2]))


px_lab2 <-  bquote(paste('CC'['hipp']," Space Dim 1"))
py_lab2 <-  bquote(paste('CC'['hipp']," Space Dim 2"))
theta_lab2 <- bquote(paste('CC'['hipp']," Position \u03B8"))


neurosphere.proj.p <- plotScatCC(neurosphere.o, dimred = "hippP", x_lab = px_lab2, y_lab = py_lab2, title = str_c(metadata(neurosphere.o)$dataname, " projected by mHippNPC"))
hipp.proj.p <- plotScatCC(hipp.o, dimred = "hippP", x_lab = px_lab2, y_lab = py_lab2, title = str_c(metadata(hipp.o)$dataname, " projected by mHippNPC"))
mRetina.proj.p <- plotScatCC(mRetina.o, dimred = "hippP", x_lab = px_lab2, y_lab = py_lab2, title = str_c(metadata(mRetina.o)$dataname, " projected by mHippNPC"))


plotTwoTheta <- function(sce.o, x_lab, y_lab, x_var, y_var, color.var = "CCStage", color.name = "CC Stage", colors.v = NULL, labels.v = NULL, title = NULL) {
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " (n=", ncol(sce.o), ")")
	if (is.null(colors.v)) colors.v <- ccColors.v
	if (is.null(labels.v)) labels.v <- ccLabels.v
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	
	tmp.df <- data.frame(x = colData(sce.o)[, x_var], y = colData(sce.o)[, y_var], color = colData(sce.o)[, color.var])
	
	if (any(is.na(tmp.df$color))) {
		tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
		scale_color <- scale_color_manual(values = c(colors.v, "grey"), name = color.name, labels =  c(labels.v, "NA"), limits =   c(labels.v, "NA"))
	} else {
		scale_color <- scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = labels.v)
	}
	
	p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size,  alpha = point.alpha) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size,  alpha = point.alpha) +
		scale_color +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( y = y_lab, x = x_lab, title = title) +
		annotate(geom = "text", x = 2, y = 5, size = 3, hjust = 0.5, vjust = 0.5,
						 label = str_c("Circular correlation \n \u03C1=", sprintf("%.3f", Directional::circ.cor1(tmp.df$x, tmp.df$y , rads = T)[1])), parse = FALSE) +
		scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi , labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2 )* pi) +
		scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi , labels = str_c(seq(0, 2, 0.5), "\u03C0"), limits = c(0, 2)* pi)
	p
}

neurosphere.theta.p <- plotTwoTheta(neurosphere.o, x_lab = theta_lab, y_lab = theta_lab2, x_var = "tricyclePosition", y_var = "tricyclePosition2")
hipp.theta.p <- plotTwoTheta(hipp.o, x_lab = theta_lab, y_lab = theta_lab2, x_var = "tricyclePosition", y_var = "tricyclePosition2")
mRetina.theta.p <- plotTwoTheta(mRetina.o, x_lab =theta_lab, y_lab = theta_lab2, x_var = "tricyclePosition", y_var = "tricyclePosition2")


mp <- plot_grid(neurosphere.proj.p + theme(legend.position = c(0, 1),
																					 legend.justification = c(0, 1),
																					 legend.key = element_blank(),
																					 legend.key.size = unit(7, "pt")), 
								hipp.proj.p + theme(legend.position = "none"),
								mRetina.proj.p + theme(legend.position = "none"),
								neurosphere.theta.p + theme(legend.position = "none"),
								hipp.theta.p + theme(legend.position = "none"),
								mRetina.theta.p + theme(legend.position = "none"), 
								nrow = 2, ncol = 3, labels = c("auto"), align = "hv", axis = "lrtb")

save_plot(here::here("figs", "sfigs", "sfig.diffRef.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 3, device = cairo_pdf)











