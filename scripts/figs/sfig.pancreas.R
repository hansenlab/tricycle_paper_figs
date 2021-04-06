rm(list=ls())
source(here::here("scripts/utils.R"))

endo.o <- qread(here::here("data/plotdata/endo.qs"))
endo_full.o <- qread(here::here("data/endo.qs"))

metadata(endo.o)$point.size <- 4.01
metadata(endo_full.o) <- metadata(endo.o)


### ductal
GO_Ductal.o <- getGO(endo_full.o[, endo_full.o$cell_type == "Ductal"], row.id = rownames(endo_full.o), id.type = "SYMBOL", runSeuratBy =  NULL)
plot(reducedDim(GO_Ductal.o, "PCA.s"))

ref_Ductal.m <- attr(reducedDim(GO_Ductal.o, "PCA.s"), "rotation")

GO_Ngn3LEp.o <- getGO(endo_full.o[,endo_full.o$cell_type == "Ngn3 low EP"], row.id = rownames(endo_full.o), id.type = "SYMBOL", runSeuratBy =  NULL, seed = 400)
plot(reducedDim(GO_Ngn3LEp.o, "PCA.s"))
GO_Ngn3HEp.o <- getGO(endo_full.o[, endo_full.o$cell_type == "Ngn3 high EP"], row.id = rownames(endo_full.o), id.type = "SYMBOL", runSeuratBy =  NULL)
plot(reducedDim(GO_Ngn3HEp.o, "PCA.s"))

GO_Ngn3LEp.o <- project_cycle_space(GO_Ngn3LEp.o, ref.m = ref_Ductal.m, exprs_values = "log.s", name = "ductal_proj")
GO_Ngn3HEp.o <- project_cycle_space(GO_Ngn3HEp.o, ref.m = ref_Ductal.m, exprs_values = "log.s", name = "ductal_proj")

plot(reducedDim(GO_Ngn3LEp.o, "ductal_proj"))
plot(reducedDim(GO_Ngn3HEp.o, "ductal_proj"))




ductal.pca.p <- plotScatCC(sce.o = GO_Ductal.o, dimred = "PCA.s", 
													 x_lab = str_c("PC1", sprintf(" (%i%%)", round(attr(reducedDim(GO_Ductal.o, "PCA.s"), "percentVar")[1]))),
													 y_lab = str_c("PC2", sprintf(" (%i%%)", round(attr(reducedDim(GO_Ductal.o, "PCA.s"), "percentVar")[2]))),
													 title = str_c("PCA of Ductal cells (n=", ncol(GO_Ductal.o), ")"))
Ngn3LEp.pca.p <- plotScatCC(sce.o = GO_Ngn3LEp.o, dimred = "PCA.s", 
													 x_lab = str_c("PC1", sprintf(" (%i%%)", round(attr(reducedDim(GO_Ngn3LEp.o, "PCA.s"), "percentVar")[1]))),
													 y_lab = str_c("PC2", sprintf(" (%i%%)", round(attr(reducedDim(GO_Ngn3LEp.o, "PCA.s"), "percentVar")[2]))),
													 title = str_c("PCA of Ngn3LEp cells (n=", ncol(GO_Ngn3LEp.o), ")"))
Ngn3HEp.pca.p <- plotScatCC(sce.o = GO_Ngn3HEp.o, dimred = "PCA.s", 
														x_lab = str_c("PC1", sprintf(" (%i%%)", round(attr(reducedDim(GO_Ngn3HEp.o, "PCA.s"), "percentVar")[1]))),
														y_lab = str_c("PC2", sprintf(" (%i%%)", round(attr(reducedDim(GO_Ngn3HEp.o, "PCA.s"), "percentVar")[2]))),
														title = str_c("PCA of Ngn3HEp cells(n=", ncol(GO_Ngn3HEp.o), ")"))
Ngn3LEp.proj.p <- plotScatCC(sce.o = GO_Ngn3LEp.o, dimred = "ductal_proj", 
														x_lab = bquote(paste('CC'['ductal']," Space Dim 1")),
														y_lab = bquote(paste('CC'['ductal']," Space Dim 2")),
														title = str_c("Ngn3LEp projected by Ductal reference"))
Ngn3HEp.proj.p <- plotScatCC(sce.o = GO_Ngn3HEp.o, dimred = "ductal_proj", 
														x_lab = bquote(paste('CC'['ductal']," Space Dim 1")),
														y_lab = bquote(paste('CC'['ductal']," Space Dim 2")),
														title = str_c("Ngn3HEp projection by Ductal reference"))


mp <- plot_grid(ductal.pca.p + theme(legend.position = "none"), 
								Ngn3LEp.pca.p + theme(legend.position = "none"), 
								Ngn3HEp.pca.p + theme(legend.position = "none"), 
								get_legend(ductal.pca.p),
								Ngn3LEp.proj.p + theme(legend.position = "none"), 
								Ngn3HEp.proj.p + theme(legend.position = "none"), 

								nrow = 2, ncol = 3,  label_size = 10, labels = c("a", "b", "c", " ", "d", "e"), align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.pancreas.pdf"), mp,
					base_height = 2, base_width = 2*1.2, nrow = 2, ncol = 3, device = cairo_pdf)



