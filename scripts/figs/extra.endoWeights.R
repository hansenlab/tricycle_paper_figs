rm(list=ls())
source(here::here("scripts/utils.R"))

library(SingleCellExperiment)
library(scuttle)
library(tricycle)
library(ggrepel)

data(neuroRef)
endo.o <- qread(here::here("data/plotdata/endo.qs"))
endo_full.o <- qread(here::here("data/endo.qs"))

metadata(endo_full.o) <- metadata(endo.o)

endo_rotation.m <- attr(reducedDim(endo.o), "rotation")

int_genes.v <- intersect(rownames(endo_rotation.m), neuroRef$symbol)

### plot scatter plot
tmp.df <- data.frame(neuro_pc1 = neuroRef$pc1.rot[match(int_genes.v, neuroRef$symbol)],
										 neuro_pc2 = neuroRef$pc2.rot[match(int_genes.v, neuroRef$symbol)],
										 endo_pc1 = endo_rotation.m[int_genes.v, 1],
										 endo_pc2 = endo_rotation.m[int_genes.v, 2],
										 gene = int_genes.v)

tmp1.df <- tmp.df %>% dplyr::filter((abs(neuro_pc1) > 0.1) | (abs(endo_pc1) > 0.1))
tmp2.df <- tmp.df %>% dplyr::filter((abs(neuro_pc2) > 0.1) | (abs(endo_pc2) > 0.1))
																		
pc1.p <- ggplot(tmp.df, aes(x = neuro_pc1, y = endo_pc1, label = gene)) +
	geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point( size = 0.5,  alpha = 0.4,  shape = 1) +
	geom_text_repel(data = tmp1.df, size = 1.3, segment.size = 0.2, max.overlaps = 25) +
	labs( y = "mPancreas PC1 weights", x = "mNeurosphere PC1 weights", title = str_c("PC1 weights of overlapped genes(n=", length(int_genes.v), ")")) 

pc2.p <- ggplot(tmp.df, aes(x = neuro_pc2, y = endo_pc2, label = gene)) +
	geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point( size = 0.5,  alpha = 0.4,  shape = 1) +
	geom_text_repel(data = tmp2.df, size = 1.3, segment.size = 0.2, max.overlaps = 25) +
	labs( y = "mPancreas PC2 weights", x = "mNeurosphere PC2 weights", title = str_c("PC2 weights of overlapped genes(n=", length(int_genes.v), ")")) 


mp <- plot_grid(pc1.p, pc2.p,
								nrow = 1, ncol = 2, label_size = 10, labels = "auto", align = "hv", axis = "tblr")

save_plot(here::here("figs", "response_figs", "extra.endoWeights.pdf"), mp,
					base_height = 2, base_width = 2*1.5, nrow = 1, ncol = 2, device = cairo_pdf)



