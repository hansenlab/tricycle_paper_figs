rm(list=ls())
source(here::here("scripts/utils.R"))

library(ggrepel)

neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
hipp.o <- qread(here::here("data/plotdata/hipp.qs"))


### 
rotation_neurosphere.m <- metadata(neurosphere.o)$rotation
rotation_hipp.m <- metadata(hipp.o)$rotation
int_rotaGenes.v <- intersect(rownames(rotation_neurosphere.m), rownames(rotation_hipp.m))

tmp.df <- data.frame(pc1.neu = rotation_neurosphere.m[int_rotaGenes.v, 1],
										 pc2.neu = rotation_neurosphere.m[int_rotaGenes.v, 2],
										 pc1.hipp = rotation_hipp.m[int_rotaGenes.v, 1],
										 pc2.hipp = rotation_hipp.m[int_rotaGenes.v, 2])
tmp.df$Gene <- rowData(neurosphere.o)$Gene[match(int_rotaGenes.v, rownames(neurosphere.o))]
tmp.df %<>% mutate(color1 = (abs(pc1.neu) > 0.1) | (abs(pc1.hipp) > 0.1))
tmp.df %<>% mutate(color2 = (abs(pc2.neu) > 0.1) | (abs(pc2.hipp) > 0.1))
tmp1.df <- tmp.df %>% dplyr::filter((abs(pc1.neu) > 0.1) | (abs(pc1.hipp) > 0.1))
tmp2.df <- tmp.df %>% dplyr::filter((abs(pc2.neu) > 0.1) | (abs(pc2.hipp) > 0.1))


pc1.scat.p <- ggplot(tmp.df, aes(x = pc1.neu, y = pc1.hipp, label = Gene)) +
	geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_abline(slope = 1, intercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point(aes(color = color1), size = 0.5,  alpha = 0.4,  shape = 1) +
	geom_text_repel(data = tmp1.df, size = 1.3, segment.size = 0.2, max.overlaps = 25) +
	scale_color_manual(values = c("black", "red"), name = " ", guide = FALSE) +
	labs( y = "mHippNPC PC1 weights", x = "mNeurosphere PC1 weights", title = str_c("PC1 weights of overlapped genes(n=", length(int_rotaGenes.v), ")")) +
	annotate(geom = "text", x = 0, y = -.22, label = str_c("PCC == ", format(cor(tmp.df$pc1.neu, tmp.df$pc1.hipp), digits = 2)), parse = TRUE, size = 3)

pc2.scat.p <- ggplot(tmp.df, aes(x = pc2.neu, y = pc2.hipp, label = Gene)) +
	geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_abline(slope = 1, intercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point(aes(color = color1), size = 0.5,  alpha = 0.4,  shape = 1) +
	geom_text_repel(data = tmp1.df, size = 1.3, segment.size = 0.2, max.overlaps = 25) +
	scale_color_manual(values = c("black", "red"), name = " ", guide = FALSE) +
	labs( y = "mHippNPC PC2 weights", x = "mNeurosphere PC2 weights", title = str_c("PC2 weights of overlapped genes(n=", length(int_rotaGenes.v), ")")) +
	annotate(geom = "text", x = 0, y = -.22, label = str_c("PCC == ", format(cor(tmp.df$pc2.neu, tmp.df$pc2.hipp), digits = 2)), parse = TRUE, size = 3)



### respective weights scat
tmp.df <- data.frame(pc1 = rotation_neurosphere.m[, 1], pc2 = rotation_neurosphere.m[, 2], Gene = rowData(neurosphere.o)$Gene[match(rownames(rotation_neurosphere.m), rownames(neurosphere.o))])
tmp.df %<>% mutate(color = (abs(pc1) > 0.1) | (abs(pc2) > 0.1))
tmp1.df <- tmp.df %>% dplyr::filter((abs(pc1) > 0.1) | (abs(pc2) > 0.1))

neurosphere_rotation.scat.p <- ggplot(tmp.df, aes(x = pc1, y = pc2, label = Gene)) +
	geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point(aes(color = color), size = 0.5,  alpha = 0.4,  shape = 1) +
	geom_text_repel(data = tmp1.df, size = 1, segment.size = 0.25, max.overlaps = 25) +
	scale_color_manual(values = c("black", "red"), name = " ", guide = FALSE) +
	labs( y = "mNeurosphere PC2 weights", x = "mNeurosphere PC1 weights", title = str_c("mNeurosphere (n=", nrow(tmp.df), ")"))

tmp.df <- data.frame(pc1 = rotation_hipp.m[, 1], pc2 = rotation_hipp.m[, 2], Gene = rowData(hipp.o)$Gene[match(rownames(rotation_hipp.m), rownames(hipp.o))])
tmp.df %<>% mutate(color = (abs(pc1) > 0.1) | (abs(pc2) > 0.1))
tmp1.df <- tmp.df %>% dplyr::filter((abs(pc1) > 0.1) | (abs(pc2) > 0.1))

hipp_rotation.scat.p <- ggplot(tmp.df, aes(x = pc1, y = pc2, label = Gene)) +
	geom_hline(yintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_vline(xintercept = 0, size = 0.5, color = "grey", alpha = 0.8, linetype = "dashed") +
	geom_point(aes(color = color), size = 0.5,  alpha = 0.4, shape = 1) +
	geom_text_repel(data = tmp1.df, size = 1, segment.size = 0.25, max.overlaps = 25) +
	scale_color_manual(values = c("black", "red"), name = " ", guide = FALSE) +
	labs( y = "mHippNPC PC2 weights", x = "mHippNPC PC1 weights", title = str_c("mHippNPC (n=", nrow(tmp.df), ")"))



mp <- plot_grid(neurosphere_rotation.scat.p, hipp_rotation.scat.p, pc1.scat.p, pc2.scat.p,
							 nrow = 2, ncol = 2, label_size = 10, labels = "auto", align = "hv", axis = "tblr")

save_plot(here::here("figs", "sfigs", "sfig.weightsScat.pdf"), mp,
					base_height = 2, base_width = 2*1.5, nrow = 2, ncol = 2, device = cairo_pdf)
