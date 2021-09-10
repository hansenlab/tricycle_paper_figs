rm(list=ls())
source(here::here("scripts/utils.R"))

library(ggrepel)

neurosphere.o <- qread(here::here("data/neurosphere.qs"))
hipp.o <- qread(here::here("data/hipp.qs"))
endo.o <- qread(here::here("data/endo.qs"))



### get cell cycle genes
cyclebase.df <- read_tsv(here::here("data/cyclebase3.0_genes.tsv"))


require(org.Hs.eg.db)
cyclebase_genes.v <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = cyclebase.df$gene, column = "SYMBOL", keytype = "ENSEMBLPROT", multiVals = "first")

go_genes.v <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", 
																		columns = "SYMBOL")[, "SYMBOL"]
diff_genes.v <- na.omit(setdiff(cyclebase_genes.v, cycle.anno))

require(org.Mm.eg.db)
diff_ensemle.v <- na.omit(unlist(AnnotationDbi::mapIds(org.Mm.eg.db, keys =  str_to_title(diff_genes.v),  keytype= "SYMBOL", column = "ENSEMBL", multiVals = "list")))
	
go_ensembl_mouse.v <- 	AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL", keys="GO:0007049", 
																						 columns = "ENSEMBL")[, "ENSEMBL"]
	
all_ensembl_mouse.v <- union(go_ensembl_mouse.v, diff_ensemle.v)
	
go_genes_mouse.v <- 	AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL", keys="GO:0007049", 
																						 columns = "SYMBOL")[, "SYMBOL"]
all_genes_mouse.v <- union(go_genes_mouse.v, str_to_title(diff_genes.v))


### get R2
neurosphere_logs.m <- assay(neurosphere.o, "log.s")
neurosphere_r2.v <- unlist(bplapply(seq_len(nrow(neurosphere.o)), function(i) {
	fit_periodic_loess(theta.v = neurosphere.o$tricyclePosition, y = neurosphere_logs.m[i, ])$rsquared
}, BPPARAM = MulticoreParam(workers = 10)))
	
type_neurosphere.v <- rowData(neurosphere.o)$Accession %in% all_ensembl_mouse.v

tmp.df <- data.frame(r2 = neurosphere_r2.v, type = type_neurosphere.v, gene = rowData(neurosphere.o)$Gene)

neurosphere.p <- ggplot(tmp.df, aes(x = type, y = r2)) +
	geom_boxplot(outlier.shape = 1, outlier.size = 0.5, size = 0.2, width = 0.8) +
	scale_x_discrete(name = "Whether in cell cycle gene list", labels = str_c(c("False", "True"), "\n(n=", table(tmp.df$type), ")")) +
	labs( y =  bquote(paste( "R"^"2")),  title = str_c( "mNeurosphere (n=", ncol(neurosphere.o), ")")) +
	geom_text(data = tmp.df %>% dplyr::filter(! `type`) %>% arrange(desc(`r2`)) %>% top_n(3, `r2`),
						aes(x = 1.1, y = r2, label = gene), inherit.aes = FALSE, hjust = 0, size = 3)
	
	

	
	

## hipp
hipp_logs.m <- assay(hipp.o, "log.s")
hipp_r2.v <- unlist(bplapply(seq_len(nrow(hipp.o)), function(i) {
	fit_periodic_loess(theta.v = hipp.o$tricyclePosition, y = hipp_logs.m[i, ])$rsquared
}, BPPARAM = MulticoreParam(workers = 10)))

type_hipp.v <- rowData(hipp.o)$Accession %in% all_ensembl_mouse.v

tmp.df <- data.frame(r2 = hipp_r2.v, type = type_hipp.v, gene = rowData(hipp.o)$Gene)

hipp.p <- ggplot(tmp.df, aes(x = type, y = r2)) +
	geom_boxplot(outlier.shape = 1, outlier.size = 0.5, size = 0.2, width = 0.8) +
	scale_x_discrete(name = "Whether in cell cycle gene list", labels = str_c(c("False", "True"), "\n(n=", table(tmp.df$type), ")")) +
	labs( y =  bquote(paste( "R"^"2")),  title = str_c( "mHippNPC (n=", ncol(hipp.o), ")")) +
	geom_text(data = (tmp.df %>% dplyr::filter(! `type`) %>% arrange(desc(`r2`)) %>% top_n(3, `r2`))[1, ],
						aes(x = 1.1, y = r2 + 0.05, label = gene), inherit.aes = FALSE, hjust = 0, size = 2.8) +
	geom_text(data = (tmp.df %>% dplyr::filter(! `type`) %>% arrange(desc(`r2`)) %>% top_n(3, `r2`))[3, ],
						aes(x = 1.1, y = r2 - 0.01, label = gene), inherit.aes = FALSE, hjust = 0, size = 2.8) +
	geom_text(data = (tmp.df %>% dplyr::filter(! `type`) %>% arrange(desc(`r2`)) %>% top_n(3, `r2`))[2, ],
						aes(x = 0.9, y = r2, label = gene), inherit.aes = FALSE, hjust = 1, size = 2.8)


## endo.o
endo_logs.m <- assay(endo.o, "log.s")
endo_r2.v <- unlist(bplapply(seq_len(nrow(endo.o)), function(i) {
	fit_periodic_loess(theta.v = endo.o$tricyclePosition, y = endo_logs.m[i, ])$rsquared
}, BPPARAM = MulticoreParam(workers = 10)))

type_endo.v <- rownames(endo.o) %in% all_genes_mouse.v

tmp.df <- data.frame(r2 = endo_r2.v, type = type_endo.v, gene = rownames(endo.o))

endo.p <- ggplot(tmp.df, aes(x = type, y = r2)) +
	geom_boxplot(outlier.shape = 1, outlier.size = 0.5, size = 0.2, width = 0.8) +
	scale_x_discrete(name = "Whether in cell cycle gene list", labels = str_c(c("False", "True"), "\n(n=", table(tmp.df$type), ")")) +
	labs( y =  bquote(paste( "R"^"2")),  title = str_c( "mPancreas (n=", ncol(endo.o), ")")) +
	geom_text(data = tmp.df %>% dplyr::filter(! `type`) %>% arrange(desc(`r2`)) %>% top_n(3, `r2`),
						aes(x = 1.1, y = r2, label = gene), inherit.aes = FALSE, hjust = 0, size = 3)


r2.l <- list(neurosphere_r2.v, hipp_r2.v, endo_r2.v)
qsave(r2.l, here::here("data/plotdata/r2.qs"))

mp <- plot_grid(
	neurosphere.p, hipp.p, endo.p,
	nrow = 1, ncol = 3, label_size = 10, labels = c("a", "b", "c"), align = "hv", axis = "tblr")

save_plot(here::here("figs", "response_figs", "response.cyclegener2.pdf"), mp,
					base_height = 2, base_width = 2 * 1.2, nrow = 1, ncol = 3, device = cairo_pdf)



