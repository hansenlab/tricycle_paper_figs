rm(list=ls())
source(here::here("scripts/utils.R"))

### read in sce
filename <- "mESC"
sce.o <- qread(here::here("data", str_c(filename, ".qs")))


dataname <- "mESC"
species <- "mouse"
point.size <- 4.01
point.alpha <- 0.6
gene <- rowData(sce.o)$Gene
GENE <- toupper(gene)
ensembl <- rownames(sce.o)



###
GO.o <- getGO(sce.o, row.id = ensembl, id.type = c("ENSEMBL"), species = species, ncomponents = 20, seed = 100,
							runSeuratBy = NULL)
go.pca.m <- reducedDim(GO.o, "PCA.s")
reducedDim(sce.o, "go.pca") <- go.pca.m

### extract genes
sce.o$top2a <- assay(sce.o, "log.s")[which(GENE == "TOP2A"), ]
sce.o$smc4 <- assay(sce.o, "log.s")[which(GENE == "SMC4"), ]

### get pct0
sce.o$pct0 <- getPct0(sce.o, GENE)

### run cyclone
sce.o$cyclone <- runCyclone(sce.o = sce.o, gname = ensembl, assay.type = "spliced", species = species, id.type = "ensembl")

### run seurat cc
sce.o$SeuratCC <- runSeuratCC(sce.o = sce.o, gname = gene, assay.type = "spliced", species = species, id.type = "name")

### get 5 stage assignment
sce.o <- estimate_cycle_stage(sce.o, gname = gene, batch.v = sce.o$sample, exprs_values = "log.s", gname.type = "SYMBOL", species = species)

### get theta and our projection
sce.o <- estimate_cycle_position(sce.o, gname = gene, gname.type = "SYMBOL", species = species, exprs_values = "log.s")

### run revelio
metadata(sce.o)$revelio <- runRevelio(sce.o, GENENAME.v = GENE, colnames.v =  seq_len(ncol(sce.o)))

### run peco
training.l <- qread(here::here("data/peco_train.qs"))
metadata(sce.o)$peco <- runPeco(sce.o, training.l, GENE)

### runreCAT
metadata(sce.o)$reCAT <- runreCAT(sce.o)

### save big data
qsave(sce.o, file =  here::here("data", str_c(filename, ".qs")))



### make plotting data
names(colData(sce.o))
sel.col <- c("TotalUMIs", "top2a", "smc4", "pct0", "tricyclePosition", "cyclone", "SeuratCC", "CCStage", "sample", "cell_type", "age", "celltype",  "type", "stage")
sel.col[!(sel.col %in% names(colData(sce.o)))]

colData.df <- colData(sce.o)[, sel.col[(sel.col %in% names(colData(sce.o)))]]
names(colData.df)[names(colData.df) == "type"] <- "stage"
colData.df$stage <- factor(colData.df$stage, levels = c("G1", "S", "G2M"))

names(metadata(sce.o))
metadata.l <- list(revelio = metadata(sce.o)$revelio, peco = metadata(sce.o)$peco, reCAT = metadata(sce.o)$reCAT,
									 dataname = dataname, species = species, point.size = point.size, point.alpha = point.alpha)

reducedDims.l <- list(go.pca = go.pca.m, umap = reducedDim(sce.o, "UMAP.s"), tricycleEmbedding = reducedDim(sce.o, "tricycleEmbedding"))

new.o <- SingleCellExperiment(reducedDims = reducedDims.l, colData = colData.df, metadata = metadata.l)


### save small data
qsave(new.o, file = here::here("data", "plotdata", str_c(filename, ".qs")))














