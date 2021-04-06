rm(list=ls())
source(here::here("scripts/utils.R"))

### read in sce
filename <- "hfIntestineSub"
sce.o <- qread(here::here("data", str_c(filename, ".qs")))


dataname <- "hfIntestineSub"
species <- "human"
point.size <- 3.01
point.alpha <- 0.6
gene <- rowData(sce.o)$Gene
GENE <- toupper(gene)
ensembl <- rowData(sce.o)$Accession



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
sce.o$cyclone <- runCyclone(sce.o = sce.o, gname = ensembl, assay.type = "counts", species = species, id.type = "ensembl")

### run seurat cc
colnames(sce.o) <- seq_len(ncol(sce.o))
sce.o$SeuratCC <- runSeuratCC(sce.o = sce.o, gname = gene, assay.type = "counts", species = species, id.type = "name")

### get 5 stage assignment
sce.o <- estimate_cycle_stage(sce.o, gname = gene, batch.v = NULL, exprs_values = "log.s", gname.type = "SYMBOL", species = species)

### get theta and our projection
sce.o <- estimate_cycle_position(sce.o, gname = gene, gname.type = "SYMBOL", species = species, exprs_values = "log.s")

### run revelio
metadata(sce.o)$revelio <- runRevelio(sce.o, GENENAME.v = GENE, colnames.v =  seq_len(ncol(sce.o)), assay.type = "counts")

### run peco
training.l <- qread(here::here("data/peco_train.qs"))
metadata(sce.o)$peco <- runPeco(sce.o, training.l, GENE, assay.type = "counts")

### runreCAT
metadata(sce.o)$reCAT <- runreCAT(sce.o, assay.type = "counts")

### save big data
qsave(sce.o, file =  here::here("data", str_c(filename, ".qs")))



### make plotting data
names(colData(sce.o))
sel.col <- c("TotalUMIs", "top2a", "smc4", "pct0", "tricyclePosition", "cyclone", "SeuratCC", "CCStage", "sample", "cell_type", "age", "celltype",  "type", "stage")
sel.col[!(sel.col %in% names(colData(sce.o)))]

colData.df <- colData(sce.o)[, sel.col[(sel.col %in% names(colData(sce.o)))]]

names(metadata(sce.o))
metadata.l <- list(revelio = metadata(sce.o)$revelio, peco = metadata(sce.o)$peco, reCAT = metadata(sce.o)$reCAT,
									 dataname = dataname, species = species, point.size = point.size, point.alpha = point.alpha)

reducedDims.l <- list(go.pca = go.pca.m,  tricycleEmbedding = reducedDim(sce.o, "tricycleEmbedding"))

new.o <- SingleCellExperiment(reducedDims = reducedDims.l, colData = colData.df, metadata = metadata.l)


### save small data
qsave(new.o, file = here::here("data", "plotdata", str_c(filename, ".qs")))














