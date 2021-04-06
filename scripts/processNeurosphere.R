rm(list=ls())
source(here::here("scripts/utils.R"))

### read in sce
filename <- "neurosphere"
sce.o <- qread(here::here("data", str_c(filename, ".qs")))


dataname <- "mNeurosphere"
species <- "mouse"
point.size <- 2.01
point.alpha <- 0.6
gene <- rowData(sce.o)$Gene
GENE <- toupper(gene)
ensembl <- rowData(sce.o)$Accession

sce.o$sample[sce.o$sample == "AX1"] <- "WT1"
sce.o$sample[sce.o$sample == "AX2"] <- "WT2"

### rerun PCA and UMAP
set.seed(1000)
corrected.m <- seuratIntegrate(count.m = assay(sce.o, "spliced"), batch = sce.o$sample)
reducedDim(sce.o, "matched.PCA.s") <- calculatePCA(corrected.m, ncomponents = 30)
set.seed(1000)
sce.o <- runUMAP(sce.o, dimred = "matched.PCA.s", external_neighbors = FALSE, pca = 30,  name = "matched.UMAP.s")

###
GO.o <- getGO(sce.o, row.id = ensembl, id.type = c("ENSEMBL"), species = species, ncomponents = 20, seed = 100,
													 runSeuratBy = "sample")
go.pca.m <- reducedDim(GO.o, "matched.PCA.s")
reducedDim(sce.o, "go.pca") <- go.pca.m

### extract genes
sce.o$top2a <- assay(sce.o, "log.s")[which(GENE == "TOP2A"), ]
sce.o$smc4 <- assay(sce.o, "log.s")[which(GENE == "SMC4"), ]

### get pct0
sce.o$pct0 <- getPct0(sce.o, GENE)

### run cyclone
sce.o$cyclone <- runCyclone(sce.o = sce.o, gname = ensembl, assay.type = "spliced", species = species, id.type = "ensembl")

### run seurat cc
sce.o$SeuratCC <- runSeuratCC(sce.o = sce.o,gname = gene, assay.type = "spliced", species = species, id.type = "name")

### get 5 stage assignment
sce.o <- estimate_cycle_stage(sce.o, gname = gene, batch.v = sce.o$sample, exprs_values = "log.s", gname.type = "SYMBOL", species = species)

### get theta and our projection
sce.o <- estimate_cycle_position(sce.o, gname = gene, gname.type = "SYMBOL", species = species, exprs_values = "log.s")

### run revelio
metadata(sce.o)$revelio <- runRevelio(sce.o, GENENAME.v = GENE, colnames.v = colnames(sce.o))

### run peco
training.l <- qread(here::here("data/peco_train.qs"))
metadata(sce.o)$peco <- runPeco(sce.o, training.l, GENE)

### data is too big for reCAT


### save big data
qsave(sce.o, file =  here::here("data", str_c(filename, ".qs")))



### make plotting data
sel.col <- c("TotalUMIs", "top2a", "smc4", "pct0", "tricyclePosition", "cyclone", "SeuratCC", "CCStage", "sample", "cell_type")
sel.col[!(sel.col %in% names(colData(sce.o)))]
colData.df <- colData(sce.o)[, sel.col[(sel.col %in% names(colData(sce.o)))]]

rowData.df <- rowData(sce.o)[, c("Accession", "Gene")]
	
metadata.l <- list(revelio = metadata(sce.o)$revelio, peco = metadata(sce.o)$peco ,
									 dataname = dataname, species = species, point.size = point.size, point.alpha = point.alpha,
									 rotation = metadata(sce.o)$rotation.m)

reducedDims.l <- list(go.pca = go.pca.m, umap = reducedDim(sce.o, "matched.UMAP.s"), tricycleEmbedding = reducedDim(sce.o, "tricycleEmbedding"))

new.o <- SingleCellExperiment(reducedDims = reducedDims.l, colData = colData.df, rowData = rowData.df ,metadata = metadata.l)


### save small data
qsave(new.o, file = here::here("data", "plotdata", str_c(filename, ".qs")))














