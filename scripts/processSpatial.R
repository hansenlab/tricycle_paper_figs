rm(list=ls())
source(here::here("scripts/utils.R"))


library(scater)
library(scran)
library(SingleCellExperiment)

raw <- Matrix::readMM(here::here("data", "spatial/filtered_feature_bc_matrix", "matrix.mtx.gz")) %>% as("dgCMatrix")
colData <- readr::read_tsv(here::here("data", "spatial/filtered_feature_bc_matrix", "barcodes.tsv.gz"), col_names = F)
rowData <- readr::read_tsv(here::here("data", "spatial/filtered_feature_bc_matrix", "features.tsv.gz"), col_names = F)
names(colData) <- "barcode"
names(rowData) <- c("Accession", "Gene", "type")


## read in spatial coords
spatial <-  readr::read_csv(here::here("data", "spatial/spatial", "tissue_positions_list.csv"), col_names = F)
names(spatial) <- c("barcode", "under", "out", "position", "spatialx", "spatialy")

colData <- cbind(colData, spatial[match(colData$barcode, spatial$barcode), 2:6])

sce.o <- SingleCellExperiment(assays = list(spliced = raw), rowData = rowData[, 1:2], colData = colData)

colnames(sce.o) <- sce.o$barcode
rownames(sce.o) <- rowData(sce.o)$Accession



mito_genes.idx <- which(grepl("^MT-", rowData(sce.o)$Gene))
sce.o$TotalUMIs <- colSums(assay(sce.o, 'spliced'))
sce.o$n_gene <- colSums(assay(sce.o, 'spliced') > 0)
sce.o$percent_mito <- colSums(assay(sce.o[mito_genes.idx, ], 'spliced'))/sce.o$TotalUMIs
sce.o$DoubletScore <- log10(scran::doubletCells(assay(sce.o, 'spliced')) + 1)

outliers.df <- tibble(TotalUMIs = isOutlier(sce.o$TotalUMIs, nmads = 3, type = "lower", log = TRUE, batch = sce.o$new_CellType),
											n_gene_s = isOutlier(sce.o$n_gene, nmads = 3, type = "lower", log = TRUE, batch = sce.o$new_CellType),
											percent_mito = isOutlier(sce.o$percent_mito, nmads = 3, type = "higher", log = FALSE, batch = sce.o$new_CellType), 
											DoubletScore = isOutlier(sce.o$DoubletScore, nmads = 3, type = "higher", log = FALSE, batch = sce.o$new_CellType)
) %>% add_column(discard = rowAnys(as.matrix(.)))
metadata(sce.o) <- list(outliers = outliers.df)

sce.o <- sce.o[, !outliers.df$discard]

sce.o <- sce.o[!grepl("^MT-", rowData(sce.o)$Gene), ]



# remove lowly expressed genes in each sample
sce.o <- sce.o[which(rowSums(assay(sce.o, 'spliced')) > 20), ]


assay(sce.o, "log.s") <-  normalizeCounts(assay(sce.o, 'spliced'), log= TRUE)
assay(sce.o, 'spliced') <- normalizeCounts(assay(sce.o, 'spliced'), log=FALSE)


# get HVG 
hvg.s.df <- modelGeneVar(sce.o, assay.type = "log.s", BPPARAM = MulticoreParam(workers = 12L))
hvg.s.v <- getTopHVGs(hvg.s.df, n = 2000)
str(hvg.s.v)

# Dimensionality reduction.
set.seed(1000)
sce.o <- runPCA(sce.o, exprs_values = "log.s", subset_row = hvg.s.v, ncomponents = 30, name = "PCA.s")
sce.o <- runUMAP(sce.o, dimred = 'PCA.s', exprs_values = "log.s", external_neighbors = FALSE, pca = 30, name = "UMAP.s")



filename <- dataname <- "spatial"
species <- "human"
point.size <- 3.01
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
#sce.o$cyclone <- runCyclone(sce.o = sce.o, gname = gene, assay.type = "spliced", species = species, id.type = "name")

### run seurat cc
sce.o$SeuratCC <- runSeuratCC(sce.o = sce.o, gname = gene, assay.type = "spliced", species = species, id.type = "name")

### get 5 stage assignment
sce.o <- estimate_cycle_stage(sce.o, gname = gene, batch.v = rep(1, ncol(sce.o)), exprs_values = "log.s", gname.type = "SYMBOL", species = species)

### get theta and our projection
sce.o <- estimate_cycle_position(sce.o, gname = gene, gname.type = "SYMBOL", species = species, exprs_values = "log.s")

### run revelio
metadata(sce.o)$revelio <- runRevelio(sce.o, GENENAME.v = GENE, colnames.v =  str_c(sce.o$age, "_", seq_len(ncol(sce.o))))


metadata(sce.o) <- c(metadata(sce.o), list(dataname = dataname, species = species, point.size = point.size, point.alpha = point.alpha))


### save big data
qsave(sce.o, file =  here::here("data", str_c(filename, ".qs")))






### make plotting data
new.o <- sce.o
assays(new.o) <- list()


### save small data
qsave(new.o, file = here::here("data", "plotdata", str_c(filename, ".qs")))





