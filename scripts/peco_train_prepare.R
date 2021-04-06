rm(list=ls())
library(peco)
library(org.Hs.eg.db)

data(training_human)

GENENAME.v <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(training_human$predict.yy), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
idx <- which(!is.na(GENENAME.v))
training_human$predict.yy <- training_human$predict.yy[idx, ]
training_human$cellcycle_function <- training_human$cellcycle_function[idx]
training_human$sigma <- training_human$sigma[idx, , drop = FALSE]
training_human$pve <- training_human$pve[idx, , drop = FALSE]

GENENAME.v <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(training_human$predict.yy), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

rownames(training_human$predict.yy) <- GENENAME.v
names(training_human$cellcycle_function) <- GENENAME.v
rownames(training_human$sigma) <- GENENAME.v
rownames(training_human$pve) <- GENENAME.v

qsave(training_human, file = here::here("data/peco_train.qs"))
