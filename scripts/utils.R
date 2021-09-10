options(stringsAsFactors = F)

library("magrittr")
library("matrixStats")
library("Matrix")
library("RColorBrewer")
library("scales")
library("tidyverse")
library("cowplot")
library("ggbeeswarm")
library("viridis")
library("colorspace")
library("zeallot")
library("sparseMatrixStats")
library("qs")
library("ggridges")
library("scattermore")
library("SingleCellExperiment")
library("scater")
library("BiocParallel")
library("circular")
library("here")
library("tricycle")



### set ggplot theme
.new_theme <- theme_bw(base_size = 8) + 
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_line(size = 0.25),
		#axis.line.x = element_blank(),
		#axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
		axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
		axis.text.y = element_text(size = 6),
		axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5),
		axis.title.y = element_text(size = 7),
		axis.title.x = element_text(size = 7),
		plot.title = element_text(face = "plain", size = 8, hjust = 0.5),
		plot.margin = unit(c(5, 5, 6, 6), "pt"),
		legend.background = element_blank(),
		legend.position = "right",
		legend.justification = c(0.5, 0.5),
		legend.key.size = unit(7, "pt"),
		legend.key = element_blank(),
		legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
		legend.title = element_text(size = 5),
		legend.title.align = 0.5,
		legend.text.align = 0.5,
		legend.text = element_text(size = 5),
		strip.text = element_text(color = "white"),
		strip.background = element_rect(fill = "black", color = "black"))

theme_set(.new_theme)

theme_main <- theme_bw(base_size = 8) + 
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_line(size = 0.25),
		#axis.line.x = element_blank(),
		#axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
		axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
		axis.text.y = element_text(size = 6),
		axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5),
		axis.title.y = element_text(size = 8),
		axis.title.x = element_text(size = 8),
		plot.title = element_text(face = "plain", size = 9, hjust = 0.5),
		plot.margin = unit(c(5, 5, 6, 6), "pt"),
		legend.background = element_blank(),
		legend.position = "right",
		legend.justification = c(0.5, 0.5),
		legend.key.size = unit(7, "pt"),
		legend.key = element_blank(),
		legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
		legend.title = element_text(size = 5),
		legend.title.align = 0.5,
		legend.text.align = 0.5,
		legend.text = element_text(size = 5),
		strip.text = element_text(color = "white"),
		strip.background = element_rect(fill = "black", color = "black"))



### help function for annotate
.percent_range <- function(value.v, percent=0.9) {
	value.v <- na.omit(value.v)
	min(value.v) + percent*diff(range(value.v))
}


### set cc colors and labels
ccColors.v <- c('#B2627C', '#F29360', '#FCEA64', '#86BBD8', '#8159ba')
ccLabels.v <- c("G1.S", "S","G2", "G2.M", "M.G1")

cc3Colors.v <- c("#82A93F", '#ac4343', "#1F83B4")
cc3Labels.v <- c("G1", "S", "G2M")

px_lab <-  bquote(paste('CC'['ns']," Space Dim 1"))
py_lab <-  bquote(paste('CC'['ns']," Space Dim 2"))

theta_lab <- bquote(paste('CC'['ns']," Position \u03B8"))

### getPct0
getPct0 <- function(sce.o, GENE) {
ref.df <- tricycle::neuroRef
data.m <- assay(sce.o, "log.s")[GENE %in% ref.df$SYMBOL, ] 
pct0.v <- colMeans(data.m < 0.01)
return(pct0.v)
}

### seurat integration
seuratIntegrate <- function(count.m, batch, nfeatures = 2000) {
	require(Seurat)
	seurat.o <- CreateSeuratObject(counts = count.m)
	seurat.o[["batch"]] <- batch
	
	seurat.list <- SplitObject(seurat.o, split.by = "batch")
	
	for (i in 1:length(seurat.list)) {
		seurat.list[[i]] <- NormalizeData(seurat.list[[i]], verbose = TRUE)
		seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]], selection.method = "vst", 
																						 nfeatures = nfeatures, verbose = TRUE)
	}
	
	seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:30)
	seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)
	
	corrected.m <- seurat.integrated@assays$integrated@data
	
	return(corrected.m)
}

### get GO
getGO <- function(sce.o, row.id, id.type = c("ENSEMBL", "SYMBOL"), species = c("mouse", "human"), ncomponents = 20, seed = 1000,
									runSeuratBy = "sample") {
	require(scater)
	species <- match.arg(species)
	id.type <- match.arg(id.type)
	if (species == "mouse") {
		require(org.Mm.eg.db)
		db.o <- "org.Mm.eg.db"
	} else {
		require(org.Hs.eg.db)
		db.o <- "org.Hs.eg.db"
	}
	if (nrow(sce.o) != length(row.id)) stop("row.id does not match nrow sce.o!")
	cycle.anno <- AnnotationDbi::select(get(db.o), keytype="GOALL", keys="GO:0007049", 
																			columns=id.type)[, id.type]
	GO.o <- sce.o[row.id %in% cycle.anno, ]
	### filter velocity genes
	# altExp(GO.o, "velocity") <- altExp(GO.o,  "velocity")[which(rownames(altExp(GO.o, "velocity")) %in% rownames(GO.o)), ]
	
	set.seed(seed)
	GO.o <- runPCA(GO.o, exprs_values = "log.s",  name = "PCA.s", ncomponents = ncomponents)
	GO.o <- runUMAP(GO.o, dimred = "PCA.s", exprs_values = "log.s", external_neighbors = FALSE,  pca = ncomponents,  name = "UMAP.s")
	# GO.o <- runPCA(GO.o, exprs_values = "log.u",  name = "PCA.u", ncomponents = ncomponents)
	# GO.o <- runUMAP(GO.o, dimred = "PCA.u", exprs_values = "log.u", external_neighbors = FALSE,  pca = ncomponents,  name = "UMAP.u")
	
	### merge samples
	if (!is.null(runSeuratBy)) {
		set.seed(seed)
		corrected.m <- seuratIntegrate(count.m = assay(GO.o, "spliced"), batch = colData(GO.o)[, runSeuratBy], nfeatures = 500)
		reducedDim(GO.o, "matched.PCA.s") <- calculatePCA(corrected.m, ncomponents = ncomponents)
		metadata(GO.o)$cc.corrected.m <- corrected.m
	}
	return(GO.o)
}



### cyclone
runCyclone <- function(sce.o, gname, assay.type = 'spliced', species = c("mouse", "human"), id.type = c("ensembl", "name")) {
	require(scran)
	species <- match.arg(species)
	if (species == "mouse") {
		pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
		AnnotationDb <- org.Mm.eg.db::org.Mm.eg.db
	} else {
		pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
		AnnotationDb <- org.Hs.eg.db::org.Hs.eg.db
	}
	if (id.type == "name") {
		pairs <- lapply(pairs, function(x) {
			x$first <- suppressMessages(AnnotationDbi::mapIds(AnnotationDb, keys = x$first, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
			x$second <- suppressMessages(AnnotationDbi::mapIds(AnnotationDb, keys = x$second, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
			retain.idx <- which(!rowAnyNAs(as.matrix(x)))
			return(x[retain.idx, ])
		})
		
	}
	cyclone.o <- cyclone(sce.o, pairs = pairs, gene.names = gname, assay.type = assay.type)
	out <- factor(cyclone.o$phases, levels = c("G1", "S", "G2M"))
	return(out)
}

### seurat cc
runSeuratCC <- function(sce.o, gname, assay.type = 'spliced', species = c("mouse", "human"), id.type = c("ensembl", "name")) {
	require(Seurat)
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes
	species <- match.arg(species)
	if (species == "mouse") {
		AnnotationDb <- org.Mm.eg.db::org.Mm.eg.db
	} else {
		AnnotationDb <- org.Hs.eg.db::org.Hs.eg.db
	}
	
	if (id.type == "ensembl") {
		gname <- suppressMessages(AnnotationDbi::mapIds(AnnotationDb, keys = x$first, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
	}
	counts.m <- assay(sce.o, assay.type)
	rownames(counts.m) <- toupper(gname)
	
	### remove NA gene names
	idx <- which(!is.na(rownames(counts.m)))
	o <- CreateSeuratObject(counts = counts.m[idx, ])
	o <- NormalizeData(o)
	o <- FindVariableFeatures(o, selection.method = "vst")
	o<- ScaleData(o, features = rownames(o))
	o <- CellCycleScoring(o, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
	return(factor(o$Phase, levels = c("G1", "S", "G2M")))
}

### revelio
runRevelio <- function(sce.o, GENENAME.v, colnames.v, assay.type = "spliced", ccPhaseAssignBasedOnIndividualBatches = TRUE, ...) {
	require(Revelio)
	counts.m <- as.matrix(assay(sce.o, assay.type))
	rownames(counts.m) <- toupper(GENENAME.v)
	colnames(counts.m) <- colnames.v
	### remove NA gene names
	idx <- which(!is.na(rownames(counts.m)))
	
	Revelio.o <- createRevelioObject(rawData = counts.m[idx, ],
																	 cyclicGenes = revelioTestData_cyclicGenes, ccPhaseAssignBasedOnIndividualBatches = ccPhaseAssignBasedOnIndividualBatches, ...)
	
	Revelio.o <- getCellCyclePhaseAssignInformation(dataList = Revelio.o)
	Revelio.o <- getPCAData(dataList = Revelio.o, boolPlotResults = FALSE)
	Revelio.o <- getOptimalRotation(dataList = Revelio.o, boolPlotResults = FALSE)
	dc.m <- t(Revelio.o@transformedData$dc$data[1:2, ])
	return(list(dc = dc.m, cc = Revelio.o@cellInfo$ccPhase))
}

### run peco
runPeco <- function(sce.o, training.l, GENE, assay.type = "spliced", FitGene.v = c("CDK1", "UBE2C", "TOP2A", "HIST1H4E", "HIST1H4C", "UBC", "H4C3", "H4C5")) {
	require(peco)
	counts.m <- as.matrix(assay(sce.o, assay.type))
	rownames(counts.m) <- toupper(GENE)
	### remove NA gene names
	idx <- which(!is.na(rownames(counts.m)))
	sce.o <- SingleCellExperiment(assays = list(counts = counts.m[idx, ])) %>% data_transform_quantile(ncores = 10)
	
	genes.v <- intersect(rownames(sce.o), rownames(training.l$predict.yy))
	
	pred.o <-
		cycle_npreg_outsample(Y_test = sce.o[genes.v, ],
													sigma_est = training.l$sigma[genes.v,],
													funs_est = training.l$cellcycle_function[genes.v],
													method.trend = "trendfilter", get_trend_estimates = FALSE, ncores = 10)
	
	theta_predict <- colData(pred.o$Y)$cellcycle_peco
	names(theta_predict) <- rownames(colData(pred.o$Y))
	g.v <- FitGene.v[FitGene.v %in% rownames(pred.o$Y)]
	yy_input <- assay(pred.o$Y,"cpm_quantNormed")[g.v,, drop = FALSE]
	if (nrow(yy_input) == 1) {
		Y_ordered <- yy_input[, names(theta_predict), drop = FALSE]
		ord <- order(theta_predict)
		theta_ordered <- theta_predict[ord]
		y_g <- Y_ordered[, ord, drop = TRUE]
		fit_g <- fit_trendfilter(yy = y_g, polyorder = 2)
		fun_g <- approxfun(x = as.numeric(theta_ordered), 
											 y = as.numeric(fit_g$trend.yy), rule = 2)
		mu_g <- fit_g$trend.yy
		sigma_g <- sqrt(sum((y_g - mu_g)^2)/length(theta_predict))
		predict.yy <- matrix(fit_g$trend.yy, nrow = 1)
		colnames(predict.yy) <- colnames(Y_ordered)
		rownames(predict.yy) <- rownames(Y_ordered)
		sigma <- data.frame(sigma = sigma_g)
		rownames(sigma) <- rownames(Y_ordered)
		pve <- data.frame(pve = fit_g$pve)
		rownames(pve) <- rownames(Y_ordered)
		cellcycle_function <- list(fun_g)
		names(cellcycle_function) <- rownames(Y_ordered)
		fit_cyclic <- list(predict.yy = predict.yy, 
				 cellcycle_peco_ordered = theta_ordered, 
				 cellcycle_function = cellcycle_function, sigma = sigma, 
				 pve = pve)
	} else {
		fit_cyclic <- fit_cyclical_many(Y = yy_input, theta = theta_predict, ncores = 10)
	}
	

	return(list(theta_predict = theta_predict, g.v = g.v, yy_input = yy_input, fit_cyclic = fit_cyclic))
}

### run reCAT
runreCAT <- function(sce.o, gname = NULL, assay.type = 'log.s') {
	require(reCAT)
	data <- as.matrix(assay(sce.o, assay.type))
	if (!is.null(gname)) rownames(data) <- gname
	test <- get_test_exp(data)
	cycle_scores <- get_score(t(test))
	# About 40 minutes
	idx <- which(rowSums(test) > 0)
	step2_ordIndex <- rep(NA, ncol(data))
	step2_ordIndex[idx] <- get_ordIndex(test[idx, ], threadnum = 10, step_size = 2, clust_method = NULL)
	
	return(list(scores = cycle_scores, ordIndex = step2_ordIndex ))
}



plotEmbScat <- function(sce.o, dimred, 
											 color_by,
											 colors.v,
											 labels.v,
											 color.name,
												x_lab = NULL,
												y_lab = NULL,
												facet_by = NULL,
												facet_labels = NULL,
												title = NULL) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " (n=", ncol(sce.o), ")")
	if(is.null(x_lab)) x_lab <- px_lab
	if(is.null(y_lab)) y_lab <- py_lab
	
	emb.m <- reducedDim(sce.o, dimred)

	x_lim <- range(emb.m[, 1]) + c(diff(range(emb.m[, 1])) * c(-0.05, 0.05))
	y_lim <- range(emb.m[, 2]) + c(diff(range(emb.m[, 2])) * c(-0.05, 0.05))
	
	tmp.df <- data.frame(x = emb.m[, 1],
											 y = emb.m[, 2],
											 color = colData(sce.o)[, color_by])
	if (any(is.na(tmp.df$color))) {
		tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
		scale_color <- scale_color_manual(values = c(colors.v, "grey"), name = color.name, labels =  c(labels.v, "NA"), limits =   c(labels.v, "NA"))
	} else {
		scale_color <- scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = labels.v)
	}
	
	
	scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size,  alpha = point.alpha) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size,  alpha = point.alpha) +
		scale_color +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( y = y_lab, x = x_lab, title = title) +
		xlim(x_lim) + ylim(y_lim)
	
	if (!is.null(facet_by)) {
		if (is.null(facet_labels)) facet_labels <- levels(factor(colData(sce.o)[, facet_by]))
		tmp.df$facet <- factor(colData(sce.o)[, facet_by])
		lp <- lapply(seq_len(nlevels(factor(tmp.df$facet))), function(idx) {
			p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
				geom_scattermore(data = tmp.df %>% dplyr::filter(`facet` != levels(factor(tmp.df$facet))[idx]), pointsize = point.size,  alpha = 0.4, color = "gray90", show.legend = FALSE) +
				geom_scattermore(data = tmp.df %>% dplyr::filter(`facet` == levels(factor(tmp.df$facet))[idx]), pointsize = point.size,  alpha = point.alpha) +
				scale_color +
				guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
				labs( y = y_lab, x = x_lab, title = str_c(facet_labels[idx], " (n=", sum(tmp.df$facet == levels(factor(tmp.df$facet))[idx], na.rm = TRUE), ")")) +
				xlim(x_lim) + ylim(y_lim)
			return(p)
		})
		return(c(list(scat.p), lp))
	}
	return(scat.p)
}

plotEmbScatViridis <- function(sce.o, dimred, 
															 color_by,
															facet_by = NULL,
															facet_labels = NULL,
															color.name,
															title = NULL,
															x_lab = NULL,
															y_lab = NULL) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " (n=", ncol(sce.o), ")")
	if(is.null(x_lab)) x_lab <- px_lab
	if(is.null(y_lab)) y_lab <- py_lab
	
	emb.m <- reducedDim(sce.o, dimred)
	
	x_lim <- range(emb.m[, 1]) + c(diff(range(emb.m[, 1])) * c(-0.05, 0.05))
	y_lim <- range(emb.m[, 2]) + c(diff(range(emb.m[, 2])) * c(-0.05, 0.05))
	
	tmp.df <- data.frame(x = emb.m[, 1],
											 y = emb.m[, 2],
											 color = colData(sce.o)[, color_by])
	
	scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(pointsize = point.size,  alpha = point.alpha) +
		scale_color_viridis(name = color.name, limits = range(tmp.df$color, na.rm = TRUE)) + 
		labs( y = y_lab, x = x_lab) +
		ggtitle(title) +
		xlim(x_lim) + ylim(y_lim)
	
	if (!is.null(facet_by)) {
		if (is.null(facet_labels)) facet_labels <- levels(factor(colData(sce.o)[, facet_by]))
		tmp.df$facet <- colData(sce.o)[, facet_by]
		lp <- lapply(seq_len(nlevels(factor(tmp.df$facet))), function(idx) {
			p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
				geom_scattermore(data = tmp.df %>% dplyr::filter(`facet` != levels(factor(tmp.df$facet))[idx]), pointsize = point.size,  alpha = 0.4, color = "gray90", show.legend = FALSE) +
				geom_scattermore(data = tmp.df %>% dplyr::filter(`facet` == levels(factor(tmp.df$facet))[idx]), pointsize = point.size,  alpha = point.alpha) +
				scale_color_viridis(name = color.name, limits = range(tmp.df$color, na.rm = TRUE)) + 
				labs( y = y_lab, x = x_lab, title = facet_labels[idx]) +
				xlim(x_lim) + ylim(y_lim)
			return(p)
		})
		return(c(list(scat.p), lp))
	}
	return(scat.p)
}


plotEmbScatCyclic <- function(sce.o, dimred, 
														 color_by = "tricyclePosition",
														 x_lab = NULL,
														 y_lab = NULL,
														 facet_by = NULL,
														 facet_labels = NULL,
														 color.name = NULL,
														 title = NULL,
														 emb.name = c("PCA", "UMAP", "FUCCI", "FUCCI2"),
														 hue.colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"),
														 hue.n = 500,
														 plot.legend = FALSE) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " (n=", ncol(sce.o), ")")
	if(is.null(x_lab)) x_lab <- px_lab
	if(is.null(y_lab)) y_lab <- py_lab
	
	emb.m <- reducedDim(sce.o, dimred)
	
	x_lim <- range(emb.m[, 1]) + c(diff(range(emb.m[, 1])) * c(-0.05, 0.05))
	y_lim <- range(emb.m[, 2]) + c(diff(range(emb.m[, 2])) * c(-0.05, 0.05))
	
	tmp.df <- data.frame(x = emb.m[, 1],
											 y = emb.m[, 2],
											 color = colData(sce.o)[, color_by])
	
	scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(pointsize = point.size,  alpha = point.alpha) +
		scale_color_gradientn(name = color.name, limits = range(0, 2 * pi), 
													breaks = seq(from = 0, to = 2 * pi, length.out = hue.n) ,
													colors = hue.colors, 
													guide = FALSE) + 
		labs( y = y_lab, x = x_lab) +
		ggtitle(title) +
		xlim(x_lim) + ylim(y_lim)
	
	if (!is.null(facet_by)) {
		if (is.null(facet_labels)) facet_labels <- levels(factor(tmp.df$facet))
		tmp.df$facet <- colData(sce.o)[, facet_by]
		lp <- lapply(seq_len(nlevels(factor(tmp.df$facet))), function(idx) {
			p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
				geom_scattermore(data = tmp.df %>% dplyr::filter(`facet` != levels(factor(tmp.df$facet))[idx]), pointsize = point.size,  alpha = 0.4, color = "gray90", show.legend = FALSE) +
				geom_scattermore(data = tmp.df %>% dplyr::filter(`facet` == levels(factor(tmp.df$facet))[idx]), pointsize = point.size,  alpha = point.alpha) +
				scale_color_gradientn(name = color.name, limits = range(0, 2 * pi), 
															breaks = seq(from = 0, to = 2 * pi, length.out = hue.n) ,
															colors = hue.colors, 
															guide = FALSE) +
				labs( y = y_lab, x = x_lab, title = facet_labels[idx]) +
				xlim(x_lim) + ylim(y_lim)
			return(p)
		})
		return(c(list(scat.p), lp))
	}
	return(scat.p)
}

plotScatCC <- function(sce.o, dimred, x_lab, y_lab, title = NULL, col.name = "CCStage", color.name = "SchwabeCC", colors.v = NULL, labels.v = NULL) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if(is.null(colors.v)) colors.v <- ccColors.v
	if(is.null(labels.v)) labels.v <- ccLabels.v
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " (n=", ncol(sce.o), ")")
	
	emb.m <- reducedDim(sce.o, dimred)
	tmp.df <- data.frame(pc1 = emb.m[, 1], pc2 = emb.m[, 2]) %>% add_column(cc = colData(sce.o)[, col.name])
	if (any(is.na(tmp.df$color))) {
		tmp.df$cc <- fct_explicit_na(tmp.df$cc, na_level = "NA") %>% fct_relevel("NA", after = Inf)
		scale_color <- scale_color_manual(values = c(colors.v, "grey"), name = color.name, labels =  c(labels.v, "NA"), limits =   c(labels.v, "NA"))
	} else {
		scale_color <- scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = labels.v)
	}
	
	cc.p <- ggplot(tmp.df, aes(x = pc1, y = pc2 , color = cc)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`cc` == "NA"), pointsize = point.size,  alpha = point.alpha) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`cc` != "NA"), pointsize = point.size,  alpha = point.alpha) +
		scale_color +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( y = y_lab, 
					x = x_lab, 
					title = title) 
	
	return(cc.p)
}


plotCirclePlotProjection <- function(sce.o, r, label.x, label.y,  dimred = NULL, x_lab = NULL, y_lab = NULL, col.name = "CCStage", color.name = "SchwabeCC", colors.v = NULL, labels.v = NULL,  label.size = 3, title = NULL) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " (n=", ncol(sce.o), ")")
	if(is.null(colors.v)) colors.v <- ccColors.v
	if(is.null(labels.v)) labels.v <- ccLabels.v
	if(is.null(x_lab)) x_lab <- px_lab
	if(is.null(y_lab)) y_lab <- py_lab
	if(is.null(dimred)) dimred <- "tricycleEmbedding"
	
	emb.m <- reducedDim(sce.o, dimred)
	tmp.df <- data.frame(pc1 = emb.m[, 1], pc2 = emb.m[, 2]) %>% add_column(color = colData(sce.o)[, col.name])
	if (any(is.na(tmp.df$color))) {
		tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
		scale_color <- scale_color_manual(values = c(colors.v, "grey"), name = color.name, labels =  c(labels.v, "NA"), limits =   c(labels.v, "NA"))
	} else {
		scale_color <- scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = labels.v)
	}
	xlim <- c(min(-r, min(tmp.df$pc1)), max(r, max(tmp.df$pc1)))
	ylim <- c(min(-r, min(tmp.df$pc2)), max(r, max(tmp.df$pc2)))
	
	cc.p <- ggplot(tmp.df, aes(x = pc1, y = pc2 , color = color)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size,  alpha = point.alpha) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size,  alpha = point.alpha) +
		scale_color +
		#scale_shape_manual(name = "Cellcycle", values = c(16, 17, 18)) +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		annotate("point", x = 0, 
						 y = 0, shape = 16, size = 0.8, alpha = 0.8) + 
		annotate("segment", x = 0, xend =  r ,
						 y = 0, yend = 0, linetype = "dashed", alpha = 0.8) + 
		annotate("segment", x = 0, xend =  cos(pi/6) * r,
						 y = 0, yend =  sin(pi/6) * r,
						 linetype = "dashed", alpha = 0.8) +
		annotate("text",  alpha = 0.8, x = label.x, y = label.y, label = "\u03B8", parse = TRUE, size = label.size) +
		xlim(xlim) + ylim(ylim) +
		labs( y = y_lab, 
					x = x_lab, 
					title = title) 
	
	return(cc.p)
}

plotCirclePlotProjection2 <- function(sce.o, r, label.x, label.y,  dimred = NULL, x_lab = NULL, y_lab = NULL, col.name = "tricyclePosition", color.name = "tricycle", colors.v = NULL, hue.n = 500, label.size = 3, title = NULL) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " (n=", ncol(sce.o), ")")
	if(is.null(colors.v)) colors.v <- colors.v <- c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9")
	
	if(is.null(x_lab)) x_lab <- px_lab
	if(is.null(y_lab)) y_lab <- py_lab
	if(is.null(dimred)) dimred <- "tricycleEmbedding"
	
	emb.m <- reducedDim(sce.o, dimred)
	tmp.df <- data.frame(pc1 = emb.m[, 1], pc2 = emb.m[, 2]) %>% add_column(color = colData(sce.o)[, col.name])
	scale_color <- scale_color_gradientn(name = color.name, limits = range(0, 2 * pi), 
																			 breaks = seq(from = 0, to = 2 * pi, length.out = hue.n) ,
																			 colors = colors.v, 
																			 guide = FALSE)
	
	xlim <- c(min(-r, min(tmp.df$pc1)), max(r, max(tmp.df$pc1)))
	ylim <- c(min(-r, min(tmp.df$pc2)), max(r, max(tmp.df$pc2)))
	
	cc.p <- ggplot(tmp.df, aes(x = pc1, y = pc2 , color = color)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size,  alpha = point.alpha) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size,  alpha = point.alpha) +
		scale_color +
		annotate("point", x = 0, 
						 y = 0, shape = 16, size = 0.8, alpha = 0.8) + 
		annotate("segment", x = 0, xend =  r ,
						 y = 0, yend = 0, linetype = "dashed", alpha = 0.8) + 
		annotate("segment", x = 0, xend =  cos(pi/6) * r,
						 y = 0, yend =  sin(pi/6) * r,
						 linetype = "dashed", alpha = 0.8) +
		annotate("text",  alpha = 0.8, x = label.x, y = label.y, label = "\u03B8", parse = TRUE, size = label.size) +
		xlim(xlim) + ylim(ylim) +
		labs( y = y_lab, 
					x = x_lab, 
					title = title) 
	
	return(cc.p)
}

plotLoess <- function(sce.o, col.name, col.outname = NULL, title = NULL, x_val = "tricyclePosition", x_lab = NULL, color.var = "CCStage", color.name = "SchwabeCC", colors.v = NULL, labels.v = NULL, log2.trans = FALSE, y_lab = NULL, addR2 = FALSE, r2size = 2) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (is.null(col.outname)) {
		if (metadata(sce.o)$species == "mouse") col.outname <- str_to_title(col.name)
		if (metadata(sce.o)$species == "human") col.outname <- str_to_upper(col.name)
	}
	
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " ", col.outname, " (n=", ncol(sce.o), ")")
	if (is.null(colors.v)) colors.v <- ccColors.v
	if (is.null(labels.v)) labels.v <- ccLabels.v
	if (is.null(x_lab)) x_lab <- theta_lab
	if (is.null(y_lab)) y_lab <- bquote(paste('log'['2'],'(expression of ', .(col.outname), ")"))
	theta.v <- colData(sce.o)[, x_val]
	y <- colData(sce.o)[, col.name]
	if (log2.trans) y <- log2(y)
	tmp.df <- data.frame(theta = theta.v, color = colData(sce.o)[, color.var], y = y)
	if (any(is.na(tmp.df$color))) {
		tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
		scale_color <- scale_color_manual(values = c(colors.v, "grey"), name = color.name, labels =  c(labels.v, "NA"), limits =   c(labels.v, "NA"))
	} else {
		scale_color <- scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = labels.v)
	}
	
	loess.l <- fit_periodic_loess(theta.v = theta.v, y = y)
	
	p <- ggplot(data = tmp.df , aes(x = theta , y = y, color = color)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size, alpha = point.alpha) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size, alpha = point.alpha) +
		geom_path(data = loess.l$pred.df, aes(x = x , y = y), linetype = "dashed", color = "black", size = 0.8, alpha = 0.6, inherit.aes = FALSE) +
		scale_color +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( y = y_lab , x = x_lab, title = title) +
		scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels =c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi) 
	if (addR2) p <- p + annotate(geom = "text", x = 0, y = .percent_range(tmp.df$y, 1), size = r2size, hjust = 0, vjust = 1,
															 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE)
	
	return(p)
}

plotLoess2 <- function(sce.o, col.name, col.outname = NULL, title = NULL, x_val = "tricyclePosition", x_lab = NULL, color.var = "tricyclePosition", color.name = "tricycle", colors.v = NULL, hue.n = 500, log2.trans = FALSE, y_lab = NULL, addR2 = FALSE, r2size = 2) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (is.null(col.outname)) {
		if (metadata(sce.o)$species == "mouse") col.outname <- str_to_title(col.name)
		if (metadata(sce.o)$species == "human") col.outname <- str_to_upper(col.name)
	}
	
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " ", col.outname, " (n=", ncol(sce.o), ")")
	if (is.null(colors.v)) colors.v <- c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9")
	if (is.null(x_lab)) x_lab <- theta_lab
	if (is.null(y_lab)) y_lab <- bquote(paste('log'['2'],'(expression of ', .(col.outname), ")"))
	theta.v <- colData(sce.o)[, x_val]
	y <- colData(sce.o)[, col.name]
	if (log2.trans) y <- log2(y)
	tmp.df <- data.frame(theta = theta.v, color = colData(sce.o)[, color.var], y = y)
	scale_color <- scale_color_gradientn(name = color.name, limits = range(0, 2 * pi), 
																			 breaks = seq(from = 0, to = 2 * pi, length.out = hue.n) ,
																			 colors = colors.v, 
																			 guide = FALSE)
	
	loess.l <- fit_periodic_loess(theta.v = theta.v, y = y)
	
	p <- ggplot(data = tmp.df , aes(x = theta , y = y, color = color)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size, alpha = point.alpha) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size, alpha = point.alpha) +
		geom_path(data = loess.l$pred.df, aes(x = x , y = y), linetype = "dashed", color = "black", size = 0.8, alpha = 0.6, inherit.aes = FALSE) +
		scale_color +
		labs( y = y_lab , x = x_lab, title = title) +
		scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels =c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi) 
	if (addR2) p <- p + annotate(geom = "text", x = 0, y = .percent_range(tmp.df$y, 1), size = r2size, hjust = 0, vjust = 1,
															 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE)
	
	return(p)
}

plotLoess3 <- function(sce.o, col.name, col.outname = NULL, title = NULL, x_val = "tricyclePosition", x_lab = NULL, color.var = "tricyclePosition", color.name = "tricycle", colors.v = NULL, hue.n = 500, log2.trans = FALSE, y_lab = NULL, addR2 = FALSE, r2size = 2) {
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (is.null(col.outname)) {
		if (metadata(sce.o)$species == "mouse") col.outname <- str_to_title(col.name)
		if (metadata(sce.o)$species == "human") col.outname <- str_to_upper(col.name)
	}
	
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " ", col.outname, " (n=", ncol(sce.o), ")")
	if (is.null(colors.v)) colors.v <- c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9")
	if (is.null(x_lab)) x_lab <- theta_lab
	if (is.null(y_lab)) y_lab <- bquote(paste('log'['2'],'(expression of ', .(col.outname), ")"))
	theta.v <- colData(sce.o)[, x_val]
	y <- colData(sce.o)[, col.name]
	if (log2.trans) y <- log2(y)
	tmp.df <- data.frame(theta = theta.v, color = colData(sce.o)[, color.var], y = y)
	scale_color <- scale_color_gradientn(name = color.name, limits = range(0, 2 * pi), 
																			 breaks = seq(from = 0, to = 2 * pi, length.out = hue.n) ,
																			 colors = colors.v, 
																			 guide = FALSE)
	
	loess.l <- fit_periodic_loess(theta.v = theta.v, y = y)
	
	p <- ggplot(data = tmp.df , aes(x = theta , y = y)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size, alpha = point.alpha, color = "grey") +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size, alpha = point.alpha, color = "grey") +
		geom_path(data = loess.l$pred.df, aes(x = x , y = y), linetype = "dashed", color = "black", size = 0.8, alpha = 1, inherit.aes = FALSE) +
		labs( y = y_lab , x = x_lab, title = title) +
		scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi, labels =c(0, str_c(seq(0.5, 2, 0.5), "\u03C0")), limits = c(0, 2) * pi) 
	if (addR2) p <- p + annotate(geom = "text", x = 0, y = .percent_range(tmp.df$y, 1), size = r2size, hjust = 0, vjust = 1,
															 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE)
	
	return(p)
}


plotGeneRidge <- function(sce.o, col.name, col.outname = NULL, title = NULL,  x_lab = NULL, color.var = "CCStage", color.name = "SchwabeCC", colors.v = NULL, labels.v = NULL, log2.trans = FALSE, y_lab = NULL) {
	require(ggridges)
	point.size <- metadata(sce.o)$point.size
	point.alpha <- metadata(sce.o)$point.alpha
	if (metadata(sce.o)$species == "mouse") {
		if (col.name == "top2a") {
			col.outname <- "Top2A"
		}  else {
			col.outname <- str_to_title(col.name)
		}
	} else {
		col.outname <- str_to_upper(col.name)
	} 
	
	
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " ", col.outname, " (n=", ncol(sce.o), ")")
	if (is.null(colors.v)) colors.v <- ccColors.v
	if (is.null(labels.v)) labels.v <- ccLabels.v
	if (is.null(x_lab)) x_lab <- bquote(paste('log'['2'],'(expression of ', .(col.outname), ")"))

	x <- colData(sce.o)[, col.name]
	if (log2.trans) y <- log2(y)
	tmp.df <- data.frame(color = colData(sce.o)[, color.var], x = x)
	if (any(is.na(tmp.df$color))) {
		tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
		scale_color <- scale_color_manual(values = c(colors.v, "grey"), name = color.name, labels =  c(labels.v, "NA"), limits =   c(labels.v, "NA"))
		scale_fill <- scale_fill_manual(values = c(colors.v, "grey"), name = color.name, labels =  c(labels.v, "NA"), limits =   c(labels.v, "NA"))
	} else {
		scale_color <- scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = labels.v)
		scale_fill <- scale_fill_manual(values = colors.v, name = color.name, labels = labels.v, limits = labels.v)
	}
	
	p <- ggplot(data = tmp.df , aes(x = x, y = color, fill = color, color = color)) +
		geom_density_ridges(scale = 3, rel_min_height = 0.01, alpha = 0.6) +
		scale_color + scale_fill +
		labs( y = y_lab , x = x_lab, title = title) 
	
	return(p)
}


plotThetaDen <- function(sce.o, color.var = "CCStage", color.name = "SchwabeCC", colors.v = NULL, labels.v = NULL, title = NULL, bw = 30, scale = FALSE, addP = FALSE, adjustbw = TRUE) {
	if (is.null(title)) title <- metadata(sce.o)$dataname
	if (is.null(colors.v)) colors.v <- ccColors.v
	if (is.null(labels.v)) labels.v <- ccLabels.v
	
	theta <- sce.o$tricyclePosition
	color <- factor(colData(sce.o)[, color.var], levels = labels.v)
	
	if (addP) {
		title <- str_c("aov Pval: ", format.pval(aov.circular(circular(theta), color)$p.value[1], digits = 3))
	}
	if (adjustbw ) {
		if (ncol(sce.o) < 2000) bw <- 20
		if (ncol(sce.o) < 1000) bw <- 10
		
	}
	
	d <- density(circular(theta), bw = bw)
	all.df <- data.frame(x = as.numeric(d$x), y = d$y)
	
	
	if (scale) {
		alpha <- table(color) / sum(!is.na(color), na.rm = TRUE)
	} else {
		alpha <- rep(1, nlevels(color))
	}
	
	stage.df <- do.call(rbind, lapply(seq_len(nlevels(color)), function(idx) {
		d <- density(circular(theta[which(color == levels(color)[idx])]), bw = bw)
		return(data.frame(x = as.numeric(d$x), y = d$y * alpha[idx], color = levels(color)[idx]))
	}))
	stage.df$color <- factor(stage.df$color, levels = levels(color))
	max.v <- max(all.df$y, stage.df$y)
	
	
	l.den.p <- ggplot(stage.df, aes(x = x  , y = y )) +
		geom_path(aes(color = color), size = 0.5, alpha = 0.5) +
		geom_path(data = all.df, size = 0.5, alpha = 0.5, color = "black", linetype = "dashed") +
		# geom_jitter(data = tmp.df, aes(x = theta1 * pi + pi / 2, y = tmp.l[[1]]$max * 0.5, color = cc), size = 0.1,  alpha = 0.5, shape = 16) + 
		scale_color_manual(values = colors.v, name = color.name, labels =  labels.v) +
		scale_x_continuous(limits = c(0, 2 * pi),
											 breaks =  c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
											 labels = str_c(c(0, 0.5, 1, 1.5, 2), "\u03C0"),
											 name = "\u03B8") +
		labs(title = title, y = "Density") 
	
	cir.den.p <- ggplot(stage.df, aes(x = x  , y = y + max.v )) +
		geom_path(aes(color = color), size = 0.5, alpha = 0.5) +
		geom_path(data = all.df, size = 0.5, alpha = 0.5, color = "black", linetype = "dashed") +
		# geom_jitter(data = tmp.df, aes(x = theta1 * pi + pi / 2, y = tmp.l[[1]]$max * 0.5, color = cc), size = 0.1,  alpha = 0.5, shape = 16) + 
		scale_color_manual(values = colors.v, name = color.name, labels =  labels.v) +
		coord_polar(theta = "x", start = - pi / 2, direction = -1, clip = "on") +
		scale_x_continuous(limits = c(0, 2 * pi),
											 breaks =  c(0, pi / 2, pi, 3 * pi / 2),
											 labels = str_c(c(0, 0.5, 1, 1.5), "\u03C0"),
											 name = "") +
		scale_y_continuous(limits = c(0, max.v * 2),
											 breaks =  c(max.v, max.v * 1.5, max.v * 2),
											 labels = c("0", format(max.v * c(0.5, 1), digits = 3)),
											 name = "Density") +
		labs(title = title) 
	
	return(list(linear = l.den.p, circular = cir.den.p))
	
}



plotSilhouetteBox <- function(sce.o, title = NULL,  color.var = "CCStage", color.name = "SchwabeCC", colors.v = NULL, labels.v = NULL) {
	require(cluster)
	if (is.null(colors.v)) colors.v <- ccColors.v
	if (is.null(labels.v)) labels.v <- ccLabels.v
	
	theta.v <- sce.o$tricyclePosition
	silhouette.v <- cluster::silhouette(as.integer(factor(colData(sce.o)[, color.var], levels = labels.v)),
																			dist.circular(circular(theta.v), method = "angularseparation"))[, 3]
	
	
	tmp.df <- data.frame(silhouette = silhouette.v, color = factor(colData(sce.o)[, color.var], levels = labels.v))
	if (is.null(title)) title <- str_c( metadata(sce.o)$dataname, " Mean:", format(mean(tmp.df$silhouette), digits = 2))
	
	p <- ggplot(tmp.df, aes(x = color, y = silhouette, color = color)) +
		geom_boxplot(outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.8, fill = "white") +
		geom_hline(yintercept = mean(tmp.df$silhouette), linetype = "dashed", size = 0.8) +
		scale_color_manual(values = colors.v, name = color.name, labels = labels.v) +
		# geom_quasirandom(width = 0.08, size = 0.2, shape = 16, alpha = 0.5, dodge.width = 0.8) +
		scale_x_discrete(name = color.name, labels = str_c(levels(tmp.df$color), "(n=", table(tmp.df$color), ")")) +
		labs( y =  "Silhouette index",  title = title) +
		ylim(c(-1, 1)) 
		#annotate(geom = "text", x = 2, y = mean(tmp.df$silhouette), label = str_c("Mean:", format(mean(tmp.df$silhouette), digits = 2)), size = 3, hjust = 0.5, vjust = - 0.8)
	return(p)
}




### for velocity figure
projectPCA <- function(new.m, rotation.m, rownames.v = NULL) {
	if (!is.null(rownames.v)) rownames(new.m) <- toupper(rownames.v)
	genes <- intersect(rownames(new.m), rownames(rotation.m))
	rotation.m <- rotation.m[genes, ]
	new.m <- new.m[genes, ]
	return(scale(t(new.m), center = T, scale = F) %*% rotation.m)
}

computeVelocityOnGrid <- function(X_emb, V_emb, x1range, x2range, density = 1, smooth = 0.5, adjust = TRUE, n_neighbors = NULL, min_mass = 1, cutoff_perc = 5) {
	# remove invalid cells
	valid.idx <- is.finite(rowSums(X_emb) + rowSums(V_emb))
	X_emb <- X_emb[valid.idx, ]
	V_emb <- V_emb[valid.idx, ]
	
	# prepare gird
	X_1 <- seq(from = x1range[1] - 0.01 * abs(x1range[1] - x1range[2]), to = x1range[2] + 0.01 * abs(x1range[1] - x1range[2]), length.out = as.integer(50 * density))
	X_2 <- seq(from = x2range[1] - 0.01 * abs(x2range[1] - x2range[2]), to = x2range[2] + 0.01 * abs(x2range[1] - x2range[2]), length.out = as.integer(50 * density))
	meshgrid.lm <- pracma::meshgrid(X_1, X_2)
	X_grid <- cbind(as.vector(t(meshgrid.lm$X)), as.vector(t(meshgrid.lm$Y)))
	
	# get neighbors in X_emb space
	if (is.null(n_neighbors)) n_neighbors <- as.integer(nrow(X_emb) / 50)
	neighbors.l <- BiocNeighbors::queryKmknn(X = X_emb, query = X_grid, k = n_neighbors, get.index=TRUE, get.distance=TRUE)
	
	scale <- mean(c(X_1[2] - X_1[1], X_2[2] - X_2[1])) * smooth
	weight <- dnorm(neighbors.l$distance, sd = scale)
	p_mass <- rowSums(weight)
	
	
	V_grid <- cbind(sapply(seq_len(nrow(X_grid)), function(i) sum(V_emb[neighbors.l$index[i, ], 1] * weight[i, ])) / pmax(p_mass, 1),
									sapply(seq_len(nrow(X_grid)), function(i) sum(V_emb[neighbors.l$index[i, ], 2] * weight[i, ])) / pmax(p_mass, 1))
	
	if (adjust) {
		mass <- sqrt(rowSums(V_grid ^ 2))
		min_mass <- 10 ^ (min_mass - 6) # default min_mass = 1e-5
		if (min_mass > (max(mass) * 0.9)) min_mass <- max(mass) * 0.9
		cutoff <- mass < min_mass
		length = sapply(seq_len(nrow(X_grid)), function(i) mean(abs(V_emb[neighbors.l$index[i, ], 1]))) + sapply(seq_len(nrow(X_grid)), function(i) mean(abs(V_emb[neighbors.l$index[i, ], 2])))
		cutoff <- cutoff | (length < quantile(length, cutoff_perc / 100))
		V_grid[cutoff, 1] <- NA
		
	} else {
		min_mass <- min_mass * quantile(p_mass, 0.99) / 100
		V_grid[p_mass < min_mass, 1] <- NA
	}
	
	
	return(list(original = data.frame(X_emb), grid = data.frame(x = X_grid[, 1], y = X_grid[, 2], dx = V_grid[, 1], dy = V_grid[, 2])))
}

