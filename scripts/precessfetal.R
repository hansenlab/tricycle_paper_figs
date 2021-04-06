rm(list=ls())
source(here::here("scripts/utils.R"))

cell.lp  <- lapply(c("Intestine", "Kidney", "Pancreas", "Stomach"), function(o) {
	sce.o <- qread(here::here("data/fetal", str_c(o, ".qs")))
	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
	return(colData(sce.o))
})

nuclei.lp <- lapply(c("Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart", "Liver", "Lung", "Muscle", "Placenta", "Spleen", "Thymus"), function(o) {
	sce.o <- qread(here::here("data/fetal/nuclei", str_c(o, ".qs")))
	sce.o <- estimate_cycle_position(sce.o, exprs_values = "log.s", gname.type = "SYMBOL", species = "human")
	return(colData(sce.o))
})

colData.df <- rbind(do.call(rbind, cell.lp), do.call(rbind, nuclei.lp))

global_umap.df <- qread(here::here("data/globalumap_df.qs"))
dim(global_umap.df)
dim(colData.df)
sum(colData.df$sample %in% global_umap.df$sample)


colData.df$Global_umap_1 <- global_umap.df$Global_umap_1[match(colData.df$sample, global_umap.df$sample)]
colData.df$Global_umap_2 <- global_umap.df$Global_umap_2[match(colData.df$sample, global_umap.df$sample)]


### save small data
qsave(colData.df, file = here::here("data", "plotdata", str_c("fetal_colData.qs")))














