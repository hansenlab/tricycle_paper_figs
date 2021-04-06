rm(list=ls())
source(here::here("scripts/utils.R"))


neurosphere.o <- qread(here::here("data/plotdata/neurosphere.qs"))
neurosphere_full.o <- qread(here::here("data/neurosphere.qs"))

### neurosphere
theta.v <- neurosphere.o$tricyclePosition
ref.m <- metadata(neurosphere.o)$rotation

AssessStability <- function(data.m, ref.m, theta0.v, n = 30, lb = 50, length.out = 9, BPPARAM = SerialParam(), seed.base = 160) {
	ref.m <- ref.m[, 1:2]
	n.gene <- length(intersect(rownames(data.m), rownames(ref.m)))
	if (n.gene < 200) stop("Too few genes")
	message(str_c("Number of overlapped genes:", n.gene))
	
	theta0.v <- estimate_cycle_position(data.m, ref.m = ref.m)
	sample.gn <- seq(from = n.gene - lb, to = lb, length.out = length.out)
	
	out.m <- do.call(cbind, lapply(sample.gn, function(g.n) {
		r.v <- unlist(bplapply(seq_len(n), function(i) {
			set.seed(i + seed.base)
			idx <- sample(x = seq_len(n.gene), size = g.n, replace = FALSE)
			theta <- estimate_cycle_position(data.m, ref.m = ref.m[idx, ]) 
			return(Directional::circ.cor1(theta0.v, theta, rads = T)[1])
		}, BPPARAM = BPPARAM))
		return(r.v)
	}))
	colnames(out.m) <- as.character(sample.gn)
	return(out.m)
}

removegeneRho.m <- AssessStability(data.m = assay(neurosphere_full.o, "log.s"),
																	 ref.m = ref.m, 
																	 theta0.v = theta.v,
																	 BPPARAM = MulticoreParam(workers = 12))


### figure
tmp.df <- data.frame(removegeneRho.m) %>% gather(key = "g.n", value = "rho", 1:9)
tmp.df$g.n <- factor(tmp.df$g.n, levels = str_c("X",  seq(from = 450, to = 50, length.out = 9)))

Rho_neurospheres.box.p <- ggplot(tmp.df, aes(x = g.n, y = rho)) +
	geom_boxplot(fill = "white", outlier.shape=NA, outlier.size = 0.2, size = 0.2, width = 0.8) +
	geom_quasirandom(width = 0.2, size = 0.6, shape = 16, alpha = 0.5) +
	labs( y = "Circular correlation \u03C1", 
				x = "Number of genes retianed in the mNeurosphere reference", 
				title = str_c("Tolerance of removing genes \nfrom the mNeuroshoere reference")) +
	scale_x_discrete(labels = as.character(seq(from = 450, to = 50, length.out = 9))) +
	scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), limits = c(min(tmp.df$rho, 0), max(tmp.df$rho, 1)))

save_plot(here::here("figs", "sfigs", "sfig.removeGenesRho.pdf"), Rho_neurospheres.box.p,
					base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 1.2, device = cairo_pdf)




