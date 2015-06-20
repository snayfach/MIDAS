# LIBRARIES & FUNCTIONS
library('RColorBrewer')

shannon_ent <- function(x){
	y <- x[x>0]
	return(-sum(y * log(y,2)))
}

# PATHS
outpath <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/genome_cluster_prop_known.pdf"
group_by <- "population"

# READ: sample metadata
sample_map <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/sample_metadata.txt"
metadata <- read.csv(sample_map, sep='\t', header=T, row.names=1)

# GET: matrix of marker-gene relative coverage
indir <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/microbe_cnv_out"
runs <- list.files(indir)
inpath <- paste(c(indir, runs[1], 'genome_clusters.abundance'), collapse='/')
indata <- read.csv(inpath, sep='\t', header=T)
my_mat <- matrix(nrow=length(runs), ncol=nrow(indata))
colnames(my_mat) <- sort(indata$cluster_id)
rownames(my_mat) <- runs
for (run in list.files(indir)){
	inpath <- paste(c(indir, run, 'genome_clusters.abundance'), collapse='/')
	indata <- read.csv(inpath, sep='\t', header=T, row.names=1)
	my_mat[run, ] <- indata[colnames(my_mat), 'prop_cov']
}

# GET & PLOT: sum of relative coverage per sample
prop_known  <- apply(my_mat, 1, sum)
categories <- metadata[rownames(my_mat), group_by]
color_pallete <- brewer.pal(length(unique(categories)), 'Spectral')
pdf(outpath)
boxplot(
	prop_known ~ categories,
	col = color_pallete,
	ylab='Proportion Known Species',
	cex.lab=1.25,
	cex.axis=1.25
)

# GET: matrix of marker-gene relative abundance
indir <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/microbe_cnv_out"
runs <- list.files(indir)
inpath <- paste(c(indir, runs[1], 'genome_clusters.abundance'), collapse='/')
indata <- read.csv(inpath, sep='\t', header=T)
my_mat <- matrix(nrow=length(runs), ncol=nrow(indata))
colnames(my_mat) <- sort(indata$cluster_id)
rownames(my_mat) <- runs
for (run in list.files(indir)){
	inpath <- paste(c(indir, run, 'genome_clusters.abundance'), collapse='/')
	indata <- read.csv(inpath, sep='\t', header=T, row.names=1)
	my_mat[run, ] <- indata[colnames(my_mat), 'rel_abun']
}

# GET & PLOT: shannon entropy per sample
shannon  <- apply(my_mat, 1, shannon_ent)
boxplot(
	shannon ~ categories,
	col = color_pallete,
	ylab='Shannon Diversity',
	cex.lab=1.25,
	cex.axis=1.25
)
dev.off()
