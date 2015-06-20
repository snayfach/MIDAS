#



# libraries
library('RColorBrewer')

prevalence <- function(x){
	return(length(which(is.finite(x) & x > 1.0)))
}

pancoverage <- function(x){
	return(length(which(x>0.25))/length(x))
}

# PATHS
indir <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/microbe_cnv_out"
outpath <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/genome_cluster_prevalence.pdf"

# READ: sample metadata
group_by <- "population"
sample_map <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/sample_metadata.txt"
metadata <- read.csv(sample_map, sep='\t', header=T, row.names=1)

# READ: cluster metadata
cluster_map <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/genome_clustering2/clusters_combined/top_30/0.035/cluster_annotations.txt"
annotations <- read.csv(cluster_map, sep='\t', header=F, row.names=1)
names(annotations) <- c('name', 'freq', 'size')

# GET: marker-gene coverage and fraction of pangenes present
runs <- list.files(indir)
row <- 0
df <- as.data.frame(matrix(nrow=0, ncol=4))
names(df) <- c('run', 'genome_cluster', 'coverage', 'fraction_present')
for (run in runs){
	cov_dir <- paste(c(indir, run, 'coverage'), collapse='/')
	for (file in list.files(cov_dir)){
		# read data
		inpath <- paste(cov_dir, file, sep='/')
		indata <- read.csv(gzfile(inpath), sep='\t', header=T)
		names(indata) <- c('pangene', 'coverage', 'copy_number')
		indata <- with(indata, indata[coverage > 0, ])
		# store data
		row <- row + 1
		df[row, 'run'] <- run
		df[row, 'genome_cluster'] <- strsplit(file, '\\.')[[1]][1]
		df[row, 'coverage'] <- mean(indata$coverage/indata$copy_number)
		df[row, 'fraction_present'] <- pancoverage(indata$copy_number)
	}
}

# GET: cluster prevalence by population
populations <- sort(unique(metadata[,group_by]))
clusters <- sort(unique(df$genome_cluster))
counts <- matrix(nrow=length(populations), ncol=length(clusters))
colnames(counts) <- clusters
rownames(counts) <- populations
for (pop in populations){
	x <- df[df$run %in% row.names(metadata[metadata[,group_by]==pop,]), ]
	for (cluster in clusters){
		counts[pop, cluster] <- prevalence(x[x$genome_cluster == cluster, 'coverage'])
	}
}
counts <- counts[,order(-colSums(counts))]

# WRITE: cluster prevalence by population
outpath <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/cluster_prevalence.txt"
write.table(t(counts), outpath, sep='\t', quote=T, col.names=T, row.names=T)

# PLOT: side-by-side boxplots of % of genome-cluster pan-genes that are present by population
outpath <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/pangenome_coverage.pdf"
pdf(outpath, width=8, height=9)
par(mar=c(5.1,1,1,8))
x <- df[df$genome_cluster %in% colnames(counts)[1:5], ] # subset dataframe by top 10 most prevalent clusters

genome_clusters <- x$genome_cluster
fraction_present <- x$fraction_present
group_labels <- metadata[x$run, group_by]
color_palette <- brewer.pal(length(populations), 'Spectral')
colors <- rep(color_palette, 5)

boxplot(
	fraction_present ~ group_labels + genome_clusters,
	col=colors,
	at=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24),
	xlab='Fraction Pangenome Covered',
	cex.axis=1.5, cex.lab=2,
	horizontal=TRUE,
	yaxt='n'
)

abline(h=c(5,10,15,20), col='gray')

axis(
	side=4,
	at=c(3,8,13,18,23),
	labels=gsub(" ", "\n", annotations[sort(unique(genome_clusters)), 'name']),
	las=2, cex.axis=1.0,
	tick=FALSE
)

dev.off()


# PLOT: boxplot of top 25 most prevalent genome-clusters; colored by population
outpath <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/cluster_prevalence.pdf"
pdf(outpath, width=26, height=12)

par(mar=c(28,8,1,1))

colors <- brewer.pal(length(populations), 'Spectral')
nplot <- 25
space <- 0.4
cluster_ids <- colnames(counts)[seq(nplot)]
cluster_names <- annotations[cluster_ids,'name']
cluster_sizes <- annotations[cluster_ids,'size']
cluster_labels <- paste(cluster_names, paste(paste(' (',cluster_sizes, sep=''),')',sep=''), sep='')

barplot(
	counts[,seq(nplot)],
	col=colors,
	xaxt='n',
	cex.axis=2,
	cex.lab=2.5,
	ylab='Genome-Cluster Prevalence\n(# samples with >1x coverage)',
	space=space
)

axis(
	side=1,
	at=seq(1, nplot+(nplot*space), 1+space),
	labels=cluster_labels,
	las=2,
	cex.axis=2
	)

legend(
	'topright',
	legend=as.character(populations),
	fill=colors,
	title='Host Population',
	cex=2.5
)

boxplot(
	log(coverage,10) ~ genome_cluster,
	data = df[df$genome_cluster %in% cluster_ids, ],
	ylab = 'Marker Gene Coverage',
	cex.axis=2, cex.lab=2.5,
	col='gray', xaxt='n', yaxt='n', frame=F
)

axis(
	side=2,
	at=seq(-1,2),
	labels=10^seq(-1,2),
	cex.axis=2
)

abline(h=seq(-1,2), col='lightgray', lty='dotted')

dev.off()




