# libraries
library('RColorBrewer')
library('gplots')

# command-line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
inpath <- args[1]
outdir <- args[2]
min_cov <- 1
cn_cutoff <- 0.25

# pangenome profile
profile <- read.csv(gzfile(inpath), sep='\t', header=T, row.names=1, check.names=F)
profile <- profile[profile$marker_cov > min_cov, ]
if (nrow(profile) == 0) quit()

# sample metadata
sample_map <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/sample_metadata.txt"
metadata <- read.csv(sample_map, sep='\t', header=T, row.names=1)
group_by <- "population"

# cluster metadata
cluster_map <- "/mnt/data/work/pollardlab/snayfach/projects/strain_variation/genome_clustering2/clusters_combined/top_30/0.035/cluster_annotations.txt"
annotations <- read.csv(cluster_map, sep='\t', header=F, row.names=1)
names(annotations) <- c('name', 'freq', 'size')
cluster_id <- strsplit(basename(inpath), '\\.')[[1]][1]
cluster_name <- annotations[cluster_id, 'name']

# get row colors
categories <- sort(unique(metadata[,group_by]))
color_pallete <- brewer.pal(length(categories), 'Spectral')
colors <- c()
for (element in metadata[row.names(profile), group_by]){
	colors <- c(colors, color_pallete[which(categories == element)])
}

# threshold pangenome profile
# add a tiny bit of jitter so heatmap.2 doesn't freak out
z <- ifelse(profile[,2] > cn_cutoff, runif(1, 0.999, 1.0), runif(1, 0, 0.001))
for (i in seq(3,ncol(profile))){
	z <- cbind(z, ifelse(profile[,i] > cn_cutoff, runif(1, 0.999, 1.0), runif(1, 0, 0.001)))
}

# plot pangenome profile
outname <- paste(c(cluster_id, 'pangenome_profile', 'png'), collapse='.')
outpath <- paste(outdir, outname, sep='/')
png(outpath, width = 480*2, height = 480)
par(cex.main=2)

x <- heatmap.2(
    z,
    trace='none',
    col='gray.colors',
    labCol='',
	labRow='',
	adjRow=c(22,22),
    margins=c(3,6),
    ylab='',
    xlab='',
    main=cluster_name,
	RowSideColors=colors,
	cexRow=3.0, cexCol=3.0, cex.main=3.0,
	density.info='none',
	keysize=1
	)
mtext(
	text=paste('Samples (N=',paste(nrow(profile),')',sep=''), sep=''),
	side=4,
	line=-1,
	cex=1.5
	)
mtext(
	text=paste('Genes (N=',paste(ncol(profile),')',sep=''), sep=''),
	side=1,
	line=3.75,
	cex=1.5
	)
legend('topright',
    legend = categories,
    fill = color_pallete,
	bg='white',
	title='Host Population'
)

dev.off()

# plot marker gene coverage
# rows sorted the same as in heatmap
outname <- paste(c(cluster_id, 'marker_coverage', 'png'), collapse='.')
outpath <- paste(outdir, outname, sep='/')
png(outpath, width = 480, height = 480*3)
par(mar=c(6,6,2,5))
barplot(
	log(profile[rev(x$rowInd), 1],10),
	horiz = TRUE,
	xlab='',
	cex.axis=3, cex.lab=3.5,
	)
mtext(
	text='Log10 Marker Coverage',
	side=1,
	line=4,
	cex=3.5
)
dev.off()




