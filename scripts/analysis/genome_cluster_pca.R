# libraries
library('RColorBrewer')

# functions
drop_const_col <- function(df){
    my_col  <- c()
    for (i in seq(1, ncol(df))){
        if (var(df[, i]) > 0) my_col <- c(my_col, i)
    }
    return(df[, my_col])
}

# command-line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
indir <- args[1] # Directory containing microbe_cnv output directories
sample_map <- args[2] # File mapping directory names to sample metadata
outpath <- args[3] # Path to output pdf
metric <- args[4]
group_by <- args[5]

# sample metadata
metadata <- read.csv(sample_map, sep='\t', header=T, row.names=1)

# list runs
runs <- list.files(indir)

# init matrix
inpath <- paste(c(indir, runs[1], 'genome_clusters.abundance'), collapse='/')
indata <- read.csv(inpath, sep='\t', header=T)
my_mat <- matrix(nrow=length(runs), ncol=nrow(indata))
colnames(my_mat) <- sort(indata$cluster_id)
rownames(my_mat) <- runs

# populate matrix
for (run in list.files(indir)){
	inpath <- paste(c(indir, run, 'genome_clusters.abundance'), collapse='/')
	indata <- read.csv(inpath, sep='\t', header=T, row.names=1)
	my_mat[run, ] <- indata[colnames(my_mat), metric]
}

# perform pca
my_pca  <- prcomp(drop_const_col(my_mat), retx=TRUE, center=TRUE, scale=TRUE)
prop_var1 <- as.character(round(summary(my_pca)$importance[2,1], 2))
prop_var2 <- as.character(round(summary(my_pca)$importance[2,2], 2))
prop_var3 <- as.character(round(summary(my_pca)$importance[2,3], 2))

# get point colors
categories <- sort(unique(metadata[,group_by]))
color_pallete <- brewer.pal(length(categories), 'Spectral')
colors <- c()
for (element in metadata[row.names(my_mat), group_by]){
	colors <- c(colors, color_pallete[which(categories == element)])
}

# open pdf and plot PCs
pdf(outpath)

par(mfrow=c(2,2), xpd=TRUE, mar=c(5.1, 4.1, 1, 1))

pt.cex <- 1.25

plot(my_pca$x[,'PC1'], my_pca$x[,'PC2'],
	col='black', bg=colors, pch=21, cex=pt.cex,
	xlab=paste(c('PC1 (', prop_var1, ')'), collapse=''), ylab=paste(c('PC2 (', prop_var2, ')'), collapse='')
)
plot(my_pca$x[,'PC1'], my_pca$x[,'PC3'],
	col='black', bg=colors, pch=21, cex=pt.cex,
	xlab=paste(c('PC1 (', prop_var1, ')'), collapse=''), ylab=paste(c('PC3 (', prop_var3, ')'), collapse='')
)
plot(my_pca$x[,'PC2'], my_pca$x[,'PC3'],
	col='black', bg=colors, pch=21, cex=pt.cex,
	xlab=paste(c('PC2 (', prop_var2, ')'), collapse=''), ylab=paste(c('PC3 (', prop_var3, ')'), collapse='')
)
plot(
	0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', bty='n'
)
legend('center',
	legend=categories,
	pt.bg=color_pallete, border='black',
	cex=pt.cex*1.5, pch=21, ncol=1, bty='n'
)
dev.off()
