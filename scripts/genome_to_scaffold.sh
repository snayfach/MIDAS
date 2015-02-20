cd /mnt/data/work/pollardlab/snayfach/projects/strain_variation/cnv_detection/develop_software/genome_clusters
i=0
for CLUSTER in `ls`; do
  zcat $CLUSTER/$CLUSTER.bed.gz | cut -f1,5 | sort -ur | gzip > $CLUSTER/$CLUSTER.genome_to_scaffold.gz
  i=$(($i + 1))
  echo $i
done