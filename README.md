


### Examples

export PATH=$PATH:$WORK/projects/strain_variation/cnv_detection/develop_software

microbe_cnv.py \
-1 example/SRR413772_1.fastq.gz \
-2 example/SRR413772_2.fastq.gz \
-D genome_clusters \
-o example \
--reads 1000 \
--abun 0.09 \
--map \
--verbose