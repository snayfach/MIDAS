#!/bin/bash
#
#$ -S /bin/bash
#$ -l arch=linux-x64 
#$ -l scratch=8G
#$ -l mem_free=8G
#$ -cwd
#$ -l h_rt=24:00:00
#$ -o /dev/null
#$ -e /dev/null
#$ -r y
#$ -t 1-60

# create tempdir on scratch
export TMPDIR=/scratch
SCRATCH_DIR=`mktemp -d`
cd $SCRATCH_DIR

# get reference paths
PROJECT=/pollard/shattuck0/snayfach/projects/strain_variation
BOWTIE2_BUILD=${PROJECT}/cnv_detection2/microbe_cnv/lib/bowtie2-2.2.4/bowtie2-build
PATRIC=/pollard/shattuck0/snayfach/genomes/PATRIC3/ftp.patricbrc.org/patric2/patric3/genomes

# loop over indexes
BATCH_SIZE=100
START=$(($(($BATCH_SIZE * $SGE_TASK_ID)) - $BATCH_SIZE))
for SUB_ID in `seq $BATCH_SIZE`; do

# get index
INDEX=$(($START + $SUB_ID))

# select genome cluster
CLUSTERS=${PROJECT}/genome_clustering2/clusters_combined/top_30/0.035/genome_clusters.txt
CLUSTER_ID=`sed -n ${INDEX}p ${CLUSTERS} | cut -f1`

# check whether output already exists
if [ -d ${PROJECT}/cnv_detection2/microbe_cnv/genome_clusters/${CLUSTER_ID} ]; then continue; fi

# cat genome sequences
GENOMES=`cut -f2 ${PROJECT}/pan_genomics2/nr_genomes/${CLUSTER_ID}/representatives.txt`
for GENOME in $GENOMES; do
FNA=${PATRIC}/${GENOME}/${GENOME}.fna
cat $FNA >> ${CLUSTER_ID}.fna
done

# build btdb for population of genomes
mkdir -p ${CLUSTER_ID}/btdb
$BOWTIE2_BUILD ${CLUSTER_ID}.fna ${CLUSTER_ID}/btdb/${CLUSTER_ID} >& /dev/null

# get path to representative genome
CENTROIDS=${PROJECT}/genome_clustering2/clusters_combined/top_30/0.035/cluster_to_centroid.txt
GENOME=`grep ^$CLUSTER_ID$'\t' $CENTROIDS | cut -f2`
PATRIC=/pollard/shattuck0/snayfach/genomes/PATRIC3/ftp.patricbrc.org/patric2/patric3/genomes

# build btdb for representative genome
mkdir -p ${CLUSTER_ID}/btdb_rep
$BOWTIE2_BUILD ${PATRIC}/${GENOME}/${GENOME}.fna ${CLUSTER_ID}/btdb_rep/${CLUSTER_ID} >& /dev/null

# manage out
DEST=${PROJECT}/cnv_detection2/microbe_cnv/genome_clusters
concurrent cp -r ${CLUSTER_ID} ${DEST}
rm -r ${CLUSTER_ID} ${CLUSTER_ID}.fna

done

# cleanup
rm -r $SCRATCH_DIR


