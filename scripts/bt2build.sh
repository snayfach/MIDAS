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
#$ -tc 200
#$ -t 1-4202

# create tempdir on scratch
TMPDIR=/scratch
SCRATCH_DIR=`mktemp -d`
cd $SCRATCH_DIR

# select genome cluster
PROJECT=/pollard/shattuck0/snayfach/projects/strain_variation
CLUSTERS=${PROJECT}/genome_clustering/clusters/top_20/0.035/genome_clusters.txt
CLUSTER_ID=`sed -n ${SGE_TASK_ID}p ${CLUSTERS} | cut -f1`

# check if output already exists
if [ -d ${PROJECT}/read_mapping/bt2dbs/${CLUSTER_ID} ]; then exit; fi

# copy genomes locally
GENOMES=`cut -f2 ${PROJECT}/pan_genomics/nr_genomes/0.05/${CLUSTER_ID}/representatives.txt`
PATRIC=/pollard/shattuck0/snayfach/genomes/PATRIC/ftp.patricbrc.org/patric2/genomes
ID_TO_DIR=${PROJECT}/genome_clustering/db/dir_to_genome_id.txt
for GENOME in $GENOMES; do
DIR=`grep $'\t'${GENOME}$ $ID_TO_DIR | cut -f1`
FNA=${PATRIC}/${DIR}/${DIR}.fna
concurrent cat $FNA >> ${CLUSTER_ID}.fna
done

# build btdb
mkdir ${CLUSTER_ID}
BOWTIE2_BUILD=${PROJECT}/shotgun_isolates/lib/bowtie2-2.2.4/bowtie2-build
$BOWTIE2_BUILD ${CLUSTER_ID}.fna ${CLUSTER_ID}/${CLUSTER_ID} >& /dev/null

# manage out
DEST=${PROJECT}/read_mapping/bt2dbs
concurrent cp -r ${CLUSTER_ID} ${DEST}

# cleanup
rm -r $SCRATCH_DIR


