#!/bin/bash
#
#SBATCH --job-name=autosomal
#SBATCH --mem=10000
#SBATCH --cpus-per-task=8
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --output=autosomal-%j.out
#SBATCH --error=autosomal-%j.err


function display_help() {
    echo " "
    echo "autosomal.sh: Map the aDNA reads against human genome hg19 without mitochondria,"
    echo " "
    echo "Usage: autosomal.sh -r reads.fastq.gz -o output_directory"
		echo " "
    echo "   -r, --reads              Fastq reads used for assembly"
		echo "   -n, --name               sample name to be used as a prefix"
		echo "   -o, --output_dir          path to output directory where results are to be stored"
    echo "   -h, --help                 Show this message."
    echo " "
    echo " "
    exit 1
}

output_dir="autosomal_"$(date +%Y%m%d_%H%M%S)
name="sample"

while [ "$1" != "" ]; do
    case $1 in
        -r | --reads )      shift
                            reads=$1
                            ;;
				-n | --name )      shift
					                 name=$1
													 ;;
				-o | --output_dir ) shift
														output_dir=$1
														;;
        -h | --help )       display_help
                            exit
                            ;;
        * )                 display_help
                            exit 1
    esac
    shift
done


mkdir $output_dir


#First, remove duplicate reads
module load seqkit
seqkit rmdup -s -j 8 $reads | seqkit seq -m 25 > $output_dir/dedup.reads.fastq

cd $output_dir

#Next, map against mt genome

module load bowtie2
bowtie2  -x /scratch/sarhan/Refs/human/auto -U dedup.reads.fastq -S $name.auto.sam  --no-unal --threads 8

rm dedup.reads.fastq

module load samtools
samtools view -@ 8 -Sbq 30 $name.auto.sam > $name.auto.bam
samtools sort -@ 8 $name.auto.bam > $name.auto.sorted.bam
samtools index -@ 8 $name.auto.sorted.bam

rm $name.auto.sam
rm $name.auto.bam

#Then, calculate the mapping statistics
unset DISPLAY
module load qualimap
qualimap bamqc -bam $name.auto.sorted.bam -outdir $name.auto.qualimap

#mapDamage and rescaling
module unload python2

module load mapdamage
mapDamage -i $name.auto.sorted.bam -r /scratch/sarhan/Refs/human/auto.fasta \
-d $name.mapDamage --rescale --rescale-out $name.auto.rescaled.bam
samtools index -@ 8 $name.auto.rescaled.bam
