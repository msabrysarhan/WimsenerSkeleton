#!/bin/bash
#
#SBATCH --job-name=mito
#SBATCH --mem=10000
#SBATCH --cpus-per-task=8
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --output=mito-%j.out
#SBATCH --error=mito-%j.err


function display_help() {
    echo " "
    echo "mitochondria.sh: Map the aDNA reads against human mitochondria genome,"
    echo "check the damage pattern, and generate vcf file for haplogroup assignment "
    echo " "
    echo "Usage: mitochondria.sh -r reads.fastq.gz -o output_directory"
		echo " "
    echo "   -r, --reads              Fastq reads used for assembly"
		echo "   -n, --name               sample name to be used as a prefix"
		echo "   -o, --output_dir          path to output directory where results are to be stored"
    echo "   -h, --help                 Show this message."
    echo " "
    echo " "
    exit 1
}

output_dir="mitochondria_"$(date +%Y%m%d_%H%M%S)
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
bowtie2  -x /scratch/sarhan/Refs/human/rCRS -U dedup.reads.fastq -S $name.mt.sam  --no-unal --threads 8

rm dedup.reads.fastq

module load samtools
samtools view -@ 8 -Sbq 30 $name.mt.sam > $name.mt.bam
samtools sort -@ 8 $name.mt.bam > $name.mt.sorted.bam
samtools index -@ 8 $name.mt.sorted.bam

rm $name.mt.sam
rm $name.mt.bam

conda deactivate
#Then, calculate the mapping statistics
unset DISPLAY
module load qualimap
qualimap bamqc -bam $name.mt.sorted.bam -outdir $name.mt.qualimap

#mapDamage and rescaling
module unload python2

module load mapdamage
mapDamage -i $name.mt.sorted.bam -r /scratch/sarhan/Refs/human/rCRS.fasta \
-d $name.mapDamage --rescale --rescale-out $name.mt.rescaled.bam
samtools index -@ 8 $name.mt.rescaled.bam


# create vcf file
module load bcftools
samtools mpileup -t AD -B -u -q 30 -f /scratch/sarhan/Refs/human/rCRS.fasta $name.mt.rescaled.bam | bcftools call -m --ploidy 1 -O v -o $name.mt.rescaled.vcf

#Schmutzi_check
module load schmutzi
mkdir $name.mt.schmutzi
contDeam.pl --length 10 --library double --out $name.mt.schmutzi/$name.mt $name.mt.sorted.bam
schmutzi.pl --notusepredC --uselength --ref /scratch/sarhan/Refs/human/rCRS.fasta $name.mt.schmutzi/$name.mt \
/apps/schmutzi/20171024/alleleFreqMT/eurasian/freqs/ $name.mt.sorted.bam
log2fasta -q 30 $name.mt.schmutzi/$name.mt_final_endo.log > $name.schmutzi.fasta
