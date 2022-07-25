#!/bin/bash
#
#SBATCH --job-name=MAGs_damage
#SBATCH --mem=50000
#SBATCH --cpus-per-task=8
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --constraint=array-8core
#SBATCH --output=MAGs_damage-%j.out
#SBATCH --error=MAGs_damage-%j.err


#Check for ancient DNA damage pattern

fasta_dir=$1
reads=$2
sample=$3
damage_results=$4


fasta_dir=/proj/IcemanInstitute/Sabry/wimsener_cave/combined_assembly/stone.megahit/stone.output/viruses
reads=/proj/IcemanInstitute/Sabry/wimsener_cave/combined_assembly/stone_1.qc.fastq.gz
sample=stone

cd $fasta_dir
mkdir damage_results



conda deactivate

module load bowtie2
module load samtools


for genome in $(ls *.fasta)
do
bowtie2-build $genome $genome
bowtie2  -x $genome -U $reads -S ${genome}_${sample}.sam --no-unal --threads 8
samtools view -@ 8 -Sbq 30 ${genome}_${sample}.sam > ${genome}_${sample}.bam
samtools sort -@ 8 ${genome}_${sample}.bam > ${genome}_${sample}.sorted.bam
samtools index -@ 8 ${genome}_${sample}.sorted.bam

java -jar /home/user/sarhan/DamageProfiler-1.1-java11.jar -i ${genome}_${sample}.sorted.bam -o ../damage_results/${genome}_${sample}_damage

rm ${genome}_${sample}.sam
rm ${genome}_${sample}.bam
rm ${genome}_${sample}.sorted.bam
rm ${genome}_${sample}.sorted.bam.bai
rm ${genome}.rev.2.bt2
rm ${genome}.rev.1.bt2
rm ${genome}.4.bt2
rm ${genome}.3.bt2
rm ${genome}.2.bt2
rm ${genome}.1.bt2

done

cd ../damage_results/


paste *damage/5pCtoT_freq.txt | tail -25 | cut -f$(seq -s, 2 2 $(ls * | wc -l))

paste *damage/3pGtoA_freq.txt | tail -25 | cut -f$(seq -s, 2 2 $(ls * | wc -l))


