#!/bin/bash
#
#SBATCH --job-name=bins2damage
#SBATCH --mem=50000
#SBATCH --cpus-per-task=8
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --constraint=array-8core
#SBATCH --output=bins2damage-%j.out
#SBATCH --error=bins2damage-%j.err



function display_help() {
    echo " "
    echo "bins2damage.sh: check for ancient DNA damage in fasta genomes independently,"
    echo "using the tool DamageProfiler"
    echo "It returns tsv file containing matrix the C>T and G>A substitutions"
    echo "Usage: bins2damage.sh -g genome_directory -r reads.fastq(gz) -x extension -o output_directory"
		echo " "
    echo "   -g, --genome_directory              Genomes-containing directory"
    echo "   -r, --reads                         reads file to be used for mapping"
		echo "   -x, --extension                     fasta files extension (Default fa)"
    echo "   -s, --sample                        sample name"
		echo "   -o, --output_dir                    path to output directory where results are to be stored (Default ./DamageProfiler_(date & time))"
    echo "   -h, --help                          Show this message."
    echo " "
    echo " "
    exit 1
}


output_dir="DamageProfiler_"$(date +%Y%m%d_%H%M%S)
extension="fa"
genome_dir="."

while [ "$1" != "" ]; do
    case $1 in
        -g | --genome_dir )      shift
                            genome_dir=$1
                            ;;
        -r | --reads )      shift
            reads=$1
                        ;;
      -s | --sample )      shift
                      sample=$1
            ;;
				-x | --extension )      shift
					                 extension=$1
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

cd $genome_dir
mkdir $output_dir


module load bowtie2
module load samtools


for genome in $(ls *.$extension)
do
bowtie2-build $genome $genome
bowtie2  -x $genome -U $reads -S ${genome}_${sample}.sam --no-unal --threads 8
samtools view -@ 8 -Sbq 30 ${genome}_${sample}.sam > ${genome}_${sample}.bam
samtools sort -@ 8 ${genome}_${sample}.bam > ${genome}_${sample}.sorted.bam
samtools index -@ 8 ${genome}_${sample}.sorted.bam

java -jar /home/user/sarhan/DamageProfiler-1.1-java11.jar -i ${genome}_${sample}.sorted.bam -o $output_dir/${genome}_${sample}_damage

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

cd $output_dir
sample=bone
ls ../*.fa | tr "\n" "\t" > $sample.damage.txt
echo -e "\n" >> $sample.damage.txt
paste */5pCtoT_freq.txt | tail -25 | cut -f1,$(seq -s, 2 2 $(ls * | wc -l))| sort -n | cut -f2- >> $sample.damage.txt
paste */3pGtoA_freq.txt | tail -25 | cut -f1,$(seq -s, 2 2 $(ls * | wc -l))| sort -n -r | cut -f2- >> $sample.damage.txt
