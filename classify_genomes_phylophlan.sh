#!/bin/bash
#
#SBATCH --job-name=classify_genome
#SBATCH --mem=20000
#SBATCH --cpus-per-task=4
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --output=classify_genome-%j.out
#SBATCH --error=classify_genome-%j.err


function display_help() {
    echo " "
    echo "classify_genomes_phylophlan.sh: classify prokaryotic genomes based on PhyloPhlAn SGB.Jul20 database"
    echo " "
    echo "Usage: classify_genomes_phylophlan.sh -g genome_directory -x extension -o output_directory"
		echo " "
    echo "   -g, --genome_directory              Genomes-containing directory (Default ./phylophlan_(date & time))"
		echo "   -x, --extension                     fasta files extension (Default fa)"
		echo "   -o, --output_dir                    path to output directory where results are to be stored"
    echo "   -h, --help                          Show this message."
    echo " "
    echo " "
    exit 1
}

output_dir="phylophlan_"$(date +%Y%m%d_%H%M%S)
extension="fa"
genome_dir="."

while [ "$1" != "" ]; do
    case $1 in
        -g | --genome_dir )      shift
                            genome_dir=$1
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

ln -s /proj/IcemanInstitute/Sabry/phylophlan_databases/ phylophlan_databases

module load metaphlan
module load mash
module load ncbiblastplus

mkdir $genome_dir/$output_dir



phylophlan_metagenomic \
    -i $genome_dir \
    -o $genome_dir/$output_dir/phylophlan \
    --nproc 1 \
    -n 4 \
    -d SGB.Jul20 \
    -e fa 2>&1 | tee $genome_dir/$output_dir/phylophlan.log


while [ $(ls $genome_dir/$output_dir/phylophlan_dists | wc -l) -lt 51320 ]
do
phylophlan_metagenomic \
    -i $genome_dir \
    -o $genome_dir/$output_dir/phylophlan \
    --nproc 1 \
    -n 4 \
    -d SGB.Jul20 \
    -e fa 2>&1 | tee $genome_dir/$output_dir/phylophlan.log
done
