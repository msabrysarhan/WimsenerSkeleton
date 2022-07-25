#!/bin/bash
#
#SBATCH --job-name=classify_genome
#SBATCH --mem=10000
#SBATCH --cpus-per-task=8
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --output=classify_genome-%j.out
#SBATCH --error=classify_genome-%j.err


function display_help() {
    echo " "
    echo "classify_genomes_gtdbtk.sh: classify prokaryotic genomes"
    echo " "
    echo "Usage: classify_genomes_gtdbtk.sh -g genome_directory -x extension -o output_directory"
		echo " "
    echo "   -g, --genome_directory              Genomes-containing directory (Default ./gtdbtk_(date & time))"
		echo "   -x, --extension                     fasta files extension (Default fa)"
		echo "   -o, --output_dir                    path to output directory where results are to be stored"
    echo "   -h, --help                          Show this message."
    echo " "
    echo " "
    exit 1
}

output_dir="gtdbtk_"$(date +%Y%m%d_%H%M%S)
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

module load gtdbtk

gtdbtk classify_wf --genome_dir $genome_dir --out_dir $genome_dir/$output_dir -x $extension --force --cpus 8 --min_perc_aa 5
