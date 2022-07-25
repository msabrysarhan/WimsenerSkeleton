#!/bin/bash
#
#SBATCH --job-name=metaphlan
#SBATCH --mem=10000
#SBATCH --cpus-per-task=4
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --output=metaphlan-%j.out
#SBATCH --error=metaphlan-%j.err


function display_help() {
    echo " "
    echo "metaphlan_full.sh: Generates metaphlan profiles for multiple inputs and returns a heatmap output"
    echo " "
    echo " "
    echo "Usage: metaphlan_full.sh -i reads.fastq.gz -o output_directory"
		echo " "
    echo "   -i, --reads              Fastq reads used for assembly, comma-separated multiple input"
		echo "   -o, --output_dir          path to output directory where heatmap is to be stored (Default metaphlan_output)"
    echo "   -n, --name                prefix for the heatmap without_extention (Default metaphlan_heatmap)"
    echo "   -h, --help                 Show this message."
    echo " "
    echo " "
    exit 1
}

output_dir="metaphlan_output_"$(date +%Y%m%d_%H%M%S)
name="metaphlan_heatmap"

while [ "$1" != "" ]; do
    case $1 in
        -i | --reads )      shift
                            reads=$1
                            ;;
        -n | --name )    shift
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

inputs=$(echo $reads|cut -d ',' --output-delimiter=$'\n' -f1-)

module load metaphlan

for file in $inputs
do
metaphlan --input_type fastq --nproc 4 --bowtie2out $output_dir/$file.bowtie2 $file > $output_dir/$file.metaphlan \
--add_viruses --read_min_len 25
done

rm $output_dir/*.bowtie2

merge_metaphlan_tables.py $output_dir/*.metaphlan > $output_dir/merged_abundance_table.txt

grep -E "s__|clade" $output_dir/merged_abundance_table.txt | sed 's/^.*s__//g' | cut -f1,3- | sed -e 's/clade_name/body_site/g' > $output_dir/merged_abundance_table_species.txt

python /apps/metaphlan/3.0.1/lib/python3.8/site-packages/metaphlan/utils/hclust2/hclust2.py \
-i $output_dir/merged_abundance_table_species.txt \
-o $output_dir/$name.pdf \
--f_dist_f euclidean \
--s_dist_f euclidean \
--cell_aspect_ratio 0.5 \
-l \
--flabel_size 10 \
--slabel_size 10 \
--max_flabel_len 100 \
--max_slabel_len 100 \
--minv 0.1 \
--dpi 300
