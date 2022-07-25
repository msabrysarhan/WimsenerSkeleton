#!/bin/bash
#
#SBATCH --job-name=contigs2bins
#SBATCH --mem=50000
#SBATCH --cpus-per-task=8
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --constraint=array-8core
#SBATCH --output=contigs2bins-%j.out
#SBATCH --error=contigs2bins-%j.err

#Just consider installing conda in your home dir, then install concoct using this:

#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda create -n concoct_env python=3 concoct


#Assembly
#The input file should be the fasta file resulted from the assembly

function display_help() {
    echo " "
    echo "contigs2bins: Convert assembled contigs into possible MAGs."
    echo " "
    echo " "
    echo "Usage: contigs2bins.sh -r reads.fastq.gz -c contigs.fasta -s sample_name -o output_directory"
		echo " "
    echo "   -r, --reads              Fastq reads used for assembly; Forward reads in case of PE assembly"
    echo "   -R, --Reads              Reverse reads in case of PE assembly"
    echo "   -c, --contigs            Fasta file containing cotigs resulted from assembly"
    echo "   -s, --sample_name         sample ID or name to be used as prefix"
		echo "   -o, --output_dir          path to output directory where genomes are to be stored"
    echo "   -h, --help                 Show this message."
    echo " "
    echo " "
    exit 1
}

# [ $# -eq 0 ] && { display_help ; exit 1; }

output_dir="contigs2bins_output"
Reads=0

while [ "$1" != "" ]; do
    case $1 in
        -r | --reads )      shift
                            reads=$1
                            ;;
      -R | --Reads )      shift
                          Reads=$1
                          ;;
        -c | --contigs )    shift
                            contigs=$1
                            ;;
		-s | --sample )     shift
                            sample=$1
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


#contigs=final.contigs.fasta
#sample="stone"
#reads=/proj/IcemanInstitute/Sabry/wimsener_cave/combined_assembly/stone_1.qc.fastq.gz

mkdir $output_dir
mkdir $output_dir/viruses
mkdir $output_dir/prokaryotes

module load seqkit

seqkit seq -m 1000 $contigs > $output_dir/final.contigs.fa


cd $output_dir

/home/user/sarhan/VirFinder.R final.contigs.fa

sed 's/"$//' final.contigs.viral.contigs.names | sed 's/^"//' > final.contigs.viral.names

rm final.contigs.viral.contigs.names



seqkit grep -n -f final.contigs.viral.names final.contigs.fa > viral.contigs.fasta

#Install checkV
#conda install -c conda-forge -c bioconda checkv
#checkv download_database ./

conda activate base

export CHECKVDB=/home/user/sarhan/checkv-db-v1.0

checkv end_to_end viral.contigs.fasta viral.qc -t 8

numb_comp_phages=$(expr $(cat viral.qc/complete_genomes.tsv | wc -l) - 1)


if [ $numb_comp_phages -gt 0 ]
then
  comp_phage_names=$(tail -$numb_comp_phages viral.qc/complete_genomes.tsv | cut -f1)
  for i in $comp_phage_names
  do
  seqkit grep -r -p $i final.contigs.fa > viruses/vir.$i.fasta
  echo $i >> comp_phage_names.txt
  seqkit grep -v -r -f comp_phage_names.txt final.contigs.fa > non.viral.contigs.fa
  rm comp_phage_names.txt

  done
else
   cat final.contigs.fa > non.viral.contigs.fa
fi



#seqkit grep -n -f $comp_phage_names final.contigs.fa > viral.contigs.fasta


#exclude the complete viral contigs from the original contigs
#exclude contigs shorter than 1000 nt



rm final.contigs.fa
rm comp_phage_names.txt
#rm viral.contigs.fasta
rm viral.contigs.fasta
rm final.contigs.viral.names

conda activate base
conda activate concoct_env

#Mapping and depth calculation

#conda deactivate

module load bowtie2
module load samtools

bowtie2-build non.viral.contigs.fa non.viral.contigs


if [ $Reads != 0 ]
then
  bowtie2  -x non.viral.contigs -1 $reads -2 $Reads -S $sample.non.viral.contigs.sam --no-unal --threads 8
else
  bowtie2  -x non.viral.contigs -U $reads -S $sample.non.viral.contigs.sam --no-unal --threads 8
fi

#bowtie2  -x non.viral.contigs -U $reads -S $sample.non.viral.contigs.sam --no-unal --threads 8

samtools view -@ 8 -Sbq 30 $sample.non.viral.contigs.sam > $sample.non.viral.contigs.bam
samtools sort -@ 8 $sample.non.viral.contigs.bam > $sample.non.viral.contigs.sorted.bam
samtools index -@ 8 $sample.non.viral.contigs.sorted.bam

rm $sample.non.viral.contigs.sam
rm $sample.non.viral.contigs.bam

rm non.viral.contigs.rev.2.bt2
rm non.viral.contigs.rev.1.bt2
rm non.viral.contigs.2.bt2
rm non.viral.contigs.1.bt2
rm non.viral.contigs.4.bt2
rm non.viral.contigs.3.bt2



#Binning
module load metabat
module load maxbin

jgi_summarize_bam_contig_depths $sample.non.viral.contigs.sorted.bam --outputDepth $sample.depth.txt
cat  $sample.depth.txt | cut -f1,3 > $sample.abundance.txt

mv non.viral.contigs.fa non.viral.contigs.fasta
gzip non.viral.contigs.fasta

metabat -i non.viral.contigs.fasta.gz -a $sample.depth.txt -o prokaryotes/$sample.metabat -m 1500 --saveCls
run_MaxBin.pl -contig non.viral.contigs.fasta.gz -abund $sample.abundance.txt -out prokaryotes/$sample.maxbin



for i in $(ls prokaryotes/*maxbin*.fasta | rev | cut -c 7- | rev)
do
mv $i.fasta $i.fa
done


cd prokaryotes
find . -type f ! -iname "*.fa" -delete
cd ..

gunzip non.viral.contigs.fasta.gz

#conda activate base
#conda activate concoct_env

cut_up_fasta.py non.viral.contigs.fasta -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fna

concoct_coverage_table.py contigs_10K.bed $sample.non.viral.contigs.sorted.bam > concoct.coverage_table.tsv


concoct --composition_file contigs_10K.fna --coverage_file concoct.coverage_table.tsv --threads 8 -b concoct_output/
merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
#mkdir concoct_output/fasta_bins

extract_fasta_bins.py non.viral.contigs.fasta concoct_output/clustering_merged.csv --output_path prokaryotes
#In case there is a problem with conda concoct, install this "conda install mkl"



cd prokaryotes


for i in $(ls --ignore={"*maxbin*","*metabat*"} | rev | cut -c 4- | rev)
do
mv $i.fa $sample.concoct.$i.fa
done

cd ..

#Binning purification
module load dastool

/apps/dastool/1.1.0/src/Fasta_to_Scaffolds2Bin.sh -i prokaryotes -e maxbin*.fa > maxbin.contigs2bins.tsv
/apps/dastool/1.1.0/src/Fasta_to_Scaffolds2Bin.sh -i prokaryotes -e metabat*.fa > metabat.contigs2bins.tsv
/apps/dastool/1.1.0/src/Fasta_to_Scaffolds2Bin.sh -i prokaryotes -e concoct*.fa > concoct.contigs2bins.tsv


DAS_Tool -i maxbin.contigs2bins.tsv,metabat.contigs2bins.tsv,concoct.contigs2bins.tsv -l maxbin,metabat,concoct -c non.viral.contigs.fasta -o $sample.dastool -t 8 --write_bins 1 --score_threshold 0.0




#Bins quality assessment
module load checkm

checkm lineage_wf $sample.dastool_DASTool_bins $sample.dastool_DASTool_bins -x ".fa" -t 8 --reduced_tree --pplacer_threads 8 --tab_table > $sample.dastool_DASTool_bins/$sample.checkm.summary

rm -r $sample.dastool_DASTool_bins/storage
rm -r $sample.dastool_DASTool_bins/bins
rm $sample.dastool_DASTool_bins/checkm.log


checkm lineage_wf prokaryotes prokaryotes -x ".fa" -t 8 --reduced_tree --pplacer_threads 8 --tab_table > prokaryotes/$sample.checkm.summary

rm -r prokaryotes/storage
rm -r prokaryotes/bins
rm prokaryotes/checkm.log

rm $sample.dastool_DASTool_scaffolds2bin.txt
rm $sample.dastool_concoct.eval
rm $sample.dastool_metabat.eval
rm $sample.dastool_maxbin.eval
rm $sample.dastool.seqlength
rm $sample.dastool_proteins.faa.archaea.scg
rm $sample.dastool_proteins.faa.bacteria.scg
rm $sample.dastool_proteins.faa
rm concoct.contigs2bins.tsv
rm metabat.contigs2bins.tsv
rm maxbin.contigs2bins.tsv
rm concoct.coverage_table.tsv
rm contigs_10K.fna
rm contigs_10K.bed
rm $sample.depth.txt
rm $sample.abundance.txt
