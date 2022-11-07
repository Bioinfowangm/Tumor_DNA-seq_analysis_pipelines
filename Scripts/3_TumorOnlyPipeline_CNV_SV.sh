#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=30G
#SBATCH --time=84:00:00
#SBATCH --output=%x-%A-%a.out

# read arguments
args=("$@")
work_dir=${args[0]};
BCB=${args[1]};
genome=${args[2]};
SL_ID_seq=${args[3]};
Pat_ID=${args[4]};

# CNVkit
. "/c4/home/mwang13/miniconda3/etc/profile.d/conda.sh"
conda activate ck2

cnvkit.py batch $work_dir/BAM/${Pat_ID}_${BCB}_preprocessed.bam -r ./resources/SmallGenePanel_bed/myflatref/sgp_1_FlatReference.cnn -d $work_dir/CNVkit
cnvkit.py scatter $work_dir/CNVkit/${Pat_ID}_${BCB}_preprocessed.cnr -s $work_dir/CNVkit/${Pat_ID}_${BCB}_preprocessed.cns -o $work_dir/CNVkit/${Pat_ID}_${BCB}_preprocessed-scatter.pdf
cnvkit.py diagram $work_dir/CNVkit/${Pat_ID}_${BCB}_preprocessed.cnr -s $work_dir/CNVkit/${Pat_ID}_${BCB}_preprocessed.cns -o $work_dir/CNVkit/${Pat_ID}_${BCB}_preprocessed-diagram.pdf

# Delly
./softwares/delly call -g $genome -o $work_dir/Delly/${Pat_ID}_${BCB}.delly.bcf $work_dir/BAM/${Pat_ID}_${BCB}_preprocessed.bam

# missensor2
msisensor2 msi -M ./resources/msisensor2/models_hg19_GRCh37 -t $work_dir/BAM/${Pat_ID}_${BCB}_preprocessed.bam -o $work_dir/msisensor2/${Pat_ID}_${BCB}
