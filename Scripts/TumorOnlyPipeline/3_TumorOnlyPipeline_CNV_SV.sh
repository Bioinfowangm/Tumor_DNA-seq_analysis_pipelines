#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=30G
#SBATCH --time=84:00:00
#SBATCH --output=%x-%A-%a.out

# CNVkit
. "/c4/home/mwang13/miniconda3/etc/profile.d/conda.sh"
conda activate ck2

cnvkit.py batch $work_dir/BAM/${sample}_preprocessed.bam -r ./resources/SmallGenePanel_bed/myflatref/sgp_1_FlatReference.cnn -d $work_dir/CNVkit
cnvkit.py scatter $work_dir/CNVkit/${sample}_preprocessed.cnr -s $work_dir/CNVkit/${sample}_preprocessed.cns -o $work_dir/CNVkit/${sample}_preprocessed-scatter.pdf
cnvkit.py diagram $work_dir/CNVkit/${sample}_preprocessed.cnr -s $work_dir/CNVkit/${sample}_preprocessed.cns -o $work_dir/CNVkit/${sample}_preprocessed-diagram.pdf
cnvkit.py export nexus-basic -o $work_dir/CNVkit/${sample}_preprocessed.nexus-basic $work_dir/CNVkit/${sample}_preprocessed.cnr

# Delly
./softwares/delly call -g $genome -o $work_dir/Delly/${sample}.delly.bcf $work_dir/BAM/${sample}_preprocessed.bam

# missensor2
msisensor2 msi -M ./resources/msisensor2/models_hg19_GRCh37 -t $work_dir/BAM/${sample}_preprocessed.bam -o $work_dir/msisensor2/${sample}
