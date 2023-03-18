. "/c4/home/mwang13/miniconda3/etc/profile.d/conda.sh"
conda activate cnvkit

cnvkit.py target ../sgp_1_baits.bed --annotate ../../refFlat.txt --split -o my_targets.bed
cnvkit.py antitarget my_targets.bed -g ../../access-5k-mappable.hg19.bed -o my_antitargets.bed
cnvkit.py reference -o sgp_1_FlatReference.cnn -f ../../../../../../Database/References/ucsc_hg19.bwa-index/genome.fa -t my_targets.bed -a my_antitargets.bed
