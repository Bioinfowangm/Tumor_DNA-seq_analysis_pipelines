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

module load CBI bwa/0.7.17 samtools/1.10
module load openjdk/1.8.0

# a) Run Fastqc
fastqc -t 2 -o $work_dir/$BCB/fastq/ \
    $work_dir/$BCB/TMP/${SL_ID_seq}_R1_001.fastq.gz \
    $work_dir/$BCB/TMP/${SL_ID_seq}_R2_001.fastq.gz

# b) Run fastp
./softwares/fastp -i $work_dir/$BCB/fastq/${SL_ID_seq}_R1_001.fastq.gz \
    -I $work_dir/$BCB/fastq/${SL_ID_seq}_R2_001.fastq.gz \
    -o $work_dir/TMP/${SL_ID_seq}_fastp_R1.fastq.gz \
    -O $work_dir/TMP/${SL_ID_seq}_fastp_R2.fastq.gz \
    --unpaired1 $work_dir/TMP/${SL_ID_seq}_fastp_unpaired.fastq.gz \
    --unpaired2 $work_dir/TMP/${SL_ID_seq}_fastp_unpaired.fastq.gz \
    -j $work_dir/TMP/${SL_ID_seq}_fastp.json \
    -h $work_dir/TMP/${SL_ID_seq}_fastp.html -c -w 5 

# c) Alignment
bwa mem \
    -t 5 $genome \
    $work_dir/TMP/${SL_ID_seq}_fastp_R1.fastq.gz \
    $work_dir/TMP/${SL_ID_seq}_fastp_R2.fastq.gz \
    -M -R "@RG\tID:${SL_ID_seq}\tPL:ILLUMINA\tPU:None\tLB:1\tSM:${Pat_ID}" | samtools sort -@ 1 - -o $work_dir/BAM/${Pat_ID}_${BCB}.bam

# d) Mark duplicates
./softwares/gatk-4.1.2.0/gatk MarkDuplicates -I $work_dir/BAM/${Pat_ID}_${BCB}.bam \
    -O $work_dir/BAM/${Pat_ID}_${BCB}_dedup.bam \
    -M $work_dir/BAM/${Pat_ID}_${BCB}.dup_metrics \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY SILENT \
    --REMOVE_DUPLICATES true \
    --TMP_DIR=$work_dir/TMP

# e) Check insert size
./softwares/gatk-4.1.2.0/gatk CollectInsertSizeMetrics \
    -I $work_dir/BAM/${Pat_ID}_${BCB}_dedup.bam \
    -O $work_dir/BAM/${Pat_ID}_${BCB}_insert_size_metrics.txt \
    -H $work_dir/BAM/${Pat_ID}_${BCB}_insert_size_histogram.pdf \
    --VALIDATION_STRINGENCY SILENT

# f) Run BQSR
./softwares/gatk-4.1.2.0/gatk BaseRecalibrator \
    -I $work_dir/BAM/${Pat_ID}_${BCB}_dedup.bam \
    -R $genome \
    --known-sites ./resources/BROAD_bundle/hapmap_3.3.hg19.sites.vcf \
    --known-sites ./resources/BROAD_bundle/1000G_phase1.indels.hg19.sites.vcf \
    --known-sites ./resources/BROAD_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    -O $work_dir/BAM/${Pat_ID}_${BCB}_BSQR.table

./softwares/gatk-4.1.2.0/gatk ApplyBQSR \
    -R $genome \
    -I $work_dir/BAM/${Pat_ID}_${BCB}_dedup.bam \
    --bqsr-recal-file $work_dir/BAM/${Pat_ID}_${BCB}_BSQR.table \
    -O $work_dir/BAM/${Pat_ID}_${BCB}_dedup_BSQR.bam

# g) Left alignment at indel site (GATK 3.6)
java -jar ./softwares/gatk_3.6.0/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R $genome \
    -known ./resources/BROAD_bundle/1000G_phase1.indels.hg19.sites.vcf \
    -known ./resources/BROAD_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    -I $work_dir/BAM/${Pat_ID}_${BCB}_dedup_BSQR.bam \
    -o $work_dir/BAM/${Pat_ID}_${BCB}_realignertargetcreator.intervals

java -Xmx30G -Djava.io.tmpdir=$work_dir/TMP -jar ./softwares/gatk_3.6.0/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R $genome \
    -targetIntervals $work_dir/BAM/${Pat_ID}_${BCB}_realignertargetcreator.intervals \
    -known ./resources/BROAD_bundle/1000G_phase1.indels.hg19.sites.vcf \
    -known ./resources/BROAD_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    -I $work_dir/BAM/${Pat_ID}_${BCB}_dedup_BSQR.bam \
    -o $work_dir/BAM/${Pat_ID}_${BCB}_preprocessed.bam

# h) Collect depth metrics
./softwares/gatk-4.1.2.0/gatk CollectHsMetrics \
    -I $work_dir/BAM/${Pat_ID}_${BCB}_preprocessed.bam \
    -O $work_dir/BAM/${Pat_ID}_${BCB}_Hs_Metrics.txt \
    -R $genome \
    --BAIT_INTERVALS ./resources/SmallGenePanel_bed/sgp_1_baits.interval_list \
    --TARGET_INTERVALS ./resources/SmallGenePanel_bed/sgp_1_targets.interval_list

# submit Step2
system("sbatch 2_Mutation_calling.sh $work_dir $BCB $genome $SL_ID_seq $Pat_ID");
# submit Step3
system("sbatch 3_CNV_SV.sh $work_dir $BCB $genome $SL_ID_seq $Pat_ID");
