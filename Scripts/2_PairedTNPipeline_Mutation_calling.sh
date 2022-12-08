#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=30G
#SBATCH --time=84:00:00
#SBATCH --output=%x-%A-%a.out

# A) Mutect2
gatk GetSampleName \
    -I /path/to/germline_bam \
    -O /path/to/sample_name.txt

    # Get normal sample name (with the following shell command):
readarray -t SM</path/to/sample_name.txt; ctl_name=${SM[0]}

gatk Mutect2
    -R /path/to/genome \
    -L /path/to/intervel_list_NoCentromereTelomere.bed \
    -I /path/to/germline_bam \
    -I /path/to/tumor_bam \
    -normal $ctl_name \
    --germline-resource /path/to/af-only-gnomad.raw.sites.hg19.vcf.gz \
    --af-of-alleles-not-in-resource 0.001 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --panel-of-normals /path/to/PoN_WXS.vcf.gz \
    --f1r2-tar-gz /path/to/f1r2.tar.gz \
    -O /path/to/somatic_m2.vcf.gz
    -bamout /path/to/bamout

gatk LearnReadOrientationModel \
    -I /path/to/f1r2.tar.gz \
    -O /path/to/read-orientation-model.tar.gz

gatk FilterMutectCalls \
    -R /path/to/genome \
    -V /path/to/somatic_m2.vcf.gz \
    --ob-priors /path/to/read-orientation-model.tar.gz \
    -O /path/to/somatic_m2_filtered.vcf.gz

gatk LeftAlignAndTrimVariants \
    -R /path/to/genome \
    -V /path/to/somatic_m2_filtered.vcf.gz \
    -O /path/to/somatic_m2_filtered_SplitMulti.vcf \
    -no-trim true --split-multi-allelics true

# B) FreeBayes
freebayes-parallel \
    /path/to/region_file 5 \
    -f /path/to/genome \
    -b /path/to/tumor_bam \
    -b /path/to/normal_bam \
    -0 --genotype-qualities --min-repeat-entropy 1 \
    --min-alternate-fraction 0.05 -m 20 -q 20 \
    --pooled-discrete \
    --pooled-continuous \
    --report-genotype-likelihood-max \
    --allele-balance-priors-off \
    >/path/to/FB_vcf

    # reorder vcf to make sure that normal is in the last column
perl checkorder.pl /path/to/FB_vcf /path/to/FB_reorder_vcf

python gmiMafAdder.py \
    /path/to/FB_reorder_vcf \
    /path/to/FB_reorder_MAF_vcf

vcffilter
    -f "QUAL > 20 & DP > 5" \
    -s /path/to/FB_reorder_MAF_vcf \
    |bcftools sort -O v -o /path/to/FB_reorder_MAF_Q20_vcf

python2 split_multiallelics.py
    /path/to/FB_reorder_MAF_Q20_vcf \
    /path/to/FB_reorder_MAF_Q20_multial_vcf

bcftools norm
    -f /path/to/genome \
    /path/to/FB_reorder_MAF_Q20_multial_vcf \
    |awk '{if (match($1,"#")) print; else {n=split($0,arr,"\t"); arr[4]=toupper(arr[4]); arr[5]=toupper(arr[5]); fix=arr[1]; for (i=2;i<=n;i++) fix=fix"\t"arr[i]; print fix}}' > /path/to/FB_reorder_MAF_Q20_multial_norm_vcf


python2 bcbio_freebayes_somatic_filter.py \
    /path/to/FB_reorder_MAF_Q20_multial_norm_vcf \
    /path/to/FB_final_vcf

# C) Strelka2 (For somatic indels only, although SNVs are also called. Strelka2 website recommends using Manta first and then Strelka2 for indels, but from my test it did not seem to add much benefit)

	# Configuration
    # "--exome" flag is recommended if analyzing exome-seq, targeted or small gene panel sequencing data
configureStrelkaSomaticWorkflow.py \
    --config=./resources/configureStrelkaSomaticWorkflow.py.ini \
    --normalBam=/path/to/normal_bam \
    --tumorBam=/path/to/tumor_bam \
    --referenceFasta=/path/to/genome \
    --exome \
    --runDir=/path/to/analysis

    # Run mutation calling step
cd /path/to/analysis; runWorkflow.py -m local -j 6

# D) Bcftools (cancer hotspot locations only)
bcftools mpileup
    -O v -f $genome \
   	-R ./resources/hotspots_coordinate_merged.bed \
    /path/to/tumor_bam \
   	-q 1 -Q 20 -d 100000 -L 100000 \
   	| bcftools call - -c -A > /path/to/hotspot_HS_vcf

bcftools norm
    -O v -o /path/to/hotspot_HS_final_vcf \
    -f /path/to/genome \
    /path/to/hotspot_HS_vcf


# E) HaplotypeCaller (Germline mutations)
gatk HaplotypeCaller
    -R /path/to/genome \
    -L /path/to/hg19_RefSeq_CDSunique_add10bp.bed \
    -I /path/to/normal_bam \
    -O /path/to/normal_vcf -ERC GVCF

gatk LeftAlignAndTrimVariants
    -O /path/to/normal_SplitMulti_vcf \
    -R /path/to/genome \
    -V /path/to/normal_vcf \
    -no-trim true --split-multi-allelics true

