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
    -I /path/to/${normal}_xx.bam \
    -I /path/to/${tumor}_xx.bam \
    -normal $normal_name \
    --germline-resource /path/to/af-only-gnomad.raw.sites.hg19.vcf.gz \
    --af-of-alleles-not-in-resource 0.001 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --panel-of-normals /path/to/PoN_WXS.vcf.gz \
    --f1r2-tar-gz /path/to/${tumor}.f1r2.tar.gz \
    -O /path/to/${tumor}.MT2.vcf.gz
    -bamout /path/to/${tumor}.MT2.bam

gatk LearnReadOrientationModel \
    -I /path/to/${tumor}.f1r2.tar.gz \
    -O /path/to/${tumor}.read-orientation-model.tar.gz

gatk FilterMutectCalls \
    -R /path/to/genome \
    -V /path/to/${tumor}.MT2.vcf.gz \
    --ob-priors /path/to/${tumor}.read-orientation-model.tar.gz \
    -O /path/to/${tumor}.MT2_filtered.vcf.gz

gatk LeftAlignAndTrimVariants \
    -R /path/to/genome \
    -V /path/to/${tumor}.MT2_filtered.vcf.gz \
    -O /path/to/${tumor}.MT2_Final.vcf \
    -no-trim true --split-multi-allelics true

table_annovar.pl \
    /path/to/${tumor}.MT2_Final.vcf \
    /path/to/Annovar/humandb --buildver hg19 \
    --vcfinput  --otherinfo  --thread 5 --remove \
    --operation g,f,f,f,f,f,f,f,f,f,f,f --protocol refGene,exac03,gnomad_exome,esp6500siv2_all,1000g2015aug_all,avsnp150,ucsf500normT,ucsf500normN,cosmic89,cbio2019jun,clinvar2019mar,ljb26_all \
    --outfile /path/to/${tumor}.MT2_Final.annovar

# B) FreeBayes
freebayes-parallel \
    ../Resources/FreeBayes_region.txt 5 \
    -f /path/to/genome \
    -b /path/to/${tumor}_xx.bam \
    -b /path/to/${normal}_xx.bam \
    -0 --genotype-qualities --min-repeat-entropy 1 \
    --min-alternate-fraction 0.05 -m 20 -q 20 \
    --pooled-discrete \
    --pooled-continuous \
    --report-genotype-likelihood-max \
    --allele-balance-priors-off \
    >/path/to/${tumor}.FB.vcf

    # reorder vcf to make sure that normal is in the last column
perl checkorder.pl /path/to/${tumor}.FB.vcf /path/to/${tumor}.FB.reorder.vcf $normal_name

python gmiMafAdder.py \
    /path/to/${tumor}.FB.reorder.vcf \
    /path/to/${tumor}.FB.reorder.MAF.vcf

vcffilter
    -f "QUAL > 20 & DP > 5" \
    -s /path/to/${tumor}.FB.reorder.MAF.vcf \
    |bcftools sort -O v -o /path/to/${tumor}.FB.reorder.MAF.Q20.vcf

python2 split_multiallelics.py
    /path/to/${tumor}.FB.reorder.MAF.Q20.vcf \
    /path/to/${tumor}.FB.reorder.MAF.Q20.multial.vcf

bcftools norm
    -f /path/to/genome \
    /path/to/${tumor}.FB.reorder.MAF.Q20.multial.vcf \
    |awk '{if (match($1,"#")) print; else {n=split($0,arr,"\t"); arr[4]=toupper(arr[4]); arr[5]=toupper(arr[5]); fix=arr[1]; for (i=2;i<=n;i++) fix=fix"\t"arr[i]; print fix}}' > /path/to/${tumor}.FB.reorder.MAF.Q20.multial.leftalign.vcf

python2 bcbio_freebayes_somatic_filter.py \
    /path/to/${tumor}.FB.reorder.MAF.Q20.multial.leftalign.vcf \
    /path/to/${tumor}.FB.reorder.MAF.Q20.multial.leftalign.bcbio.vcf

awk 'BEGIN{FS="\t"}{if($7!="REJECT")print}' ${tumor}.FB.reorder.MAF.Q20.multial.leftalign.bcbio.vcf > ${tumor}.FB_Final.vcf

#
table_annovar.pl \
    /path/to/${sample}.FB_Final.vcf \
    /path/to/Annovar/humandb --buildver hg19 \
    --vcfinput  --otherinfo  --thread 5 --remove \
    --operation g,f,f,f,f,f,f,f,f,f,f,f --protocol refGene,exac03,gnomad_exome,esp6500siv2_all,1000g2015aug_all,avsnp150,ucsf500normT,ucsf500normN,cosmic89,cbio2019jun,clinvar2019mar,ljb26_all \
    --outfile /path/to/${sample}.FB_Final.annovar

# C) Strelka2 (For somatic indels only, although SNVs are also called. Strelka2 website recommends using Manta first and then Strelka2 for indels, but from my test it did not seem to add much benefit)

	# Configuration
    # "--exome" flag is recommended if analyzing exome-seq, targeted or small gene panel sequencing data
configureStrelkaSomaticWorkflow.py \
    --config=./resources/configureStrelkaSomaticWorkflow.py.ini \
    --normalBam=/path/to/${normal}_xx.bam \
    --tumorBam=/path/to/${tumor}_xx.bam \
    --referenceFasta=/path/to/genome \
    --exome \
    --runDir=/path/to/analysis

    # Run mutation calling step
cd /path/to/analysis; runWorkflow.py -m local -j 6

# D) Bcftools (cancer hotspot locations only)
bcftools mpileup -O v -f $genome \
    -R ./resources/hotspots_coordinate_merged.bed \
    /path/to/${tumor}_xx.bam \
    -q 1 -Q 20 -d 100000 -L 100000 \
    | bcftools call - -c -A > /path/to/${tumor}.HS.vcf

bcftools norm -O v -o /path/to/${tumor}.HS_Final.vcf \
    -f $genome \
    /path/to/${tumor}.HS.vcf

table_annovar.pl \
    /path/to/${tumor}.HS_Final.vcf \
    /path/to/Annovar/humandb --buildver hg19 \
    --vcfinput  --otherinfo  --thread 5 --remove \
    --operation g,f,f,f,f,f,f,f,f,f,f,f --protocol refGene,exac03,gnomad_exome,esp6500siv2_all,1000g2015aug_all,avsnp150,ucsf500normT,ucsf500normN,cosmic89,cbio2019jun,clinvar2019mar,ljb26_all \
    --outfile /path/to/${tumor}.HS_Final.annovar

# E) HaplotypeCaller (Germline mutations)
gatk HaplotypeCaller
    -R /path/to/genome \
    -L /path/to/hg19_RefSeq_CDSunique_add10bp.bed \
    -I /path/to/${normal}_xx.bam \
    -O /path/to/${normal}.HT.vcf -ERC GVCF

gatk LeftAlignAndTrimVariants
    -O /path/to/${normal}.HT.vcf \
    -R /path/to/genome \
    -V /path/to/${normal}.HT_Final.vcf \
    -no-trim true --split-multi-allelics true

