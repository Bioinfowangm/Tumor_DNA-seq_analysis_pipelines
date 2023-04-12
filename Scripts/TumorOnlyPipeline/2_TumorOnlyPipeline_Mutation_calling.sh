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

module load CBI samtools

#Tool1: Mutect2
./softwares/gatk-4.1.2.0/gatk --java-options "-Djava.io.tmpdir=/path/TMP" Mutect2 \
    -R $genome \
    -L ./resources/intervel_list_NoCentromereTelomere.bed \
    -I /path/to/${Pat_ID}_${BCB}_preprocessed.bam \
    --germline-resource ./resources/BROAD_bundle/af-only-gnomad.raw.sites.hg19.vcf.gz \
    --af-of-alleles-not-in-resource 0.0000025 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --panel-of-normals ./resources/Mutect2_panelofnormal/PoN_WXS.vcf.gz \
    --f1r2-tar-gz /path/to/${Pat_ID}_${BCB}.f1r2.tar.gz \
    -O /path/to/${Pat_ID}_${BCB}.somatic_m2.vcf.gz
##
./softwares/gatk-4.1.2.0/gatk LearnReadOrientationModel \
    -I /path/to/${Pat_ID}_${BCB}.f1r2.tar.gz \
    -O /path/to/${Pat_ID}_${BCB}.read-orientation-model.tar.gz
##
./softwares/gatk-4.1.2.0/gatk FilterMutectCalls \
    -V /path/to/${Pat_ID}_${BCB}.somatic_m2.vcf.gz \
    -R $genome \
    --ob-priors /path/to/${Pat_ID}_${BCB}.read-orientation-model.tar.gz \
    -O /path/to/${Pat_ID}_${BCB}.somatic_m2_filtered.vcf.gz
##
./softwares/gatk-4.1.2.0/gatk LeftAlignAndTrimVariants \
    -O /path/to/${Pat_ID}_${BCB}.somatic_m2_filtered_SplitMulti.vcf \
    -R $genome \
    -V /path/to/${Pat_ID}_${BCB}.somatic_m2_filtered.vcf.gz \
    -no-trim true --split-multi-allelics true
#
table_annovar.pl \
    /path/to/${Pat_ID}_${BCB}.somatic_m2_filtered_SplitMulti.vcf \
    /path/to/Annovar/humandb --buildver hg19 \
    --vcfinput  --otherinfo  --thread 5 --remove \
    --operation g,f,f,f,f,f,f,f,f,f,f,f --protocol refGene,exac03,gnomad_exome,esp6500siv2_all,1000g2015aug_all,avsnp150,ucsf500normT,ucsf500normN,cosmic89,cbio2019jun,clinvar2019mar,ljb26_all \
    --outfile /path/to/${Pat_ID}_${BCB}.MT2.annovar


##Tool2: FreeBayes
freebayes \
    -f $genome \
    /path/to/${Pat_ID}_${BCB}_preprocessed.bam \
    --min-repeat-entropy 1 \
    -q 20 \
    --pooled-discrete --pooled-continuous \
    --genotype-qualities \
    --report-genotype-likelihood-max \
    --allele-balance-priors-off \
    --min-alternate-fraction 0.02 >/path/to/${Pat_ID}_${BCB}.FB.vcf

python ./softwares/freebayes_scripts/gmiMafAdder.py \
    /path/to/${Pat_ID}_${BCB}.FB.vcf /path/to/${Pat_ID}_${BCB}.FB.MAF.vcf

vcffilter -f "QUAL > 20 & DP > 5" -s /path/to/${Pat_ID}_${BCB}.FB.MAF.vcf \
    |bcftools sort -O v -o /path/to/${Pat_ID}_${BCB}.FB.MAF.Q20.vcf

python2 ./softwares/freebayes_scripts/split_multiallelics.py \
    /path/to/${Pat_ID}_${BCB}.FB.MAF.Q20.vcf \
    /path/to/${Pat_ID}_${BCB}.FB.MAF.Q20.multial.vcf

bcftools norm -f $genome /path/to/${Pat_ID}_${BCB}.FB.MAF.Q20.multial.vcf \
    |awk '{if (match($1,"#")) print; else {n=split($0,arr,"\t"); arr[4]=toupper(arr[4]); arr[5]=toupper(arr[5]); fix=arr[1]; for (i=2;i<=n;i++) fix=fix"\t"arr[i]; print fix}}' > /path/to/${Pat_ID}_${BCB}.FB.Final.vcf

#
table_annovar.pl \
    /path/to/${Pat_ID}_${BCB}.FB.Final.vcf \
    /path/to/Annovar/humandb --buildver hg19 \
    --vcfinput  --otherinfo  --thread 5 --remove \
    --operation g,f,f,f,f,f,f,f,f,f,f,f --protocol refGene,exac03,gnomad_exome,esp6500siv2_all,1000g2015aug_all,avsnp150,ucsf500normT,ucsf500normN,cosmic89,cbio2019jun,clinvar2019mar,ljb26_all \
    --outfile /path/to/${Pat_ID}_${BCB}.FB.annovar

##Tool3: Bcftools (cancer hotspot locations only)
bcftools mpileup -O v -f $genome \
    -R ./resources/hotspots_coordinate_merged.bed \
    /path/to/${Pat_ID}_${BCB}_preprocessed.bam \
    -q 1 -Q 20 -d 100000 -L 100000 \
    | bcftools call - -c -A > /path/to/${Pat_ID}_${BCB}.HS.vcf

bcftools norm -O v -o /path/to/${Pat_ID}_${BCB}.HS.Final.vcf \
    -f $genome \
    /path/to/${Pat_ID}_${BCB}.HS.vcf

table_annovar.pl \
    /path/to/${Pat_ID}_${BCB}.HS.vcf \
    /path/to/Annovar/humandb --buildver hg19 \
    --vcfinput  --otherinfo  --thread 5 --remove \
    --operation g,f,f,f,f,f,f,f,f,f,f,f --protocol refGene,exac03,gnomad_exome,esp6500siv2_all,1000g2015aug_all,avsnp150,ucsf500normT,ucsf500normN,cosmic89,cbio2019jun,clinvar2019mar,ljb26_all \
    --outfile /path/to/${Pat_ID}_${BCB}.HS.annovar
python filter_after_annovar.py ${Pat_ID}_${BCB}
