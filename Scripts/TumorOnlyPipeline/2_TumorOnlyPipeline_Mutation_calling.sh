#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=30G
#SBATCH --time=84:00:00
#SBATCH --output=%x-%A-%a.out

module load CBI samtools

#$sample is the name of the sample to be processed here

#Tool1: Mutect2
./softwares/gatk-4.1.2.0/gatk --java-options "-Djava.io.tmpdir=/path/TMP" Mutect2 \
    -R $genome \
    -L ./resources/intervel_list_NoCentromereTelomere.bed \
    -I /path/to/${sample}_preprocessed.bam \ 
    --germline-resource ./resources/BROAD_bundle/af-only-gnomad.raw.sites.hg19.vcf.gz \
    --af-of-alleles-not-in-resource 0.0000025 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --panel-of-normals ./resources/Mutect2_panelofnormal/PoN_WXS.vcf.gz \
    --f1r2-tar-gz /path/to/${sample}.f1r2.tar.gz \
    -O /path/to/${sample}.somatic_m2.vcf.gz
##
./softwares/gatk-4.1.2.0/gatk LearnReadOrientationModel \
    -I /path/to/${sample}.f1r2.tar.gz \
    -O /path/to/${sample}.read-orientation-model.tar.gz
##
./softwares/gatk-4.1.2.0/gatk FilterMutectCalls \
    -V /path/to/${sample}.MT2.vcf.gz \
    -R $genome \
    --ob-priors /path/to/${sample}.read-orientation-model.tar.gz \
    -O /path/to/${sample}.MT2_filtered.vcf.gz
##
./softwares/gatk-4.1.2.0/gatk LeftAlignAndTrimVariants \
    -O /path/to/${sample}.MT2_Final.vcf \
    -R $genome \
    -V /path/to/${sample}.MT2_filtered.vcf.gz \
    -no-trim true --split-multi-allelics true
#
table_annovar.pl \
    /path/to/${sample}.MT2_Final.vcf \
    /path/to/Annovar/humandb --buildver hg19 \
    --vcfinput  --otherinfo  --thread 5 --remove \
    --operation g,f,f,f,f,f,f,f,f,f,f,f --protocol refGene,exac03,gnomad_exome,esp6500siv2_all,1000g2015aug_all,avsnp150,ucsf500normT,ucsf500normN,cosmic98,cbio2019jun,clinvar_20221231,ljb26_all \
    --outfile /path/to/${sample}.MT2_Final.annovar

##Tool2: FreeBayes
freebayes \
    -f $genome \
    /path/to/${sample}_preprocessed.bam \
    --min-repeat-entropy 1 \
    -q 20 \
    --pooled-discrete --pooled-continuous \
    --genotype-qualities \
    --report-genotype-likelihood-max \
    --allele-balance-priors-off \
    --min-alternate-fraction 0.02 >/path/to/${sample}.FB.vcf

python ./softwares/freebayes_scripts/gmiMafAdder.py \
    /path/to/${sample}.FB.vcf /path/to/${sample}.FB.MAF.vcf

vcffilter -f "QUAL > 20 & DP > 5" -s /path/to/${sample}.FB.MAF.vcf \
    |bcftools sort -O v -o /path/to/${sample}.FB.MAF.Q20.vcf

python2 ./softwares/freebayes_scripts/split_multiallelics.py \
    /path/to/${sample}.FB.MAF.Q20.vcf \
    /path/to/${sample}.FB.MAF.Q20.multial.vcf

bcftools norm -f $genome /path/to/${sample}.FB.MAF.Q20.multial.vcf \
    |awk '{if (match($1,"#")) print; else {n=split($0,arr,"\t"); arr[4]=toupper(arr[4]); arr[5]=toupper(arr[5]); fix=arr[1]; for (i=2;i<=n;i++) fix=fix"\t"arr[i]; print fix}}' > /path/to/${sample}.FB_Final.vcf

#
table_annovar.pl \
    /path/to/${sample}.FB_Final.vcf \
    /path/to/Annovar/humandb --buildver hg19 \
    --vcfinput  --otherinfo  --thread 5 --remove \
    --operation g,f,f,f,f,f,f,f,f,f,f,f --protocol refGene,exac03,gnomad_exome,esp6500siv2_all,1000g2015aug_all,avsnp150,ucsf500normT,ucsf500normN,cosmic98,cbio2019jun,clinvar_20221231,ljb26_all \
    --outfile /path/to/${sample}.FB_Final.annovar

##Tool3: Bcftools (cancer hotspot locations only)
bcftools mpileup -O v -f $genome \
    -R ./resources/hotspots_coordinate_merged.bed \
    /path/to/${sample}_preprocessed.bam \
    -q 1 -Q 20 -d 100000 -L 100000 \
    | bcftools call - -c -A > /path/to/${sample}.HS.vcf

bcftools norm -O v -o /path/to/${sample}.HS_Final.vcf \
    -f $genome \
    /path/to/${sample}.HS.vcf

table_annovar.pl \
    /path/to/${sample}.HS_Final.vcf \
    /path/to/Annovar/humandb --buildver hg19 \
    --vcfinput  --otherinfo  --thread 5 --remove \
    --operation g,f,f,f,f,f,f,f,f,f,f,f --protocol refGene,exac03,gnomad_exome,esp6500siv2_all,1000g2015aug_all,avsnp150,ucsf500normT,ucsf500normN,cosmic89,cbio2019jun,clinvar2019mar,ljb26_all \
    --outfile /path/to/${sample}.HS_Final.annovar

##Finally, filter the annovar outputs of the above three tools
#Three parameters are required: the ANNOVAR file (to be filtered), the filtered output, and the datasource(MT2, FB or HS(corresponding to hotspots))
python filter_after_annovar_TumorOnly.py /path/to/${sample}.MT2_Final.annovar.hg19_multianno.txt /path/to/${sample}.MT2_Final.annovar.hg19_multianno.Filtered.txt MT2
python filter_after_annovar_TumorOnly.py /path/to/${sample}.FB_Final.annovar.hg19_multianno.txt /path/to/${sample}.FB_Final.annovar.hg19_multianno.Filtered.txt FB 
python filter_after_annovar_TumorOnly.py /path/to/${sample}.FB_Final.annovar.hg19_multianno.txt /path/to/${sample}.FB_Final.annovar.hg19_multianno.Filtered.txt HS 
