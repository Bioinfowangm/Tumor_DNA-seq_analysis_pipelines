# version = 0.9.6

#Key variables(need to update for each patient)
REFERENCE_GENOME="genome.fa"  # Genome Reference (hg19)
dbSNP="common_all_20180418.chr.vcf" # SNP from database for FACETS
TUMOR_BAM_base=" "        # prefix of tumor bam sample(not including '.bam')
TUMOR_BAM=${TUMOR_BAM_base}.bam
GERMLINE_BAM_base=" "     # prefix of normal/germline sample(not including ".bam")
GERMLINE_BAM=${GERMLINE_BAM_base}.bam
Germline_vcf=" " # Germline mutations in tumor and normal samples, can be from freebayes or other sources
Purity=" "
Diploid_logR=" "
Target_Interval=" "
Patient_Name= " "
Gender="x/y" # 'x': female; 'y': male

#${*ALL_TUMOR_BAM*}: bam of all tumor samples
#${*ALL_FEMALE_GERMLINE_BAM}: bam of all germline female samples

# 1) CNVkit for copy number analysis
    # 1.1a) The normal/germline BAMs are used to generate a reference.cnn. If there is no normal/germline for a specific sequencing run, samples from other run (with the same capture strategy) can be used. Merging more normal/germline BAMs would be optimal. Importantly, pre-built reference.cnn file can also be applied for new samples.
    # Therefore, if the reference.cnn is not available, then run (all samples at once):
cnvkit.py batch \
    ${*ALL_TUMOR_BAM*} \
    --normal ${*ALL_FEMALE_GERMLINE_BAM} \
    --targets $Target_Interval \
    --annotate ../Related_Files/refFlat.txt \
    --fasta $REFERENCE_GENOME \
    --access ../Related_Files/access-5k-mappable.hg19.bed \
    --output-reference my_reference.cnn \
    --output-dir results/ \
    --diagram --scatter

    # 1.1b) Alternative, if the reference.cnn is already available, then run:
cnvkit.py batch \
    $TUMOR_BAM \
    -r my_reference.cnn \
    -d results/
cnvkit.py scatter \
    results/${TUMOR_BAM_base}.cnr \
    -s results/${TUMOR_BAM_base}.cns \
    -o results/${TUMOR_BAM_base}-scatter.pdf
cnvkit.py diagram \
    results/${TUMOR_BAM_base}.cnr \
    -s results/${TUMOR_BAM_base}.cns \
    -o results/${TUMOR_BAM_base}-diagram.pdf

    # 1.2) After running 1.1a or 1.1b, do re-centering (optional)
cnvkit.py call ${TUMOR_BAM_base}.cnr  --purity $Purity --center-at $Diploid_logR -m clonal -x $Gender -o ${TUMOR_BAM_base}.call.cnr
cnvkit.py call ${TUMOR_BAM_base}.cns  --purity $Purity --center-at $Diploid_logR -m clonal -x $Gender -v $Germline_vcf -o ${TUMOR_BAM_base}.call.cns

# 2) FACETS(https://github.com/mskcc/facets), which is an independent tool for copy number analysis
snp-pileup -g -q14 -Q20 -P100 -r15,0 $dbSNP ${Patient_Name}_facets.csv.gz $GERMLINE_BAM $TUMOR_BAM
Rscript ./Additional_Codes/run_facets.R $Patient_Name

# 3) Delly
delly call -g $REFERENCE_GENOME \
    -x ../Related_Files/human.hg19.excl.tsv \
    -o ${Patient_Name}.bcf \
    $TUMOR_BAM $GERMLINE_BAM

delly filter -f somatic -o ${Patient_Name}.somatic.bcf -s $Samp_List ${Patient_Name}.bcf

# 4) MISsensor
#Step1: Scan genome
msisensor scan $REFERENCE_GENOME -o hg19_msisensor.txt

#Step2: run for each patient
msisensor msi \
    -d hg19_msisensor.txt \
    -t $TUMOR_BAM \
    -n $GERMLINE_BAM \
    -e ../Related_Files/hg19_RefSeq_CDSunique_add10bp.bed \
    -o ${Patient_Name}.MSIsensor
