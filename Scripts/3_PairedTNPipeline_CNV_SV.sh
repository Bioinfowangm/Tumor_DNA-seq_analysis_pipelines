# version = 0.9.6

#Key variables(need to update for each patient)
REFERENCE_GENOME="genome.fa"  # Genome Reference (hg19)
TUMOR_BAM_base=" "        # prefix of tumor bam sample(not including '.bam')
TUMOR_BAM=${TUMOR_BAM_base}.bam
Germline_vcf=" " # Germline mutations in tumor and normal samples, can be from freebayes or other sources
Purity=" "
Diploid_logR=" "
Target_Interval=" "
Gender="x/y" # 'x': female; 'y': male

#${*ALL_TUMOR_BAM*}: bam of all tumor samples
#${*ALL_FEMALE_GERMLINE_BAM}: bam of all germline female samples

# CNVkit
# The normal/germline BAMs are used to generate a reference.cnn. If there is no normal/germline for a specific sequencing run, samples from other run (with the same capture strategy) can be used. Merging more normal/germline BAMs would be optimal. Importantly, pre-built reference.cnn file can also be applied for new samples.
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

# Re-center (optional)
cnvkit.py call ${TUMOR_BAM_base}.cnr  --purity $Purity --center-at $Diploid_logR -m clonal -x $Gender -o ${TUMOR_BAM_base}.call.cnr
cnvkit.py call ${TUMOR_BAM_base}.cns  --purity $Purity --center-at $Diploid_logR -m clonal -x $Gender -v $Germline_vcf -o ${TUMOR_BAM_base}.call.cns

# Delly
delly call -g $REFERENCE_GENOME \
    -x ../Related_Files/human.hg19.excl.tsv \
    -o ${Patient_Name}.bcf \
    $TUMOR_BAM $GERMLINE_BAM


delly filter -f somatic -o ${Patient_Name}.somatic.bcf -s $Samp_List ${Patient_Name}.bcf

# MISsensor
#Step1: Scan genome
msisensor scan $REFERENCE_GENOME -o hg19_msisensor.txt

#Step2: run for each patient
msisensor msi \
    -d hg19_msisensor.txt \
    -t $TUMOR_BAM \
    -n $GERMLINE_BAM \
    -e ../Related_Files/hg19_RefSeq_CDSunique_add10bp.bed \
    -o ${Patient_Name}.MSIsensor
