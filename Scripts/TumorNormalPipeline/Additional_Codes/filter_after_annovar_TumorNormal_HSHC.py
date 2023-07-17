import re
import sys

Filein = sys.argv[1]
Fileout = sys.argv[2]

# remove common SNPs (VAF > 0.01 in human population)
def Filter_Common_SNP(myList):
    indexes = list(range(10,29))
    #indexes.append(30)
    PASS = True
    for i in indexes:
        if myList[i] != "." and float(myList[i]) > 0.01:
            PASS = False
    return(PASS)

# remove variants with mutant allele frequency < 0.02
def Filter_DP_MAF_1sample(ADs):
    Depth = ADs[0] + ADs[1]
    PASS = True
    MAF = 0
    if Depth < 10:
        PASS = False
    else:
        MAF = ADs[1]/Depth
        if MAF < 0.02:
            PASS = False
    return PASS, Depth, MAF

def Filter_DP_MAF_TN(Read_Info):
    PASS = True
    MAF = 0
    if Read_Info[3] < 10:
        PASS = False
    else:
        MAF = Read_Info[-1]/Read_Info[3]
        if MAF < 0.02:
            PASS = False
    return PASS, MAF

# process the final outputs of Bcftools and HaplotypeCaller, which both have only one sample
def Process_HSHC(Filein, Fileout):
    with open(Filein,"r") as file_H:
        with open(Fileout,"w") as fileout_H:
            for line in file_H:
                line = line.rstrip()
                Columns = line.split("\t")
                if "Func.refGene" in line:
                    Header = Columns[0:60] + ["Depth","Ref_reads","Alt_reads","MAF","VCF_Col_FILTER","VCF_Col_INFO","VCF_Col_FORMAT","VCF_Col_Details"]
                    fileout_H.write("\t".join(Header)+"\n")
                    continue
                # keep mutations in exonic regions, splicing sites and TERT promoter region
                if Columns[5] == "exonic" or Columns[5] == "splicing" or (Columns[5] == "upstream" and Columns[6] == "TERT"):
                    if "HC_Final" in Filein:
                        ADs = Columns[-1].rsplit(":")[1].rsplit(",")
                    else:
                        pattern = r"DP4=(\d+),(\d+),(\d+),(\d+)"
                        match = re.search(pattern,line)
                        ref = int(match.group(1)) + int(match.group(2))
                        alt = int(match.group(3)) + int(match.group(4))
                        ADs = [ref,alt]

                    ADs = [int(item) for item in ADs]
                    Depth = sum(ADs)

                    PASS1 = Filter_Common_SNP(Columns)
                    PASS2, Depth, MAF = Filter_DP_MAF_1sample(ADs)
                    if PASS1 == True and PASS2 == True:
                        AD_list = [Depth,ADs[0],ADs[1],MAF]
                        AD_list = [str(i) for i in AD_list]
                        newColumns = Columns[0:60] + AD_list + Columns[69:]

                        fileout_H.write("\t".join(newColumns)+"\n")

Process_HSHC(Filein, Fileout)
