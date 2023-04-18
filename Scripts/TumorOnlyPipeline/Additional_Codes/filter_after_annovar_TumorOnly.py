import re
import sys

#Recommended suffix of input name: "_Final.annovar.hg19_multianno.txt"
#Recommended suffix of output name: "_Final.annovar.hg19_multianno_Filtered.txt"
#Source = "MT2" or "FB", or "HS"
Filein = sys.argv[1]
Fileout = sys.argv[2]
Source = sys.argv[3]

# remove common SNPs (VAF > 0.01 in human population)
def Filter_Common_SNP(myList):
    indexes = list(range(10,29))
    indexes.append(30)
    PASS = True
    for i in indexes:
        if myList[i] != "." and float(myList[i]) > 0.01:
            PASS = False
    return(PASS)

# remove variants with mutant allele frequency < 0.02
def Filter_DP_MAF(ADs):
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


def Process_files(Filein, Fileout, Source):

    with open(Filein,"r") as file_H:
        with open(Fileout,"w") as fileout_H:
            for line in file_H:
                line = line.rstrip()
                Columns = line.split("\t")
                if "Func.refGene" in line:
                    Header = Columns[0:59] + ["Depth","Ref_reads","Alt_reads","MAF","VCF_Col_FILTER","VCF_Col_INFO","VCF_Col_FORMAT","VCF_Col_Sample"]
                    fileout_H.write("\t".join(Header)+"\n")
                    continue
                # keep mutations in exonic regions, splicing sites and TERT promoter region
                if Columns[5] == "exonic" or Columns[5] == "splicing" or (Columns[5] == "upstream" and Columns[6] == "TERT"):
                    if Source == "MT2":
                        ADs = Columns[-1].rsplit(":")[1].rsplit(",")
                    elif Source == "FB":
                        ADs = Columns[-1].rsplit(":")[3].rsplit(",")
                    else:
                        pattern = r"DP4=(\d+),(\d+),(\d+),(\d+)"
                        match = re.search(pattern,line)
                        ref = int(match.group(1)) + int(match.group(2))
                        alt = int(match.group(3)) + int(match.group(4))
                        ADs = [ref,alt]

                    ADs = [int(item) for item in ADs]
                    Depth = sum(ADs)

                    PASS1 = Filter_Common_SNP(Columns)
                    PASS2, Depth, MAF = Filter_DP_MAF(ADs)
                    if PASS1 == True and PASS2 == True:
                        AD_list = [Depth,ADs[0],ADs[1],MAF]
                        AD_list = [str(i) for i in AD_list]
                        newColumns = Columns[0:59] + AD_list + Columns[69:73]

                        fileout_H.write("\t".join(newColumns)+"\n")

Process_files(Filein,Fileout,Source)
