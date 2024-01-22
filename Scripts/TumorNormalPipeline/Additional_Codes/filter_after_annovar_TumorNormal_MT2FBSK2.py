import re
import sys

TumorName = sys.argv[1]
NormalName = sys.argv[2]
Filein = sys.argv[3]
Vcfin = sys.argv[4]
Fileout = sys.argv[5]

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

# process the final outputs of Mutect2 and FreeBayes, which both have two samples (both tumor and normal)
def Process_Data(TumorName,NormalName, Filein, Vcfin, Fileout):
    VcfHeader = ""
    with open(Vcfin, "r") as vcf_H:
        for line in vcf_H:
            line = line.rstrip()
            if line.startswith("#CHROM"):
                VcfHeader = line.split("\t")
                break

    with open(Filein, "r") as file_H:
        with open(Fileout,"w") as fileout_H:
            for line in file_H:
                line = line.rstrip()
                Columns = line.split("\t")
                if "Func.refGene" in line:
                    Header = Columns[0:64]+ \
                    ["Normal_Depth","Normal_Ref_reads","Normal_Alt_reads","Tumor_Depth","Tumor_Ref_reads","Tumor_Alt_reads","Tumor_MAF"]+ \
                    ["VCF_Col_" + i for i in VcfHeader[6:11]]
                    fileout_H.write("\t".join(Header)+"\n")
                    continue

                if Columns[5] == "exonic" or Columns[5] == "splicing" or (Columns[5] == "upstream" and Columns[6] == "TERT"):
                    if "MT2_Final" in Filein:
                        if VcfHeader[-1] == NormalName:
                            Normal_Depth= Columns[-1].rsplit(":")[3]
                            NRef,NAlt = Columns[-1].rsplit(":")[1].rsplit(",")
                            Tumor_Depth= Columns[-2].rsplit(":")[3]
                            TRef,TAlt = Columns[-2].rsplit(":")[1].rsplit(",")
                        else:
                            Normal_Depth= Columns[-2].rsplit(":")[3]
                            NRef,NAlt = Columns[-2].rsplit(":")[1].rsplit(",")
                            Tumor_Depth= Columns[-1].rsplit(":")[3]
                            TRef,TAlt = Columns[-1].rsplit(":")[1].rsplit(",")
                    elif "FB_Final" in Filein:
                        if Columns[-1].rsplit(":")[-1]== "0" or Columns[-2].rsplit(":")[-1]== "0":
                            continue
                        if VcfHeader[-1] == NormalName:
                            Normal_Depth= Columns[-1].rsplit(":")[2]
                            NRef = Columns[-1].rsplit(":")[3].rsplit(",")[0]
                            NAlt = Columns[-1].rsplit(":")[-3]
                            Tumor_Depth= Columns[-2].rsplit(":")[2]
                            TRef = Columns[-2].rsplit(":")[3].rsplit(",")[0]
                            TAlt = Columns[-2].rsplit(":")[-3]
                        else:
                            Normal_Depth= Columns[-2].rsplit(":")[2]
                            NRef = Columns[-2].rsplit(":")[3].rsplit(",")[0]
                            NAlt = Columns[-2].rsplit(":")[-3]
                            Tumor_Depth= Columns[-1].rsplit(":")[2]
                            TRef = Columns[-1].rsplit(":")[3].rsplit(",")[0]
                            TAlt = Columns[-1].rsplit(":")[-3]

                    elif "SK2_Final" in Filein:
                        if Columns[-1].rsplit(":")[-1]== "0" or Columns[-2].rsplit(":")[-1]== "0":
                            continue
                        Normal_Depth= Columns[-2].rsplit(":")[1]
                        NRef = Columns[-2].rsplit(":")[3].rsplit(",")[0]
                        NAlt = Columns[-2].rsplit(":")[4].rsplit(",")[0]
                        Tumor_Depth= Columns[-1].rsplit(":")[1]
                        TRef = Columns[-1].rsplit(":")[4].rsplit(",")[0]
                        TAlt = Columns[-1].rsplit(":")[4].rsplit(",")[0]

                    Read_Info = [Normal_Depth,NRef,NAlt,Tumor_Depth,TRef,TAlt]
                    Read_Info_int = [int(x) for x in Read_Info]

                    PASS1 = Filter_Common_SNP(Columns)
                    PASS2, MAF = Filter_DP_MAF_TN(Read_Info_int)

                    if PASS1 == True and PASS2 == True:
                        newColumns = Columns[0:64] + Read_Info + [str(MAF)] + Columns[73:]
                        fileout_H.write("\t".join(newColumns)+"\n")

Process_Data(TumorName,NormalName, Filein, Vcfin, Fileout)
