#!/usr/bin/env python
"""gmiMafAdder.py
Usage: gmiMafAdder.py [in.vcf] [out.vcf]
Adds mutant allele read counts, mutant allele frequencies, and locus coverage info to VCFs.
For now, only run this script on SNP_indel VCFs. Don't run on SV vcfs (pindel_sv, Delly).
"""

from __future__ import print_function
import sys

print("Running Python version", sys.version_info, file=sys.stderr)
print("Arguments given are:", *sys.argv, file=sys.stderr)

# Deal with argv
if len(sys.argv) == 3:
    if sys.argv[1] == "stdin":
        inFile = sys.stdin
    else:
        inFile = open(sys.argv[1], 'r')
    if sys.argv[2] == "stdout":
        outFile = sys.stdout
    else:
        outFile = open(sys.argv[2], 'w')
else:
    sys.exit(__doc__)


CMD = "##gmiMafAdder_command=" + " ".join(sys.argv)

ADHEAD = (
    '##INFO=<ID=ADERROR,Number=0,Type=Flag,Description="The AD field '
    'from the original variant caller VCF has an error at this locus. '
    'Please use caution when interpreting AD, GMIMUT, GMIMAF, and GMICOV.">'
    )

HEADER = (
    '##INFO=<ID=FMTERROR,Number=0,Type=Flag,Description="The FORMAT field '
    'from the original variant caller VCF has an error at this locus. '
    'Please use caution when interpreting sample FORMAT values.">\n'
    '##FORMAT=<ID=GMIMUT,Number=A,Type=Integer,Description="GMI-derived '
    'mutant read count: the number of reads observed at the variant locus '
    'that support each ALT allele.">\n'
    '##FORMAT=<ID=GMIMAF,Number=A,Type=Integer,Description="GMI-derived MAF: '
    'the percent (0-100) of reads observed at the variant locus '
    'that support each ALT allele.">\n'
    '##FORMAT=<ID=GMICOV,Number=1,Type=Integer,Description="GMI-derived '
    'coverage: total read depth at the variant locus.">'
    )


def convert_to_ints(arr):
    """Convert all elements in arr to ints."""
    converted = []
    for item in arr:
        if item == '.':
            converted.append(0)
        else:
            converted.append(int(item))
    return converted


hasAD = False
for line in iter(inFile):
    line = line.rstrip('\n\r')

    if line.startswith("##"):
        print(line, file=outFile)
        if line.startswith("##FORMAT=<ID=AD,"):
            hasAD = True
        continue
    elif line.startswith("#"):
        print(CMD, file=outFile)
        if hasAD:
            print(ADHEAD, file=outFile)
        print(HEADER, file=outFile)
        print(line, file=outFile)
        continue

    fields = [str(x) for x in line.split('\t')]
    print("Looking at variant:", *fields, file=sys.stderr)
    alts = fields[4].split(",")
    fmt = fields[8].split(":")

    printFields = list(fields)
    printFields[8] += ":GMIMUT:GMIMAF:GMICOV"

    for j in range(9, len(fields)):
        sampleValues = ['.'] * len(fmt)
        if fields[j] != '.':
            sampleValues = fields[j].split(':')

        # Found bug in haplotype caller where a variant did not have
        # sample values for all FORMAT fields.
        while len(sampleValues) < len(fmt):
            sampleValues.append('.')
            # Mutect2 intentionally does this for 2 fields.
            if (len(sampleValues) == len(fmt)-2
                and fmt[-2:] == ['SA_MAP_AF', 'SA_POST_PROB']):
                continue
            print("FMTERROR: Wrong number of FMT values",
                  "found for this variant.", file=sys.stderr)
            if not printFields[7].endswith(';FMTERROR'):
                printFields[7] += ';FMTERROR'

        if all([x in fmt for x in ("AO", "RO")]):
            # This is a FreeBayes VCF.
            # For ".:." cases, MAF=0 and cov=0
            # Sometimes the sample column is just '.' instead of '.:.:*'

            ro = sampleValues[fmt.index('RO')]
            print("Column:", j, ", ao:", sampleValues[fmt.index('AO')],
                  ", ro:", ro, file=sys.stderr)
            # make sure ro and ao are numeric
            if ro == '.':
                ro = 0
            ao = convert_to_ints(sampleValues[fmt.index('AO')].split(','))
            # sometimes ao is '.' but there is more than one alt allele
            while len(ao) < len(alts):
                ao.append(0)

            rawmut = ao
            rawcov = int(ro) + sum(ao)

        elif "AD" in fmt:
            # this is an HC, Mutect, Pindel, SID, or UGT vcf.
            # for Pindel long insertions, AD only contains info about the ALT
            # allele. For now, don't run this script on Pindel SVs.
            # ad[0] is reference, everything else is alt
            print("Column:", j, ", ad:", sampleValues[fmt.index('AD')],
                  file=sys.stderr)
            ad = convert_to_ints(sampleValues[fmt.index('AD')].split(','))

            # found errors in UGT vcfs
            while len(ad) < len(alts) + 1:
                if not printFields[7].endswith('ADERROR'):
                    print("ADERROR: wrong number of elements in AD field.",
                          file=sys.stderr)
                    printFields[7] += ';ADERROR'
                ad.append(int(0))
                sampleValues[fmt.index('AD')] += ',0'

            rawmut = ad[1:]
            rawcov = sum(ad)

        elif hasAD is True and "AD" not in fmt:
            # found bug with haplotype caller where some variants had no AD
            # information.
            print("ADERROR: Input VCF has AD field,"
                  "but no AD found for this variant.", file=sys.stderr)
            if not printFields[7].endswith('ADERROR'):
                printFields[7] += ';ADERROR'

            rawmut = ['.'] * len(alts)
            rawcov = '.'

        # final values
        gmimut = ','.join([str(x) for x in rawmut])
        gmicov = str(rawcov)
        if rawcov == 0 or rawcov == '.':
            gmimaf = ','.join([str(rawcov)] * len(alts))
        else:
            gmimaf = ','.join(["%i" % (100.0 * int(val) / rawcov + 0.5)
                               for val in rawmut])
        newSampleValues = sampleValues + [gmimut, gmimaf, gmicov]
        printFields[j] = ":".join(newSampleValues)

    print(*printFields, sep="\t", file=outFile)
