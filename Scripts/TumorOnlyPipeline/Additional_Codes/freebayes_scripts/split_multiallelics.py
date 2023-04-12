"""split_multiallelics.py
Usage: python split_multiallelics.py [in.vcf] [out.vcf]
"""

import sys

# Deal with argv
if len (sys.argv) == 3:
    if sys.argv[1] == "stdin":
        inFile = sys.stdin
    else:
        inFile = open(sys.argv[1],'r')
    if sys.argv[2] == "stdout":
        outFile = sys.stdout
    else:
        outFile = open(sys.argv[2],'w')
else:
    print >> sys.stderr, \
        __doc__,
    sys.exit (1)

print >> sys.stderr, 'Arguments given are: %s' % sys.argv

# store number types for all INFO and FORMAT fields
idToNumberType = {'info': {},
                  'format': {}}
idToDescription = {'info': {},
                   'format': {}}

for line in iter(inFile):
    line = line.rstrip('\n\r')
    #print >> sys.stderr, "Looking at line %s" % line
    if line.startswith("##"):
        print >> outFile, line
        if line.startswith("##INFO") or line.startswith("##FORMAT"):
            data = line[line.index('<')+1:-1].split(',')
            id = data[0].split('=')[1]
            number = data[1].split('=')[1]
            description = data[3].split('=')[1]
            fieldType = line[2:line.index('=')].lower()
            idToNumberType[fieldType][id] = number
            idToDescription[fieldType][id] = description
    elif line.startswith("#"):
        print >> outFile, '##split_multiallelics_command=%s' % " ".join(map(str,sys.argv))
        print >> outFile, '##INFO=<ID=SPLITMULTIALLELIC,Number=0,Type=Flag,Description="Variant represents one alternate allele from a multiallelic variant call.">'
        print >> outFile, line
        columns = line.split('\t')
        numSamples = len(columns) - 9
        print >> sys.stderr, "There are %s samples" % numSamples
    else:
        columns = map(str,line.split('\t'))
        alt = columns[4].split(',')
        if len(alt) == 1:
            print >> outFile, line
        else:
            print >> sys.stderr, "Found multiallelic variant: %s" % columns
            # read in info and format columns.
            infoKeys = []
            infoValues = []
            for item in columns[7].split(';'):
                if '=' in item:
                    key, value = item.split('=')
                else:
                    key = item
                    value = 'FLAGITEM'
                infoKeys.append(key)
                infoValues.append(map(str,value.split(',')))

            formatKeys = columns[8].split(':')
            formatValues = []
            for sampleIndex in range(0,numSamples):
                formatValues.append(map(str,columns[9+sampleIndex].split(':')))

            for altIndex in range(0,len(alt)):
                # indices for Number = G columns
                altNumber = altIndex+1
                homRefIndex = 0
                hetIndex = altNumber*(altNumber+1)/2
                homAltIndex = altNumber*(altNumber+1)/2 + altNumber

                # parse genotypes here first, because INFO/AN needs this info.
                # GT does not always have to be in the format field, but if it is, it should always be first according to VCF specification.
                sampleGenotypes = []
                if formatKeys[0] == 'GT':
                    # AN just counts the total number of non-'.' alleles in the called genotypes for the VCF line.
                    newAN = 0
                    for sampleIndex in range(0,numSamples):
                        genotype = map(str,formatValues[sampleIndex][0].split('/'))
                        if genotype == map(str,[altNumber,altNumber]):
                            sampleGenotypes.append('1/1')
                            newAN += 2
                        elif genotype == map(str,[0,altNumber]):
                            sampleGenotypes.append('0/1')
                            newAN += 2
                        elif genotype == map(str,[0,0]):
                            sampleGenotypes.append('0/0')
                            newAN += 2
                        elif str(0) in genotype:
                            # is this valid?
                            sampleGenotypes.append('0/.')
                            newAN += 1
                        elif str(altNumber) in genotype:
                            sampleGenotypes.append('./1')
                            newAN += 1
                        else:
                            sampleGenotypes.append('./.')

                # create new split record to print.
                # start with everything up to INFO column
                splitRecord = '\t'.join(map(str,columns[0:4])) + '\t' + alt[altIndex] + '\t' + '\t'.join(map(str,columns[5:7])) + '\t'

                # deal with INFO column
                for i in range(0,len(infoKeys)):
                    # deal with INFO column
                    splitRecord += infoKeys[i]
                    if idToNumberType['info'][infoKeys[i]] == 'A':
                        splitRecord += '=' + infoValues[i][altIndex] + ';'
                    elif idToNumberType['info'][infoKeys[i]] == '.' and len(infoValues[i]) == len(alt) + 1 and ('ref' in idToDescription['info'][infoKeys[i]] or 'Ref' in idToDescription['info'][infoKeys[i]]):
                        # this is a list with 1+len(alt) elements; the first is for REF and the rest are for the ALTs, in order.
                        splitRecord += '=' + infoValues[i][0] + ',' + infoValues[i][altIndex] + ';'
                    elif idToNumberType['info'][infoKeys[i]] == 'G':
                        splitRecord += '=' + ','.join(map(str,[infoValues[i][homRefIndex],infoValues[i][hetIndex],infoValues[i][homAltIndex]])) + ';'
                    elif idToNumberType['info'][infoKeys[i]] == '0':
                        splitRecord += ';'
                    elif infoKeys[i] == 'AN':
                        splitRecord += '=' + str(newAN) + ';'
                    else:
                        splitRecord += '=' + infoValues[i][0] + ';'
                splitRecord += "SPLITMULTIALLELIC"

                # add FORMAT column
                splitRecord += '\t' + columns[8]

                # sample-specific columns
                for sampleIndex in range(0,numSamples):
                    splitRecord += '\t'
                    for j in range(0,len(formatKeys)):
                        vals = map(str,formatValues[sampleIndex][j].split(','))
                        #print >> sys.stderr, "FORMAT/%s : %s" % (formatKeys[j], vals)
                        if formatKeys[j] == 'GT':
                            splitRecord += sampleGenotypes[sampleIndex] + ':'
                        elif vals == ["."]:
                            splitRecord += ".:"
                        elif idToNumberType['format'][formatKeys[j]] == 'A':
                            splitRecord += vals[altIndex] + ':'
                        elif idToNumberType['format'][formatKeys[j]] == '.' and len(vals) == len(alt) + 1  and ('ref' in idToDescription['format'][formatKeys[j]] or 'Ref' in idToDescription['format'][formatKeys[j]]):
                            # this is a list with 1+len(alt) elements; the first is for REF and the rest are for the ALTs, in order.
                            splitRecord += vals[0] + ',' + vals[altNumber] + ':'
                        elif idToNumberType['format'][formatKeys[j]] == 'G':
                            splitRecord += ','.join(map(str,[vals[homRefIndex],vals[hetIndex],vals[homAltIndex]])) + ':'
                        else:
                            splitRecord += formatValues[sampleIndex][j] + ':'
                    # clean up: remove trailing ':' from sample-specific column
                    splitRecord = splitRecord.rstrip(':')
                print >> outFile, splitRecord
