import sys
from miscFuncs import *

vcfSnpFileName1, vcfSnpFileName2, vcfAllSitesFileName1, vcfAllSitesFileName2, maskedFaFileName, qualThreshold, outVcfPrefix = sys.argv[1:]
qualThreshold = int(qualThreshold)

# SV: this function is used to initialize all the SNP locations
# SV: this is the only use of the first SNP VCF files
def addSnpLocsInVcfFile(vcfFileName, snpLocs, outGenos, quals):
    with open(vcfFileName) as vcfFile:
        for line in vcfFile:
            if not line.startswith("#"):
                arm, pos, varId, ref, alt, qual, varFilter, info, fmt = line.strip().split("\t")[:9]
                pos = int(pos)
                snpLocs[(arm, pos)] = 0
                outGenos[(arm, pos)] = []
                quals[(arm, pos)] = []

maskedGenome = readFaAsLists(maskedFaFileName, upper=True)

arms = maskedGenome.keys()
outFiles = {}
for arm in arms:
#    outFileName = outVcfPrefix + arm + ".vcf"
#    outFiles[arm] = open(outFileName, "w")
    outFiles[arm] = ""

snpLocs = {}
outGenos = {}
quals = {}
addSnpLocsInVcfFile(vcfSnpFileName1, snpLocs, outGenos, quals)
addSnpLocsInVcfFile(vcfSnpFileName2, snpLocs, outGenos, quals)
print len(snpLocs)

firstFile = True
visitedSnpCount = [0,0]
fileIndex = 0
writtenLineCount = 0
altAlleles = {}
for vcfFileName in [vcfAllSitesFileName1, vcfAllSitesFileName2]:
    with open(vcfFileName) as vcfFile:
        for line in vcfFile:
            if line.startswith("#CHROM\tPOS\t"):
                if firstFile:
                    header1 = line.strip()  # SV: stores the header of the first VCF file
                else:
                    fileIndex += 1
                    header2Samps = line.strip().split("\tFORMAT\t")
                    assert len(header2Samps) == 2
                    header2Samps = header2Samps[1]  # SV: stores the format of the second VCF file
                    for arm in arms:
                        outFiles[arm] += header1 + "\t" + header2Samps + "\n"  # SV: creates the header based on all the samples for the resulting VCF files
            elif not line.startswith("#"):
                line = line.strip().split("\t")
                arm, pos, varId, ref, alt, qual, varFilter, info, fmt = line[:9]
                pos = int(pos)
                if snpLocs.has_key((arm, pos)) and maskedGenome[arm][pos-1] != 'N':  # SV: there is a non-masked SNP at this position
                    snpLocs[(arm, pos)] += 1  # SV: stores the number of occurrences of this SNP (0, 1 or 2 if in both species)
                    if alt != ".":
                        if not altAlleles.has_key((arm, pos)):
                            altAlleles[(arm, pos)] = {}
                        altAlleles[(arm, pos)][alt] = 1  # SV: counts the number of alternate alleles for this position
                    visitedSnpCount[fileIndex] += 1
                    if qual == ".":
                        qual = 0.0
                    else:
                        qual = float(qual)
                    quals[(arm, pos)].append(qual)

                    # SV: creation of new format
                    genoFormat = fmt.split(":")
                    genoIndex = None
                    gqIndex = None
                    rgqIndex = None
                    newGenoFormat = "GT:GQ"
                    for j in range(len(genoFormat)):
                        if genoFormat[j] == "GQ":
                            gqIndex = j
                        elif genoFormat[j] == "RGQ":
                            rgqIndex = j
                        elif genoFormat[j] == "GT":
                            genoIndex = j

                    # SV: creation of entry for each individual
                    for i in range(9, len(line)):
                        gtLs = line[i].split(":")
                        assert gtLs[genoIndex] != "1/0"
                        if gtLs[genoIndex] == "./." or (len(gtLs) == 1 and gtLs[genoIndex] == "0/0"):
                            gq = 0
                        elif rgqIndex != None and gtLs[rgqIndex] != ".":
                            gq = int(gtLs[rgqIndex])
                        elif gtLs[gqIndex] != ".":
                            gq = int(gtLs[gqIndex])
                        else:
                            gq = 0
                        if gq < qualThreshold:
                            outGenos[(arm, pos)].append("./." + ":" + str(gq))
                        else:
                            outGenos[(arm, pos)].append(gtLs[genoIndex] + ":" + str(gq))
                    if not firstFile:
                        if not altAlleles.has_key((arm, pos)):
                            print "allSites vs. varSites discrepancy:", arm, pos, maskedGenome[arm][pos-1]
                        else:
                            if len(altAlleles[(arm, pos)]) != 1:
                                print arm, pos, altAlleles[(arm, pos)]
                            outLine = "\t".join(line[:4] + [altAlleles[(arm, pos)].keys()[0], str(max(quals[(arm, pos)])), ".", ".", newGenoFormat] + outGenos[(arm, pos)])
                            outFiles[arm] += outLine + "\n"
                            writtenLineCount += 1
    firstFile = False

for arm in outFiles:
    outFileName = outVcfPrefix + arm + ".vcf"
    outFile = open(outFileName, "w")
    outFile.write(outFiles[arm])
    outFile.close()

print visitedSnpCount
print writtenLineCount

for arm, pos in snpLocs:
    if not snpLocs[(arm, pos)] in [0,2]:
        print arm, pos, snpLocs[(arm, pos)]
