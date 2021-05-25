import sys
from miscFuncs import *

vcfFileName, maskedRefFileName, fastaFileName1, fastaFileName2 = sys.argv[1:]

# SV: filenames are of form /globalfs/ulb/ebe/slukiche/introgression/filet/15_vcf_to_fasta.vcf_X/gi_gq_ctgXXXX_phased.vcf.gz
# SV: on VEGA filenames are of form ./vcf_X/gi_gq_ctgXXXX_phased.vcf
arm = vcfFileName.split("_")[3]

# SV: definition of names on individuals in population 1 (GI) and population 2 (GQ)
pop1_names = ["A00387_98_GW190421167th_NovaXP_2_AAGACCGT",
              "A00387_98_GW190421167th_NovaXP_2_AGTGACCT",
              "A00387_98_GW190421167th_NovaXP_2_CAGGTTCA",
              "A00387_98_GW190421167th_NovaXP_2_CGAATTGC",
              "A00387_98_GW190421167th_NovaXP_2_TACGGTCT",
              "A00387_98_GW190421167th_NovaXP_2_TACTCCAG",
              "A00387_98_GW190421167th_NovaXP_2_TCGGATTC",
              "A00387_98_GW190421167th_NovaXP_2_TCGTCTGA",
              "A00387_98_GW190421167th_NovaXP_2_TTACCGAC",
              "A00387_98_GW190421167th_NovaXP_2_TTCCAGGT"]
pop2_names = ["A00387_98_GW190421167th_NovaXP_2_AGCCTATC",
              "A00387_98_GW190421167th_NovaXP_2_CCAGTATC",
              "A00387_98_GW190421167th_NovaXP_2_CTGTACCA",
              "A00387_98_GW190421167th_NovaXP_2_GAACGAAG",
              "A00387_98_GW190421167th_NovaXP_2_GGAAGAGA",
              "A00387_98_GW190421167th_NovaXP_2_TAGGAGCT",
              "A00387_98_GW190421167th_NovaXP_2_TCATCTCC",
              "A00387_98_GW190421167th_NovaXP_2_TCTACGCA",
              "A00387_98_GW190421167th_NovaXP_2_TTGCGAGA"]

# SV: indicates whether this is an individual from the first population
def isPop1ToKeep(name):
    return name in pop1_names

# SV: indicates whether this is an inndividual from the second population
def isPop2ToKeep(name):
    return name in pop2_names

def alleleNumToBase(alleleNum, ref, alt):
    assert alleleNum in ["0","1"]
    if alleleNum == "0":
        return ref
    else:
        return alt

def getInbredHap(genos, ref, alt):
    allele1, allele2 = genos.split("|")
    if allele1 == allele2:
        return alleleNumToBase(allele1, ref, alt)
    else:
        return 'N'

def getDiploidHaps(genos, ref, alt):
    allele1, allele2 = genos.split("|")
    return alleleNumToBase(allele1, ref, alt), alleleNumToBase(allele2, ref, alt)

def writeFastaFileWithRefAndSnpGenos(names, snpGenos, arm, maskedData, fastaFileName):
    with open(fastaFileName, "w") as outFile:
        for name in names:
            sys.stderr.write("writing %s for %s\n" %(name, arm))
            outS = ">%s\n" %(name)
            for i in range(len(maskedData[arm])):
                pos = i+1
                if snpGenos.has_key((arm, pos)) and maskedData[arm][i] != 'N':
                    outS += snpGenos[(arm, pos)][name]
                else:
                    outS += maskedData[arm][i]
                if pos % 60 == 0:
                    outS += '\n'
            outFile.write(outS + '\n')
            if pos % 60 != 0:
                outFile.write('\n')

def readSnpHapsFromPhasedSimSechVcf(vcfFileName):
    pop1Indices, pop2Indices = [], []
    pop1Names, pop2Names = [], []
    snpGenosPop1, snpGenosPop2 = {}, {}
    with open(vcfFileName) as vcfFile:
        for line in vcfFile:
            if line.startswith("#CHROM"):
                line = line.strip().split()
                for i in range(9, len(line)):
                    #TODO: the line below checks to see if we are examining a
                    # sechellia genome in our data set. See function definition
                    # at the top of the file and modify as necessary to pull
                    # out individuals in population 2 of your data set
                    if isPop1ToKeep(line[i]):
                        pop1Indices.append(i)
                        pop1Names.append(line[i]+".1") # SV: two times to handle the diploidy
                        pop1Names.append(line[i]+".2")

                    #TODO: the line below checks to see if we are examining a
                    # simulans genome in our data set. See function definition
                    # at the top of the file and modify as necessary to pull
                    # out individuals in population 1 of your data set.
                    elif isPop2ToKeep(line[i]):
                        pop2Indices.append(i)
                        pop2Names.append(line[i]+".1") # SV: two times to handle the diploidy
                        pop2Names.append(line[i]+".2")
                header = line
            elif not line.startswith("#"):
                line = line.strip().split()
                c, pos, varId, ref, alt = line[:5]
                assert ref in ['G','T','C','A'] and alt in ['G','T', 'C', 'A']
                pos = int(pos)
                snpGenosPop1[(c, pos)] = {}
                snpGenosPop2[(c, pos)] = {}
                for i in pop1Indices:
                    snpGenosPop1[(c, pos)][header[i]+".1"], snpGenosPop1[(c, pos)][header[i]+".2"] = getDiploidHaps(line[i], ref, alt)
                for i in pop2Indices:
                    snpGenosPop2[(c, pos)][header[i]+".1"], snpGenosPop2[(c, pos)][header[i]+".2"] = getDiploidHaps(line[i], ref, alt)
    return pop1Names, snpGenosPop1, pop2Names, snpGenosPop2

headersPop1, snpGenosPop1, headersPop2, snpGenosPop2 = readSnpHapsFromPhasedSimSechVcf(vcfFileName)

maskedData = readFa(maskedRefFileName, upper=True)
writeFastaFileWithRefAndSnpGenos(headersPop1, snpGenosPop1, arm, maskedData, fastaFileName1)
writeFastaFileWithRefAndSnpGenos(headersPop2, snpGenosPop2, arm, maskedData, fastaFileName2)
