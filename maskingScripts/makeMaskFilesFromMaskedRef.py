import sys,os
import miscFuncs

maskedRefFaFileName, maskFileDir, winSize, max_masked_threshold = sys.argv[1:]
winSize = int(winSize)
max_masked_threshold = float(max_masked_threshold)
os.system("mkdir -p %s" %(maskFileDir))

maskedGenome = miscFuncs.readFa(maskedRefFaFileName)

# SV: provides coordinates of masked regions on the provided sequence
def getMaskedRuns(seq, winSize):
    inRun = False
    runs = []
    masked_len = 0
    masked_percentage = 0
    for i in range(len(seq)):
        if inRun:
            if seq[i] == 'N':
                end = i
                masked_len += 1
            else:
                runs.append((start/float(winSize), end/float(winSize)))
                inRun = False
        else:
            if seq[i] == 'N':
                start = i
                end = i
                inRun = True
                masked_len += 1
            else:
                pass
    if inRun:
        runs.append((start/float(winSize), end/float(winSize)))
    masked_percentage = float(masked_len) / len(seq)
    if masked_percentage > max_masked_threshold:
        runs = None
    return runs, masked_percentage

def writeMaskedRegionsToFile(seq, maskFileName, winSize):
    assert len(seq) == winSize
    maskedRuns, masked_percentage = getMaskedRuns(seq, winSize)
    if maskedRuns != None:
        maskFileName = maskFileName + "_" + str(float(1) - masked_percentage)
        maskFile = open(maskFileName, "w")
        for s, e in maskedRuns:
            maskFile.write("0 %f %f\n" %(s, e))
        maskFile.write("//\n")
        maskFile.close()

for arm in maskedGenome:
    prevEnd = 0
    while prevEnd < len(maskedGenome[arm]):
        currStart = prevEnd+1
        currEnd = currStart + winSize-1
        if currEnd <= len(maskedGenome[arm]):
            maskFileName = maskFileDir + "/%s.%d.%d.maskedRegions" %(arm, currStart, currEnd)
            writeMaskedRegionsToFile(maskedGenome[arm][currStart-1:currEnd], maskFileName, winSize)
        prevStart = currStart
        prevEnd = currEnd
