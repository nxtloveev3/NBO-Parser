from .basicReadingFunctions import readlines, find, findExact, replacingDigit, fixEnding, finalClean

lines = readlines("CN01_NN13_S_nbo.log")
nboSumStart = find("NATURAL BOND ORBITALS (Summary):",lines)
startNLMO = find("NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:",lines)
nboSum = lines[nboSumStart[0]:startNLMO[0]]


def parseNBO(nboSum):
    c1 = [] #This part trying to find first element of each line in the file
    for line in nboSum:
        if len(line) > 0: 
            line = replacingDigit(0,line)
            c1.append(line[0])
        else:
            c1.append(line)
    posOrb = [] #Find the location of lines that have valuable infromation
    for i in range(len(c1)):
        if c1[i].isdigit():
            posOrb.append(i)
    tmpOrb = []
    for index in posOrb:
        tmpOrb.append(nboSum[index])
    posCR = find("CR",tmpOrb)
    posLP = find("LP",tmpOrb)
    posLV = find("LV",tmpOrb)
    posBD = findExact("BD",tmpOrb)
    posBDS = find("BD*(",tmpOrb)
    pos3C = findExact("3C",tmpOrb)
    pos3Cn = find("3Cn(",tmpOrb)
    pos3Cs = find("3C*(",tmpOrb)
    tabCR = fixEnding(posCR,tmpOrb)
    tabLP = fixEnding(posLP,tmpOrb)
    tabLV = fixEnding(posLV,tmpOrb)
    tabBD = fixEnding(posBD,tmpOrb)
    tabBDS = fixEnding(posBDS,tmpOrb)
    tab3C = fixEnding(pos3C,tmpOrb)
    tab3Cn = fixEnding(pos3Cs,tmpOrb)
    tab3Cs = fixEnding(pos3Cn,tmpOrb)
    tabCRF = finalClean(tabCR)
    tabLVF = finalClean(tabLV)
    tabBDF = finalClean(tabBD)
    tabBDSF = finalClean(tabBDS)

parseNBO(nboSum)