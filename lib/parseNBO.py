from basicReadingFunctions import readlines, find

lines = readlines("CN01_NN13_S_nbo.log")
nboSumStart = find("NATURAL BOND ORBITALS (Summary):",lines)
startNLMO = find("NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:",lines)
nboSum = lines[nboSumStart[0]:startNLMO[0]]


def parseNBO(nboSum):
    c1 = [] #This part trying to find first element of each line in the file
    for line in nboSum:
        if len(line) > 0: 
            line = line.split()
            if line[0].endswith("."): #Treat the number so it is easier to recognize as digit (Ex: 120. -> 120)
                line[0] = line[0][:-1]
            c1.append(line[0])
        else:
            c1.append(line)
    posOrb = [] #Find the location of lines that have valuable infromation
    for elem in c1:
        if elem.isdigit():
            posOrb.append(c1.index(elem))
    tmpOrb = []
    for index in posOrb:
        tmpOrb.append(nboSum[index])
    posCR = find("CR",tmpOrb)
    posLP = find("LP",tmpOrb)
    posLV = find("LV",tmpOrb)
    posBD = find("BD",tmpOrb)
    posBDS = find("BD*(",tmpOrb)
    

parseNBO(nboSum)