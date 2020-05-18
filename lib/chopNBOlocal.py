from basicReadingFunctions import find,findExact,readlines

def chopNBOlocal(file):
    lines = readlines(file)
    npaStart = find("Summary of Natural Population Analysis:",lines)
    pos1 = findExact("-------------------------------------------------------------------------------",lines)
    pos2 = findExact("===============================================================================",lines)
    pos11 = []
    pos22 = []
    for elem in pos1:
        if npaStart[0]<elem<npaStart[1]:
            pos11.append(elem)
    for elem in pos2:
        if npaStart[0]<elem<npaStart[1]:
            pos22.append(elem)
    npa = lines[pos11[0]+1:pos22[0]]
    ats = []
    posBadAt = []
    length = 8
    for line in npa:
        ats.append(line[0])
        count = 0
        for elem in line.split():
            count += 1
        if count != 8:
            posBadAt.append(npa.index(line))
    badAts = []
    if len(posBadAt) == 1:
        badAts.append(ats[posBadAt[0]])
    elif len(posBadAt) > 1:
        for num in range(len(posBadAt)):
            badAts.append(ats[posBadAt[num]])
    
#chopNBOlocal("CN01_NN02_T_nbo.log")
