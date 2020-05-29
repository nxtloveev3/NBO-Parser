from .basicReadingFunctions import find,readlines,findExact,finalClean
import string

def parseCMO(cmoA):
    sp = []
    for line in cmoA: 
        newline = []
        if line not in string.whitespace and isinstance(line,str):
            for elem in line.split():
                newline.append(elem)
            sp.append(newline)
    posB = findExact("MO",cmoA)
    spp = []
    for num in posB:
        spp.append(sp[num])
    cmoTab = finalClean(spp)
    return cmoTab