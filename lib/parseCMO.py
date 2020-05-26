from .basicReadingFunctions import find,readlines,findExact,finalClean
import string

lines = readlines("CN01_NN13_S_nbo.log")
startCMO = find("CMO: NBO Analysis of Canonical Molecular Orbitals",lines)
endCMO = find("Molecular Orbital Atom-Atom Bonding Character",lines)
cmoAlpha = lines[startCMO[0]:endCMO[0]]

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


parseCMO(cmoAlpha)