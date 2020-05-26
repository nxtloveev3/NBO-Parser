from basicReadingFunctions import find,readlines
lines = readlines("CN01_NN13_S_nbo.log")
startCMO = find("CMO: NBO Analysis of Canonical Molecular Orbitals",lines)
endCMO = find("Molecular Orbital Atom-Atom Bonding Character",lines)
cmoAlpha = lines[startCMO[0]:endCMO[0]]
print(cmoAlpha)

def parseCMO(file):
    return none




    