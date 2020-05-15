def readlines(file):
    f = open(file,'rb')
    lines = f.readlines()
    lines = [line.decode("utf-8").strip() for line in lines]
    return lines

class nbo(object):
    def __init__(self,file):
        lines = readlines(file)
        nboSumStart = lines.index("NATURAL BOND ORBITALS (Summary):")
        startNAO = lines.index("NATURAL POPULATIONS:  Natural atomic orbital occupancies")
        endNAO = lines.index("Summary of Natural Population Analysis:")
        #startNBOalpha = lines.index("NATURAL BOND ORBITAL ANALYSIS, alpha spin orbitals") # could not find this line
        #startNBObeta = lines.index("NATURAL BOND ORBITAL ANALYSIS, beta spin orbitals")
        #endNBO = lines.index("NHO DIRECTIONALITY AND BOND BENDING") #could not find this line
        startCMO = lines.index("CMO: NBO Analysis of Canonical Molecular Orbitals")
        endCMO = lines.index("Molecular Orbital Atom-Atom Bonding Character")
        startPert = lines.index("SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS")
        endPert = lines.index("NATURAL BOND ORBITALS (Summary):")
        startNLMO = lines.index("NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:")
        endNLMO = lines.index("NATURAL POPULATIONS:  Natural atomic orbital occupancies")
        #endNLMO2 = lines.index("NBO analysis completed") #could not find the line

        self.nboAll = lines[startNAO:endNAO]
        self.nboSumAlpha = lines[nboSumStart:startNLMO] #What exactly is the index in mathmatica


nbo("CN01_NN01_opt_opt.log_r_nbo.log")