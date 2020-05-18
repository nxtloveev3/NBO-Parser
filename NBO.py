from lib import basicReadingFunctions

# These two function extract the information needed and seperate them to different seperate files
# based on the origin file(weather it is a triplet or singlet)
def tripletFile(file):
    nboSumStart = basicReadingFunctions.find("NATURAL BOND ORBITALS (Summary):",file)
    startNAO = basicReadingFunctions.find("NATURAL POPULATIONS:  Natural atomic orbital occupancies",file)
    endNAO = basicReadingFunctions.find("Summary of Natural Population Analysis:",file)
    startNBOalpha = basicReadingFunctions.find("NATURAL BOND ORBITAL ANALYSIS, alpha spin orbitals:",file)
    startNBObeta = basicReadingFunctions.find("NATURAL BOND ORBITAL ANALYSIS, beta spin orbitals:",file)
    endNBO = basicReadingFunctions.find("NHO DIRECTIONALITY AND BOND BENDING",file)
    startCMO = basicReadingFunctions.find("CMO: NBO Analysis of Canonical Molecular Orbitals",file)
    endCMO = basicReadingFunctions.find("Molecular Orbital Atom-Atom Bonding Character",file)
    startPert = basicReadingFunctions.find("SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS",file)
    endPert = basicReadingFunctions.find("NATURAL BOND ORBITALS (Summary):",file)
    startNLMO = basicReadingFunctions.find("NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:",file)
    endNLMO = basicReadingFunctions.find("NATURAL POPULATIONS:  Natural atomic orbital occupancies",file)
    endNLMO2 = basicReadingFunctions.find("NBO analysis completed",file) 
    nboSumAlpha = file[nboSumStart[0]:startNLMO[0]]
    nboSumBeta = file[nboSumStart[1]:startNLMO[1]]
    nlmoAlpha = file[startNLMO[0]:endNLMO[0]]
    nlmoBeta = file[startNLMO[1]:endNLMO2[0]]
    naoAll = file[startNAO[0]:endNAO[0]]
    naoAlpha = file[startNAO[1]:endNAO[1]]
    naoBeta = file[startNAO[2]:endNAO[2]]
    nboAlpha = file[startNBOalpha[0]:endNBO[0]]
    nboBeta = file[startNBObeta[0]:endNBO[1]]
    cmoAlpha = file[startCMO[0]:endCMO[0]]
    cmoBeta = file[startCMO[1]:endCMO[1]]
    pertAlpha = file[startPert[0]:endPert[0]]
    pertBeta = file[startPert[1]:endPert[1]]
    return (nboSumAlpha,nboSumBeta,nlmoAlpha,nlmoBeta,naoAll,naoAlpha,naoBeta,nboAlpha,nboBeta,cmoAlpha,cmoBeta,pertAlpha,pertBeta)    
    
def singletFile(file):
    nboSumStart = basicReadingFunctions.find("NATURAL BOND ORBITALS (Summary):",file)
    startNAO = basicReadingFunctions.find("NATURAL POPULATIONS:  Natural atomic orbital occupancies",file)
    endNAO = basicReadingFunctions.find("Summary of Natural Population Analysis:",file)
    startNBOalpha = basicReadingFunctions.find("NATURAL BOND ORBITAL ANALYSIS",file)
    endNBO = basicReadingFunctions.find("NHO DIRECTIONALITY AND BOND BENDING",file)
    startCMO = basicReadingFunctions.find("CMO: NBO Analysis of Canonical Molecular Orbitals",file)
    endCMO = basicReadingFunctions.find("Molecular Orbital Atom-Atom Bonding Character",file)
    startPert = basicReadingFunctions.find("SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS",file)
    endPert = basicReadingFunctions.find("NATURAL BOND ORBITALS (Summary):",file)
    startNLMO = basicReadingFunctions.find("NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:",file)
    endNLMO2 = basicReadingFunctions.find("NBO analysis completed",file) 
    nboSum = file[nboSumStart[0]:startNLMO[0]]
    nlmo = file[startNLMO[0]:endNLMO[0]]
    nao = file[startNAO[0]:endNAO[0]]
    nbo = file[startNBOalpha[0]:endNBO[0]]
    cmo = file[startCMO[0]:endCMO[0]]
    pert = file[startPert[0]:endPert[0]]
    return (nboSum,nlmo,nao,nbo,cmo,pert)

class nbo(object):
    def __init__(self,file):
        lines = basicReadingFunctions.readlines(file)
        if "alpha spin orbitals" or "beta spin orbitals" in lines:
            self.nboSumAlpha,self.nboSumBeta,self.nlmoAlpha,self.nlmoBeta,self.naoAll,self.naoAlpha,self.naoBeta,self.nboAlpha,self.nboBeta,self.cmoAlpha,self.cmoBeta,self.pertAlpha,self.pertBeta = tripletFile(lines)
        else:
            self.nboSum,self.nlmo,self.nao,self.nbo,self.cmo,self.pert = singletFile(lines)
        print(self.nboSumAlpha)
        

nbo("CN01_NN02_T_nbo.log")