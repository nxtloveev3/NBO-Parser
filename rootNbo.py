from .lib import basicReadingFunctions as brf
import copy

# These two function extract the information needed and seperate them into different seperate files
# based on the origin file(whether it is a triplet or singlet)
def unrestricted(file):
    nboSumStart = brf.find("NATURAL BOND ORBITALS (Summary):",file)
    startNAO = brf.find("NATURAL POPULATIONS:  Natural atomic orbital occupancies",file)
    endNAO = brf.find("Summary of Natural Population Analysis:",file)
    startNBOalpha = brf.find("NATURAL BOND ORBITAL ANALYSIS, alpha spin orbitals:",file)
    startNBObeta = brf.find("NATURAL BOND ORBITAL ANALYSIS, beta spin orbitals:",file)
    nboSumEnd = brf.find("Total Lewis",file)
    endNBO = brf.find("NHO DIRECTIONALITY AND BOND BENDING",file)
    startCMO = brf.find("CMO: NBO Analysis of Canonical Molecular Orbitals",file)
    endCMO = brf.find("Molecular Orbital Atom-Atom Bonding Character",file)
    startPert = brf.find("SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS",file)
    endPert = brf.find("NATURAL BOND ORBITALS (Summary):",file)
    startNLMO = brf.find("NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:",file)
    endNLMO = brf.find("*******         Beta  spin orbitals         *******",file)
    endNLMO2 = brf.find("NBO analysis completed",file) 
    nboSumAlpha = file[nboSumStart[0]:nboSumEnd[len(nboSumEnd)-3]]
    nboSumBeta = file[nboSumStart[1]:nboSumEnd[len(nboSumEnd)-1]]
    nlmoAlpha = file[startNLMO[0]:endNLMO[0]]
    nlmoBeta = file[startNLMO[1]:endNLMO2[0]]
    naoAll = file[startNAO[0]:endNAO[0]]
    naoAlpha = file[startNAO[1]:endNAO[1]]
    naoBeta = file[startNAO[2]:endNAO[2]]
    nboAlpha = file[startNBOalpha[0]:endNBO[0]]
    nboBeta = file[startNBObeta[0]:endNBO[1]]
    cmoAlpha = file[startCMO[0]:startPert[0]]
    cmoBeta = file[startCMO[1]:startPert[1]]
    pertAlpha = file[startPert[0]:endPert[0]]
    pertBeta = file[startPert[1]:endPert[1]]
    return (nboSumAlpha,nboSumBeta,nlmoAlpha,nlmoBeta,naoAll,naoAlpha,naoBeta,nboAlpha,nboBeta,cmoAlpha,cmoBeta,pertAlpha,pertBeta)    
    
def restricted(file):
    nboSumStart = brf.find("NATURAL BOND ORBITALS (Summary):",file)
    nboSumEnd = brf.find("$END",file)
    startNAO = brf.find("NATURAL POPULATIONS:  Natural atomic orbital occupancies",file)
    endNAO = brf.find("Summary of Natural Population Analysis:",file)
    startNBOalpha = brf.find("NATURAL BOND ORBITAL ANALYSIS",file)
    endNBO = brf.find("NHO DIRECTIONALITY AND BOND BENDING",file)
    startCMO = brf.find("CMO: NBO Analysis of Canonical Molecular Orbitals",file)
    endCMO = brf.find("Molecular Orbital Atom-Atom Bonding Character",file)
    startPert = brf.find("SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS",file)
    endPert = brf.find("NATURAL BOND ORBITALS (Summary):",file)
    startNLMO = brf.find("NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:",file)
    endNLMO = brf.find("NBO analysis completed",file) 
    nboSum = file[nboSumStart[0]:nboSumEnd[0]]
    nlmo = file[startNLMO[0]:endNLMO[0]]
    nao = file[startNAO[0]:endNAO[0]]
    nbo = file[startNBOalpha[0]:endNBO[0]]
    cmo = file[startCMO[0]:endCMO[0]]
    pert = file[startPert[0]:endPert[0]]
    return (nboSum,nlmo,nao,nbo,cmo,pert)

class nbo(object):
    def __init__(self,file):
        lines = brf.readlines(file)
        triplet = False
        for line in lines: # Determines if the file is a singlet file or triplet file
            if "alpha spin orbitals" in line:
                triplet = True
                break
        if triplet == True:
            self.nboSumAlpha,self.nboSumBeta,self.nlmoAlpha,self.nlmoBeta,self.naoAll,self.naoAlpha,self.naoBeta,self.nboAlpha,self.nboBeta,self.cmoAlpha,self.cmoBeta,self.pertAlpha,self.pertBeta = unrestricted(lines)
            self.npa,self.badAts,self.badAtsF = nbo.findNpa(file)
            self.npa = nbo.replacement(self.npa,self.badAts,self.badAtsF)
            self.nboSumA = nbo.replacement(self.nboSumAlpha,self.badAts,self.badAtsF)
            self.nboSumB = nbo.replacement(self.nboSumBeta,self.badAts,self.badAtsF)
            self.nlmoA = nbo.replacement(self.nlmoAlpha,self.badAts,self.badAtsF)
            self.nlmoB = nbo.replacement(self.nlmoBeta,self.badAts,self.badAtsF)
            self.nao = nbo.replacement(self.naoAll,self.badAts,self.badAtsF)
            self.naoA = nbo.replacement(self.naoAlpha,self.badAts,self.badAtsF)
            self.naoB = nbo.replacement(self.naoBeta,self.badAts,self.badAtsF)
            self.nboA = nbo.replacement(self.nboAlpha,self.badAts,self.badAtsF)
            self.nboB = nbo.replacement(self.nboBeta,self.badAts,self.badAtsF)
            self.cmoA = nbo.replacement(self.cmoAlpha,self.badAts,self.badAtsF)
            self.cmoB = nbo.replacement(self.cmoBeta,self.badAts,self.badAtsF)
            self.pertA = nbo.replacement(self.pertAlpha,self.badAts,self.badAtsF)
            self.pertB = nbo.replacement(self.pertBeta,self.badAts,self.badAtsF)
        else:
            self.nboSum,self.nlmo,self.nao,self.nbo,self.cmo,self.pert = restricted(lines)
            self.npa,self.badAts,self.badAtsF = nbo.findNpa(file)
            self.nboSum = nbo.replacement(self.nboSum,self.badAts,self.badAtsF)
            self.nlmo = nbo.replacement(self.nlmo,self.badAts,self.badAtsF)
            self.nao = nbo.replacement(self.nao,self.badAts,self.badAtsF)
            self.nbo = nbo.replacement(self.nbo,self.badAts,self.badAtsF)
            self.cmo = nbo.replacement(self.cmo,self.badAts,self.badAtsF)
            self.pert = nbo.replacement(self.pert,self.badAts,self.badAtsF)

    #This is a method to locate the summary of natural population analysis and determine if there is 
    #incorrect character due to the Gaussian generation and find the correct way. (Ex: C125 -> C 125)
    @staticmethod 
    def findNpa(file):
        lines = brf.readlines(file)
        npaStart = brf.find("Summary of Natural Population Analysis:",lines)
        pos1 = brf.find("--------------------------------------------------------------------",lines)
        pos2 = brf.find("====================================================================",lines)
        pos11 = []
        pos22 = []
        for elem in pos1:
            if elem > npaStart[0]:
                pos11.append(elem)
                break
        for elem in pos2:
            if elem > pos11[0]:
                pos22.append(elem)
                break
        npa = lines[pos11[0]+1:pos22[0]]
        ats = []
        badAts = []
        badAtsF = []
        for line in npa:
            ats.append(line[0])
        for atom in ats:
            if not atom.isalpha:
                position = ats.index(atom)
                badAt.append(ats[position])
        for elem in badAts:
            character = ""
            number = ""
            for char in elem:
                if char.isalpha():
                    character += char
                else:
                    number += char
            result = character + " " + number
            badAtsF.append(result)
        return (npa,badAts,badAtsF)
    
    #This method reparses NPA into a list of dictionaries. Could be converted into a dataframe directly.
    def parseNPA(self) -> list:
        def helper(line):
            element = ''.join([i for i in line[0] if i.isalpha()])
            return [element, line[0].replace(element, '')] + line[1:]
        columns = ['Atom', 'No', 'Natural Charge', 'Core', 'Valence', 'Rydeberg', 'Total']
        text = [i.split() for i in self.npa]
        result = []
        for line in text:
            if len(line) != 7:
                line = helper(line)
            new = {'Atom': line[0]}
            new['No'] = int(line[1])
            for i in range(2, 7):
                new[columns[i]] = float(line[i])
            result.append(new)
        return result
    
    #This method replaces all the incorrect characters with correct ones
    @staticmethod
    def replacement(file,incorrect,correct):
        result = file.copy()
        for line in result:
            for num in range(len(incorrect)):
                if incorrect[num] in line:
                    result.replace(incorrect[num],correct[num])
        return result