from .lib import basicReadingFunctions as brf
import copy
import re

# These two function extract the information needed and seperate them into different seperate files
# based on the origin file(whether it is a triplet or singlet)
def unrestricted(file, selected=None):

    tabs = {}

    if 'nboSum' in selected:
        nboSumStart = brf.find("NATURAL BOND ORBITALS (Summary):",file)
        nboSumEnd = brf.find("Total Lewis",file)
        nboSumAlpha = file[nboSumStart[0]:nboSumEnd[len(nboSumEnd)-4]]
        nboSumBeta = file[nboSumStart[1]:nboSumEnd[len(nboSumEnd)-1]]
        tabs['nboSumAlpha'], tabs['nboSumBeta'] = nboSumAlpha, nboSumBeta

    if 'nao' in selected:
        startNAO = brf.find("NATURAL POPULATIONS:  Natural atomic orbital occupancies",file)
        endNAO = brf.find("Summary of Natural Population Analysis:",file)
        naoAll = file[startNAO[0]:endNAO[0]]
        naoAlpha = file[startNAO[1]:endNAO[1]]
        naoBeta = file[startNAO[2]:endNAO[2]]
        tabs['naoAll'], tabs['naoAlpha'], tabs['naoBeta'] = naoAll, naoAlpha, naoBeta
    
    if 'nbo' in selected:
        startNBOalpha = brf.find("NATURAL BOND ORBITAL ANALYSIS, alpha spin orbitals:",file)
        startNBObeta = brf.find("NATURAL BOND ORBITAL ANALYSIS, beta spin orbitals:",file)
        endNBO = brf.find("NHO DIRECTIONALITY AND BOND BENDING",file)
        nboAlpha = file[startNBOalpha[0]:endNBO[0]]
        nboBeta = file[startNBObeta[0]:endNBO[1]]
        tabs['nboAlpha'], tabs['nboBeta'] = nboAlpha, nboBeta

    if 'nlmo' in selected:
        startNLMO = brf.find("NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:",file)
        endNLMO = brf.find("*******         Beta  spin orbitals         *******",file)
        endNLMO2 = brf.find("NBO analysis completed",file) 
        nlmoAlpha = file[startNLMO[0]:endNLMO[0]]
        nlmoBeta = file[startNLMO[1]:endNLMO2[0]]
        tabs['nlmoAlpha'], tabs['nlmoBeta'] = nlmoAlpha, nlmoBeta
    
    if 'cmo' in selected:
        startCMO = brf.find("CMO: NBO Analysis of Canonical Molecular Orbitals",file)
        startPert = brf.find("SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS",file)
        cmoAlpha = file[startCMO[0]:startPert[0]]
        cmoBeta = file[startCMO[1]:startPert[1]]
        tabs['cmoAlpha'], tabs['cmoBeta'] = cmoAlpha, cmoBeta
    
    if 'pert' in selected:
        if 'CMO' not in selected:
            startPert = brf.find("SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS",file)
        endPert = brf.find("NATURAL BOND ORBITALS (Summary):",file)
        pertAlpha = file[startPert[0]:endPert[0]]
        pertBeta = file[startPert[1]:endPert[1]]
        tabs['pertAlpha'], tabs['pertBeta'] = pertAlpha, pertBeta
    
    return tabs
    
def restricted(file, selected=None):

    tabs = {}

    if 'nboSum' in selected:
        nboSumStart = brf.find("NATURAL BOND ORBITALS (Summary):",file)
        nboSumEnd = brf.find("$END",file)
        tabs['nboSum'] = file[nboSumStart[0]:nboSumEnd[0]]
    
    if 'nao' in selected:
        startNAO = brf.find("NATURAL POPULATIONS:  Natural atomic orbital occupancies",file)
        endNAO = brf.find("Summary of Natural Population Analysis:",file)
        tabs['nao'] = file[startNAO[0]:endNAO[0]]
    
    if 'nbo' in selected:
        startNBO = brf.find("NATURAL BOND ORBITAL ANALYSIS",file)
        endNBO = brf.find("NHO DIRECTIONALITY AND BOND BENDING",file)
        tabs['nbo'] = file[startNBO[0]:endNBO[0]]
    
    if 'nlmo' in selected:
        startNLMO = brf.find("NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:",file)
        endNLMO = brf.find("NBO analysis completed",file) 
        tabs['nlmo'] = file[startNLMO[0]:endNLMO[0]]
    
    if 'cmo' in selected:
        startCMO = brf.find("CMO: NBO Analysis of Canonical Molecular Orbitals",file)
        endCMO = brf.find("Molecular Orbital Atom-Atom Bonding Character",file)
        tabs['cmo'] = file[startCMO[0]:endCMO[0]]
    
    if 'pert' in selected:
        startPert = brf.find("SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS",file)
        endPert = brf.find("NATURAL BOND ORBITALS (Summary):",file)
        tabs['pert'] = file[startPert[0]:endPert[0]]
    
    return tabs

class nbo(object):
    def __init__(self, file, selected=None, triplet=False, determine_triplet=False):
        lines = brf.readlines(file)
        self.triplet = triplet
        if determine_triplet:
            for line in lines: # Determines if the file is a singlet file or triplet file
                if "alpha spin orbitals" in line:
                    self.triplet = True
                    break
        self.npa,self.badAts,self.badAtsF = findNpa(file)
        if 'npa' in selected:
            self.npa = nbo.replacement(self.npa,self.badAts,self.badAtsF)
        else:
            del self.npa
        if self.triplet:
            tabs = unrestricted(lines, selected=selected)
            for key in selected:
                if key != 'npa':
                    setattr(self, key+'A', nbo.replacement(tabs[key+'Alpha'],self.badAts,self.badAtsF))
                    setattr(self, key+'B', nbo.replacement(tabs[key+'Beta'],self.badAts,self.badAtsF))
                    if key == 'nao':
                        setattr(self, 'nao', nbo.replacement(tabs['naoAll'],self.badAts,self.badAtsF))
        else:
            tabs = restricted(lines, selected=selected)
            for key in selected:
                if key != 'npa':
                    setattr(self, key, nbo.replacement(tabs[key],self.badAts,self.badAtsF))
    
    #This method replaces all the incorrect characters with correct ones
    @staticmethod
    def replacement(file,incorrect,correct):
        result = file.copy()
        for line in result:
            for num in range(len(incorrect)):
                if incorrect[num] in line:
                    result.replace(incorrect[num],correct[num])
        return result

#This is a method to locate the summary of natural population analysis and determine if there is 
#incorrect character due to the Gaussian generation and find the correct way. (Ex: C125 -> C 125)
 
def findNpa(text):
    lines = brf.readlines(text)
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
    ats, badAtsF = [line[0] for line in npa], []
    badAts = [atom for atom in ats if not atom.isalpha()]
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
def parseNPA(npa, triplet=False) -> list:
    columns = ['Atom', 'No', 'Natural Charge', 'Core', 'Valence', 'Rydeberg', 'Total']
    length = 7 + int(triplet)
    if triplet:
        columns.append('Natural Spin Density')
    def helper(line):
        element = ''.join([i for i in line[0] if i.isalpha()])
        return [element, line[0].replace(element, '')] + line[1:]
    text = [i.split() for i in npa]
    result = []
    for line in text:
        if len(line) != length:
            line = helper(line)
        new = {'Atom': line[0]}
        new['No'] = int(line[1])
        for i in range(2, length):
            new[columns[i]] = float(line[i])
        result.append(new)
    return result

'''
Finds the final geometry. 
Set nbo_log=True when parsing a log file of nbo calculation. 
'''
_coord = lambda name: "(?P<" + name + ">" + r'-?\d+\.+\d*' + ")"
coord_line = "(?P<atom>" + r'[A-Z][a-z]?' +")"
coord_line += r','
coord_line += _coord('x')
coord_line += r','
coord_line += _coord('y')
coord_line += r','
coord_line += _coord('z')
coordRe = re.compile(coord_line)
coord_line_nbo = "(?P<atom>" + r'[A-Z][a-z]?' +")"
coord_line_nbo += r',0'
coord_line_nbo += r','
coord_line_nbo += _coord('x')
coord_line_nbo += r','
coord_line_nbo += _coord('y')
coord_line_nbo += r','
coord_line_nbo += _coord('z')
coordRe_nbo = re.compile(coord_line_nbo)

def parseXYZ(filename, nbo_log=False):
    with open(filename, 'r') as file:
        text = file.read().replace('\n', '').replace(' ', '')
    if nbo_log:
        r = coordRe_nbo
    else:
        r = coordRe
    def helper(x):
        for i in ['x', 'y', 'z']:
            x[i] = float(x[i])
        return x
    return [helper(j.groupdict()) for j in r.finditer(text)]