from .lib import basicReadingFunctions as brf
import copy
import re

# These two function extract the information needed and seperate them into different seperate files
# based on the origin file(whether it is a triplet or singlet)
def unrestricted(file, options=None):

    tabs = {}
    
    if 'nboSum' in options:
        start = brf.find("SUMMARY OF NATURAL BOND ORBITALS",file)
        assert(len(start)==2)
        start_alpha, start_beta = start

        end_total = brf.find("TOTAL LEWIS",file)
        end_valence = brf.find("VALENCE NON-NEWIS",file)
        end_alpha, end_beta = 0, 0
        for i in end_total:
            for j in end_valence:
                if j == i+1:
                    if i < start_beta and j < start_beta:
                        end_alpha = max(end_alpha, i)
                    else:
                        assert(i > start_beta and j > start_beta)
                        end_beta = max(end_beta, i)
        tabs['nboSumAlpha'], tabs['nboSumBeta'] = file[start_alpha:end_alpha+5], file[start_beta:end_beta+5]

    if 'nao' in options:
        startNAO = brf.find("NATURAL ATOMIC ORBITAL OCCUPANCIES",file)
        endNAO = brf.find("SUMMARY OF NATURAL POPULATION ANALYSIS",file)
        naoAll = file[startNAO[0]:endNAO[0]]
        naoAlpha = file[startNAO[1]:endNAO[1]]
        naoBeta = file[startNAO[2]:endNAO[2]]
        tabs['naoAll'], tabs['naoAlpha'], tabs['naoBeta'] = naoAll, naoAlpha, naoBeta
    
    return tabs
    
def restricted(file, options=None):

    tabs = {}

    if 'nboSum' in options:
        nboSumStart = brf.find("SUMMARY OF NATURAL BOND ORBITALS (TOTAL DENSITY):",file)
        nboSumEnd = brf.find("NATURAL LOCALIZED MOLECULAR ORBITAL ANALYSIS (TOTAL DENSITY):",file)
        tabs['nboSum'] = file[nboSumStart[0]:nboSumEnd[0]]

    if 'nao' in options:
        startNAO = brf.find("NATURAL ATOMIC ORBITAL OCCUPANCIES (TOTAL DENSITY):",file)
        endNAO = brf.find("SUMMARY OF NATURAL POPULATION ANALYSIS (TOTAL DENSITY):",file)
        tabs['nao'] = file[startNAO[0]:endNAO[0]]
    
    return tabs

class nbo(object):
    def __init__(self, file, options=None, triplet=False, determine_triplet=False):
        lines = brf.readlines(file)
        self.triplet = triplet
        if determine_triplet:
            for line in lines: # Determines if the file is a singlet file or triplet file
                if "alpha spin orbitals" in line:
                    self.triplet = True
                    break
        self.npa,self.badAts,self.badAtsF = findNpa(file)
        if 'npa' in options:
            self.npa = nbo.replacement(self.npa,self.badAts,self.badAtsF)
        else:
            del self.npa
        if self.triplet:
            tabs = unrestricted(lines, options=options)
            if len(re.compile('\W').split(options)) > 1:
                for key in options:
                    if key != 'npa':
                        setattr(self, key+'A', nbo.replacement(tabs[key+'Alpha'],self.badAts,self.badAtsF))
                        setattr(self, key+'B', nbo.replacement(tabs[key+'Beta'],self.badAts,self.badAtsF))
                        if key == 'nao':
                            setattr(self, 'nao', nbo.replacement(tabs['naoAll'],self.badAts,self.badAtsF))
            else:
                setattr(self, options, nbo.replacement(tabs[options],self.badAts,self.badAtsF))
                setattr(self, options, nbo.replacement(tabs[options],self.badAts,self.badAtsF))
                if options == 'nao':
                    setattr(self, 'nao', nbo.replacement(tabs['naoAll'],self.badAts,self.badAtsF))            
        else:
            tabs = restricted(lines, options=options)
            if len(re.compile('\W').split(options)) > 1:
                for key in options:
                    if key != 'npa':
                        setattr(self, key, nbo.replacement(tabs[key],self.badAts,self.badAtsF))
            else:
                setattr(self, options, nbo.replacement(tabs[options],self.badAts,self.badAtsF))
    
    #This method replaces all the incorrect characters with correct ones
    @staticmethod
    def replacement(file,incorrect,correct):
        result = file.copy()
        for line in result:
            for num in range(len(incorrect)):
                if incorrect[num] in line:
                    str(result).replace(incorrect[num],correct[num])
        return list(result)


def findNpa(text):
    lines = brf.readlines(text)
    npaStart = brf.find("SUMMARY OF NATURAL POPULATION ANALYSIS (TOTAL DENSITY):",lines)
    pos1 = brf.find("-----------------------------------------------------------------------------",lines)
    pos2 = brf.find("=============================================================================",lines)
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