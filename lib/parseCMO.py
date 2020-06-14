import re
from .basicReadingFunctions import namedRe, find, findExact, extractTab
import pandas as pd

def parseCMON(file):
    loc = find("Molecular Orbital Atom-Atom Bonding Character",file)
    cmon = file[:loc[0]]
    tmpFile = []
    for line in cmon:
        newline = line
        if "BD*" in line:
            newline = line.replace("BD*","ABB") #ABB stand for antibonding bonds
        tmpFile.append(newline)

    ####NonBond Regex####
    reFloat = r"-?\d+\.\d+"
    nonBondType = r"\d*[A-Z][A-Z]?[a-z]?"
    atom = r"[A-Z][a-z]?"
    cmonLine = namedRe("Coef", reFloat, before="allow", after='none')
    cmonLine += r"\*\[" + namedRe("NBO", r"\d+", before="allow", after='none') + r"\]\:"
    cmonLine += namedRe("Type", nonBondType, before="require")
    cmonLine += r"\(" + namedRe("Index", r"\d+", before="allow", after='none') + r"\)"
    cmonLine += namedRe("Atom", atom, before="allow", after='allow')
    cmonLine += namedRe("Loc", r"\d+", before="allow", after='allow')
    cmonLineRe = re.compile(cmonLine)

    ####Bond Regex####
    bondType = r"\d*[A-Z][A-Z][A-Z]?"
    cmonBondLine = namedRe("Coef", reFloat, before="allow", after='none')
    cmonBondLine += r"\*\[" + namedRe("NBO", r"\d+", before='allow', after='none') + r"\]\:"
    cmonBondLine += namedRe("Type", bondType, before="require", after='allow')
    cmonBondLine += r"\(" + namedRe("Index", r"\d+",before="allow", after='none') + r"\)"
    cmonBondLine += namedRe("Atom1", atom, before="allow", after='allow')
    cmonBondLine += namedRe("Loc1", r"\d+", before="allow", after='none')
    cmonBondLine += r"\-" + namedRe("Atom2", atom, before='allow', after='none')
    cmonBondLine += namedRe("Loc2", r"\d+", before='allow', after='none') + r"\*?"
    cmonBondLineRe = re.compile(cmonBondLine)

    ####MO Regex####
    cmonMOLine = namedRe("Label1", r"[M][O]", before='allow')
    cmonMOLine += namedRe("MO", r"\d+")
    cmonMOLine += r"\(" + namedRe("Label2", r"\w+", after='none') + r"\)\:"
    cmonMOLine += namedRe("Label3", r"\w+\s+\w+", before="require")
    cmonMOLine += r"\=" + namedRe("Energy", reFloat, before="require")
    cmonMOLine += namedRe("Unit", r"[a]\.[u]\.",after='allow')
    cmonMOLineRe = re.compile(cmonMOLine)

    result = {}
    currentMO = None
    for line in tmpFile:
        if cmonMOLineRe.search(line):
            correct = cmonMOLineRe.search(line).groupdict()
            currentMO = correct["MO"]
            energy = correct["Energy"]
            result[currentMO] = dict()
            result[currentMO]["Energy"] = energy
            result[currentMO]["wf"] = []
        elif cmonLineRe.search(line):
            correct = cmonLineRe.search(line).groupdict()
            result[currentMO]["wf"].append(correct)
        elif cmonBondLineRe.search(line):
            correct = cmonBondLineRe.search(line).groupdict()
            result[currentMO]["wf"].append(correct)
        else:
            print("parseCMON: WARNING line not recognized "+ line)

    return result


def parseCMO2(file):
    loc = find("Molecular Orbital Atom-Atom Bonding Character",file)
    cmo2 = file[loc[0]:]

    ####Bonding Regex####
    reFloat = r"-?\d+\.\d+"
    atom = r"[A-Z][a-z]?\s*\d+"
    cmo2BondLine = r"\b" + namedRe("Coeff", reFloat, before='allow')
    cmo2BondLine += namedRe('Atom1', atom, after='none') + r"\-"
    cmo2BondLine += namedRe('Atom2', atom, before='allow', after='allow') + r"\b"
    cmo2BondLineRe = re.compile(cmo2BondLine)

    ####NonBonding Regex####
    cmo2NonBondLine = namedRe("Coeff", reFloat, before='allow')
    cmo2NonBondLine += namedRe('Atom', atom, after='allow')
    cmo2NonBondLineRe = re.compile(cmo2NonBondLine)

    ####AntiBonding Regex####
    cmo2AntiBondLine = namedRe("Coeff", reFloat, before='allow')
    cmo2AntiBondLine += namedRe('Atom1', atom, after='none') + r"\-"
    cmo2AntiBondLine += namedRe('Atom2', atom, before='allow', after='allow') + r"\*"
    cmo2AntiBondLineRe = re.compile(cmo2AntiBondLine)

    ####MOLine Regex####
    cmo2MOLine = namedRe("MO", r"\d+", after="none") + r"\(" + r"\w+" + r"\)"
    cmo2MOLineRe = re.compile(cmo2MOLine)

    ####TotalLine Regex####
    cmo2TotalLine = namedRe("Bonding", reFloat, after='none') + r"\(" + r"\w" + r"\)"
    cmo2TotalLine += namedRe("NonBonding", reFloat, before='require', after='none') + r"\(" + r"\w" + r"\)"
    cmo2TotalLine += namedRe("AntiBonding", reFloat, before='require', after='none') +r"\(" + r"\w" + r"\)"
    cmo2TotalLine += r"\s+\w+"
    cmo2TotalLineRe = re.compile(cmo2TotalLine)

    result={}
    currentMO = None
    for line in cmo2:
        if cmo2MOLineRe.match(line):
            correct = cmo2MOLineRe.match(line).groupdict()
            currentMO = correct["MO"]
            result[currentMO] = dict()
            result[currentMO]["Bonding"] = []
            result[currentMO]["NonBonding"] = []
            result[currentMO]["AntiBonding"] = []
            result[currentMO]["Total"] = []
            treatline = line.split()
            newline = ""
            for elem in treatline[1:]:
                newline += elem + " "
            if cmo2BondLineRe.match(newline):
                info = cmo2BondLineRe.match(newline).groupdict()
                result[currentMO]["Bonding"].append(info)
            if cmo2AntiBondLineRe.match(newline):
                info = cmo2AntiBondLineRe.match(newline).groupdict()
                result[currentMO]["AntiBonding"].append(info)
            if cmo2NonBondLineRe.match(newline):
                info = cmo2NonBondLineRe.match(newline).groupdict()
                result[currentMO]["NonBonding"].append(info)
        elif cmo2BondLineRe.match(line):
            info = cmo2BondLineRe.match(line).groupdict()
            result[currentMO]["Bonding"].append(info)
        elif cmo2AntiBondLineRe.match(line):
            info = cmo2AntiBondLineRe.match(line).groupdict()
            result[currentMO]["AntiBonding"].append(info)
        elif cmo2NonBondLineRe.match(line):
            info = cmo2NonBondLineRe.match(line).groupdict()
            result[currentMO]["NonBonding"].append(info)
        elif cmo2TotalLineRe.match(line):
            info = cmo2TotalLineRe.match(line).groupdict()
            result[currentMO]["Total"].append(info)
        else:
            print("parseCMO2: WARNING line not recognized "+ line)

    return result

