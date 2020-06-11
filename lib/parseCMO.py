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
    unparsed = []
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
            unparsed.append(line)
            print("parseCMON: WARNING line not recognized "+ line)

    return result
    
    