from .basicReadingFunctions import namedRe, find
import re

def parseNBO(file, verbose=False):
    start = find("(Occupancy)   Bond orbital / Coefficients / Hybrids", file)
    nbo = file[start[0]:]
    tmpFile = []
    for line in nbo:
        newline = line
        if "BD*" in line:
            newline = line.replace("BD*","ABB") #ABB stand for antibonding bonds
        tmpFile.append(newline)
    
    ####NBO NonBonding Regex####
    reFloat = r"-?\d+\.\d+"
    nonBondType = r"[A-Z][A-Z]"
    atom = r"[A-Z][a-z]?\s*\d+"
    hybrids = r"[s]\(\s*\d+\.\d+\%\)([p]\s*\d+\.\d+\(\s*\d+\.\d+\%\))?([d]\s*\d+\.\d+\(\s*\d+\.\d+\%\))?"
    nboNonBondLine = namedRe("Index", r"\d+",after='none') + r"\." + r"\s+"
    nboNonBondLine += r"\(" + namedRe("Occupancy", reFloat, before='allow', after='none') + r"\)"
    nboNonBondLine += namedRe("Type", nonBondType, before='require')
    nboNonBondLine += r"\(" + namedRe("BondOrbitals", r"\d+", before='allow', after='none') + r"\)"
    nboNonBondLine += namedRe("Atom", atom, before='allow')
    nboNonBondLine += namedRe("Hybrids", hybrids,after='allow')
    nboNonBondLineRe = re.compile(nboNonBondLine)

    ####Bonding Regex####
    bondType = r"[A-Z][A-Z][A-Z]?"
    nboBondLine = namedRe("Index", r"\d+",after='none') + r"\." + r"\s+"
    nboBondLine += r"\(" + namedRe("Occupancy", reFloat, before='allow', after='none') + r"\)"
    nboBondLine += namedRe("Type", bondType, before='require',after='allow')
    nboBondLine += r"\(" + namedRe("BondOrbitals", r"\d+", before='allow', after='none') + r"\)"
    nboBondLine += namedRe("Atom1", atom, before='allow', after='none') + r"\-"
    nboBondLine += namedRe("Atom2", atom, before='allow', after='allow')
    nboBondLineRe = re.compile(nboBondLine)

    ####Hybrids(bond) Regex####
    nboBondHybridsLine = namedRe("Hybrids", hybrids, before='allow', after='allow')
    nboBondHybridsLineRe = re.compile(nboBondHybridsLine)

    
    def fix_badatom(atom):
        if ' ' in atom:
            return atom.split()
        a = ''
        i = ''
        for chr_ in atom:
            if chr_.isalpha():
                a += chr_
            else:
                i += chr_
        return [a, i]
        
    result = {}
    currentIndex = None
    for line in tmpFile:
        if nboNonBondLineRe.match(line):
            info = nboNonBondLineRe.match(line).groupdict()
            currentIndex = info["Index"]
            if verbose:
                print(currentIndex)
            occupancy = info["Occupancy"]
            atomictype = info["Type"]
            bondOrbital = info["BondOrbitals"]
            atom = fix_badatom(info["Atom"])
            hybrids = info["Hybrids"]
            result[currentIndex] = dict()
            result[currentIndex]["Occupancy"] = float(occupancy)
            result[currentIndex]["Type"] = atomictype
            result[currentIndex]["BondOrbital"] = int(bondOrbital)
            result[currentIndex]["Atom"] = atom[0]
            result[currentIndex]["Atom_ind"] = int(atom[1])
            result[currentIndex]["Hybrids"] = hybrids
        elif nboBondLineRe.match(line):
            info = nboBondLineRe.match(line).groupdict()
            currentIndex = info["Index"]
            occupancy = info["Occupancy"]
            atomictype = info["Type"]
            bondOrbital = info["BondOrbitals"]
            atom1 = fix_badatom(info["Atom1"])
            atom2 = fix_badatom(info["Atom2"])
            result[currentIndex] = dict()
            result[currentIndex]["Occupancy"] = float(occupancy)
            result[currentIndex]["Type"] = atomictype
            result[currentIndex]["BondOrbital"] = int(bondOrbital)
            result[currentIndex]["Atom1"] = atom1[0]
            result[currentIndex]["Atom1_ind"] = int(atom1[1])
            result[currentIndex]["Atom2"] = atom2[0]
            result[currentIndex]["Atom2_ind"] = int(atom2[1])
        elif nboBondHybridsLineRe.search(line):
            info = nboBondHybridsLineRe.search(line).groupdict()
            hybrids = info["Hybrids"]
            result[currentIndex]["Hybrids"] = hybrids
        else:
            if verbose:
                print("Ignoring this line:", line)
    return result