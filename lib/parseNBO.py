from .basicReadingFunctions import namedRe, find, fix_badatom
import re

def parseNBO(text, verbose=False):
    start = text.index("(Occupancy)   Bond orbital / Coefficients / Hybrids")
    text = text[start:]
    
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

    ####Coeff of NAO Regex####
    coeffN = namedRe('coeff', r'0.\d\d*', before='allow', after='none')
    coeffN += r'\*'
    coeffN += namedRe('Atom', r'[A-Z][a-z]*', before='require', after='allow')
    coeffN += namedRe('Loc', r'\d\d*', before='allow', after='allow')
    coeffNRe = re.compile(coeffN)
        
    result = {}
    currentIndex = None
    for line in text:
        if "BD*" in line:
            line = line.replace("BD*", "ABB") #ABB stand for antibonding bonds
        if nboNonBondLineRe.match(line):
            info = nboNonBondLineRe.match(line).groupdict()
            currentIndex = int(info["Index"])
            if verbose:
                print(currentIndex)
            atom = fix_badatom(info["Atom"])
            result[currentIndex] = dict()
            result[currentIndex]["Occupancy"] = float(info["Occupancy"])
            result[currentIndex]["Type"] = info["Type"]
            result[currentIndex]["BondOrbital"] = int(info["BondOrbitals"])
            result[currentIndex]["Atom"] = atom[0]
            result[currentIndex]["Loc"] = int(atom[1])
            result[currentIndex]["Hybrids"] = info["Hybrids"]
        elif nboBondLineRe.match(line):
            info = nboBondLineRe.match(line).groupdict()
            currentIndex = int(info["Index"])
            atom1 = fix_badatom(info["Atom1"])
            atom2 = fix_badatom(info["Atom2"])
            result[currentIndex] = dict()
            result[currentIndex]["Occupancy"] = float(info["Occupancy"])
            result[currentIndex]["Type"] = info["Type"]
            result[currentIndex]["BondOrbital"] = int(info["BondOrbitals"])
            result[currentIndex]["Atom1"] = atom1[0]
            result[currentIndex]["Loc1"] = int(atom1[1])
            result[currentIndex]["Atom2"] = atom2[0]
            result[currentIndex]["Loc2"] = int(atom2[1])
        elif coeffNRe.search(line):
            info = coeffNRe.search(line).groupdict()
            loc = int(info['Loc'])
            if loc == result[currentIndex]["Loc1"] and "Coeff1" not in result[currentIndex]:
                result[currentIndex]["Coeff1"] = float(info['coeff'])
            elif loc == result[currentIndex]["Loc2"] and "Coeff2" not in result[currentIndex]:
                result[currentIndex]["Coeff2"] = float(info['coeff'])
        else:
            if verbose:
                print("Ignoring this line:", line)
    return result