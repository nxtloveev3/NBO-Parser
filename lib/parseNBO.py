from .basicReadingFunctions import namedRe, find, fix_badatom
import re

def fix_info(info):
    for key in list(info.keys()):
        if key in ['Index', 'BondOrbital']:
            info[key] = int(info[key])
        elif key in ['Occupancy']:
            info[key] = float(info[key])
        elif key[:4]=='Atom':
            atom, i = fix_badatom(info[key])
            info[key] = atom
            info['Loc'+key[4:]] = int(i)

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
    nboNonBondLine += r"\(" + namedRe("BondOrbital", r"\d+", before='allow', after='none') + r"\)"
    nboNonBondLine += namedRe("Atom1", atom, before='allow')
    nboNonBondLine += namedRe("Hybrids", hybrids,after='allow')
    nboNonBondLineRe = re.compile(nboNonBondLine)

    ####Bonding Regex####
    bondType = r"[A-Z][A-Z][A-Z]?"
    nboBondLine = namedRe("Index", r"\d+",after='none') + r"\." + r"\s+"
    nboBondLine += r"\(" + namedRe("Occupancy", reFloat, before='allow', after='none') + r"\)"
    nboBondLine += namedRe("Type", bondType, before='require',after='allow')
    nboBondLine += r"\(" + namedRe("BondOrbital", r"\d+", before='allow', after='none') + r"\)"
    nboBondLine += namedRe("Atom1", atom, before='allow', after='none') + r"\-"
    nboBondLine += namedRe("Atom2", atom, before='allow', after='allow')
    nboBondLineRe = re.compile(nboBondLine)

    ####3C bond Regex####
    nbo3CBondLine = namedRe("Index", r"\d+",after='none') + r"\." + r"\s+"
    nbo3CBondLine += r"\(" + namedRe("Occupancy", reFloat, before='allow', after='none') + r"\)"
    nbo3CBondLine += namedRe("Type", r'3C[a-z]*\**', before='require',after='allow')
    nbo3CBondLine += r"\(" + namedRe("BondOrbital", r"\d+", before='allow', after='none') + r"\)"
    nbo3CBondLine += namedRe("Atom1", atom, before='allow', after='none') + r"\-"
    nbo3CBondLine += namedRe("Atom2", atom, before='allow', after='none') + r"\-"
    nbo3CBondLine += namedRe("Atom3", atom, before='allow', after='allow')
    nbo3CBondLineRe = re.compile(nbo3CBondLine)

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
            fix_info(info)
            currentIndex = int(info["Index"])
            del info["Index"]
            if verbose:
                print(currentIndex)
            result[currentIndex] = info
        elif nboBondLineRe.match(line):
            info = nboBondLineRe.match(line).groupdict()
            fix_info(info)
            currentIndex = int(info["Index"])
            del info["Index"]
            result[currentIndex] = info
        elif nbo3CBondLineRe.match(line):
            info = nbo3CBondLineRe.match(line).groupdict()
            fix_info(info)
            currentIndex = int(info["Index"])
            del info["Index"]
            result[currentIndex] = info
        elif coeffNRe.search(line):
            info = coeffNRe.search(line).groupdict()
            loc = int(info['Loc'])
            if loc == result[currentIndex]["Loc1"] and "Coeff1" not in result[currentIndex]:
                result[currentIndex]["Coeff1"] = float(info['coeff'])
            elif loc == result[currentIndex]["Loc2"] and "Coeff2" not in result[currentIndex]:
                result[currentIndex]["Coeff2"] = float(info['coeff'])
            elif loc == result[currentIndex]["Loc3"] and "Coeff3" not in result[currentIndex]:
                result[currentIndex]["Coeff3"] = float(info['coeff'])
        else:
            if verbose:
                print("Ignoring this line:", line)
    return result