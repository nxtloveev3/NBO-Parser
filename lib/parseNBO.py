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
    return info


reFloat = r"-?\d+\.\d+"
nonBondType = r"[A-Z]{2,}"
atom = r"[A-Z][a-z]?\s*\d+"

####NBO NonBonding Regex####
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
coeffN += namedRe('Atom', r'[A-Z][a-z]*', before='allow', after='allow')
coeffN += namedRe('Loc', r'\d\d*', before='allow', after='allow')
coeffNRe = re.compile(coeffN)

def parseNBO(text, verbose=False):

    start = text.index("(Occupancy)   Bond orbital / Coefficients / Hybrids")
    text = text[start:]
    
    result = {}
    currentIndex = None
    for line in text:
        if "BD*" in line:
            line = line.replace("BD*", "ABB") #ABB stands for antibonding bonds
        parsed = False
        for r in [nboNonBondLineRe, nboBondLineRe, nbo3CBondLineRe]:
            if r.match(line):
                info = fix_info(r.match(line).groupdict())
                currentIndex = int(info["Index"])
                del info["Index"]
                result[currentIndex] = info
                parsed = True
                break
        if not parsed and coeffNRe.search(line):
            info = coeffNRe.search(line).groupdict()
            loc = int(info['Loc'])
            for i in range(1,4):
                if "Loc"+str(i) in result[currentIndex]:
                    if loc == result[currentIndex]["Loc"+str(i)] and "Coeff"+str(i) not in result[currentIndex]:
                        result[currentIndex]["Coeff"+str(i)] = float(info['coeff'])
        elif not parsed:
            if verbose:
                print("Ignoring this line:", line)
    return result