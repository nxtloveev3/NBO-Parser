import re
from .basicReadingFunctions import namedRe, fix_badatom
import pandas as pd

def fix_info(info):
    for key in list(info.keys()):
        if key in ['NBO', 'Index', 'Loc', 'Loc1', 'Loc2', 'Loc3']:
            try:
                info[key] = int(info[key])
            except:
                pass
        elif key in ['Coef', 'Coeff', 'Bonding', 'NonBonding', 'AntiBonding']:
            info[key] = float(info[key])
        elif key[:4]=='Atom' and 'Loc'+key[4:] not in info:
            atom, i = fix_badatom(info[key])
            info[key] = atom
            info['Loc'+key[4:]] = int(i)
    if info['Type'] == '3C*':
        info['Type'] = '3Cs'
    return info

def parseCMON(text, verbose=False):

    try:
        loc = text.index("Molecular Orbital Atom-Atom Bonding Character")
        text = text[:loc]
    except:
        pass

    reFloat = r"-?\d+\.\d+"
    atom = r"[A-Z][a-z]?"
    loc = r"\d?\**"
    ####NonBond Regex####
    nonBondType = r"\d*[A-Z][A-Z]?[a-z]?"
    cmonLine = namedRe("Coef", reFloat, before="allow", after='none')
    cmonLine += r"\*\[" + namedRe("NBO", r"\d+", before="allow", after='none') + r"\]\:"
    cmonLine += namedRe("Type", nonBondType, before="require")
    cmonLine += r"\(" + namedRe("Index", r"\d+", before="allow", after='none') + r"\)"
    cmonLine += namedRe("Atom", atom, before="allow", after='allow')
    cmonLine += namedRe("Loc", loc, before="none", after='none')
    cmonLineRe = re.compile(cmonLine)

    ####Bond Regex####
    bondType = r"[A-Z]{2,}"
    cmonBondLine = namedRe("Coef", reFloat, before="allow", after='none')
    cmonBondLine += r"\*\[" + namedRe("NBO", r"\d+", before='allow', after='none') + r"\]\:"
    cmonBondLine += namedRe("Type", bondType, before="require", after='allow')
    cmonBondLine += r"\(" + namedRe("Index", r"\d+",before="allow", after='none') + r"\)"
    cmonBondLine += namedRe("Atom1", atom, before="allow", after='allow')
    cmonBondLine += namedRe("Loc1", loc, before="allow", after='none')
    cmonBondLine += r"\-" + namedRe("Atom2", atom, before='allow', after='none')
    cmonBondLine += namedRe("Loc2", loc, before='allow', after='none') + r"\*?"
    cmonBondLineRe = re.compile(cmonBondLine)

    ####3C Bond Regex####
    cmon3CBondLine = namedRe("Coef", reFloat, before="allow", after='none')
    cmon3CBondLine += r"\*\[" + namedRe("NBO", r"\d+", before='allow', after='none') + r"\]\:"
    cmon3CBondLine += namedRe("Type",  r'3C[a-z]*\**', before="require", after='allow')
    cmon3CBondLine += r"\(" + namedRe("Index", r"\d+",before="allow", after='none') + r"\)"
    cmon3CBondLine += namedRe("Loc1", loc, before="allow", after='allow')+r'\-'
    cmon3CBondLine += namedRe("Loc2", loc, before='allow', after='allow')+r'\-'
    cmon3CBondLine += namedRe("Loc3", loc, before='allow', after='none')
    cmon3CBondLineRe = re.compile(cmon3CBondLine)

    ####MO Regex####
    cmonMOLine = namedRe("Label1", r"[M][O]", before='allow')
    cmonMOLine += namedRe("MO", r"\d+")
    cmonMOLine += r"\(" + namedRe("Type", r"\w+", after='none') + r"\)\:"
    cmonMOLine += namedRe("Label3", r"\w+\s+\w+", before="require")
    cmonMOLine += r"\=" + namedRe("Energy", reFloat, before="require")
    cmonMOLine += namedRe("Unit", r"[a]\.[u]\.",after='allow')
    cmonMOLineRe = re.compile(cmonMOLine)

    result = {}
    currentMO = None
    HOMO = 0
    for line in text:
        if "BD*" in line:
            line = line.replace("BD*", "ABB") #ABB stands for antibonding 
        if cmonMOLineRe.search(line):
            correct = cmonMOLineRe.search(line).groupdict()
            currentMO = int(correct["MO"])
            energy = float(correct["Energy"])
            typ = correct['Type']
            if HOMO == 0 and typ == 'vir':
                HOMO = currentMO - 1
            result[currentMO] = dict()
            result[currentMO]["Energy"] = energy
            result[currentMO]["Type"] = typ
            result[currentMO]["wf"] = []
        elif cmon3CBondLineRe.search(line):
            result[currentMO]["wf"].append(fix_info(cmon3CBondLineRe.search(line).groupdict()))
        elif cmonBondLineRe.search(line):
            result[currentMO]["wf"].append(fix_info(cmonBondLineRe.search(line).groupdict()))
        elif cmonLineRe.search(line):
            result[currentMO]["wf"].append(fix_info(cmonLineRe.search(line).groupdict()))
        else:
            if verbose:
                print("parseCMON WARNING: line not recognized "+ line)

    return result, HOMO


def parseCMO2(text, verbose=False):
    loc = text.index("Molecular Orbital Atom-Atom Bonding Character")
    text = text[loc:]
    reFloat = r"-?\d+\.\d+"
    # ####Bonding Regex####
    # atom = r"[A-Z][a-z]?\s*\d+"
    # cmo2BondLine = r"\b" + namedRe("Coeff", reFloat, before='allow', after='allow')
    # cmo2BondLine += namedRe('Atom1', atom, after='none') + r"\-"
    # cmo2BondLine += namedRe('Atom2', atom, before='allow', after='allow') + r"\b"
    # cmo2BondLineRe = re.compile(cmo2BondLine)

    # ####NonBonding Regex####
    # cmo2NonBondLine = namedRe("Coeff", reFloat, before='require', after='require')
    # cmo2NonBondLine += namedRe('Atom', atom, before='allow', after='allow')
    # cmo2NonBondLineRe = re.compile(cmo2NonBondLine)

    # ####AntiBonding Regex####
    # cmo2AntiBondLine = namedRe("Coeff", reFloat, before='allow')
    # cmo2AntiBondLine += namedRe('Atom1', atom, after='none') + r"\-"
    # cmo2AntiBondLine += namedRe('Atom2', atom, before='allow', after='allow') + r"\*"
    # cmo2AntiBondLineRe = re.compile(cmo2AntiBondLine)

    ####MOLine Regex####
    cmo2MOLine = namedRe("MO", r"\d+", after="none") + namedRe("Type", r"\(" + r"\w+" + r"\)", after="none")
    cmo2MOLineRe = re.compile(cmo2MOLine)

    ####TotalLine Regex####
    cmo2TotalLine = namedRe("Bonding", reFloat, after='none') + r"\(" + r"\w" + r"\)"
    cmo2TotalLine += namedRe("NonBonding", reFloat, before='require', after='none') + r"\(" + r"\w" + r"\)"
    cmo2TotalLine += namedRe("AntiBonding", reFloat, before='require', after='none') +r"\(" + r"\w" + r"\)"
    cmo2TotalLine += r"\s+\w+"
    cmo2TotalLineRe = re.compile(cmo2TotalLine)

    result = {}
    currentMO = None
    SOMO = 0
    for line in text:
        if cmo2MOLineRe.match(line):
            correct = cmo2MOLineRe.match(line).groupdict()
            currentMO = int(correct["MO"])
            typ = correct['Type']
            if SOMO == 0 and typ == '(v)':
                SOMO = currentMO - 1
            result[currentMO] = {'Type': typ} 
        elif cmo2TotalLineRe.search(line):
            info = fix_info(cmo2TotalLineRe.search(line).groupdict())
            for key in info:
                if key not in result[currentMO]:
                    result[currentMO][key] = info[key]
        elif verbose:
            print("parseCMO2 WARNING: line not recognized "+ line)
    return result, SOMO