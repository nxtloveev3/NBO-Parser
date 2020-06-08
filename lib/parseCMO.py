import re
from .basicReadingFunctions import namedRe, find, findExact, extractTab
import pandas as pd

def parseCMON(file):
    loc = find("Molecular Orbital Atom-Atom Bonding Character",file)
    cmon = file[:loc[0]]
    pos = find("MO",cmon)
    spp = []
    for num in pos:
        spp.append(cmon[num+1])
        spp.append(cmon[num+2])
    posCR = find("CR",spp)      #Finding the location of different type
    posLP = find("LP",spp)
    posLV = find("LV",spp)
    posBD = findExact("BD",spp)
    posBDS = find("BD*(",spp)
    pos3C = findExact("3C",spp)
    pos3Cn = find("3Cn(",spp)
    pos3Cs = find("3C*(",spp)
    tabCR = extractTab(spp,posCR)   #Putting the same type in one location
    tabLP = extractTab(spp,posLP)
    tabLV = extractTab(spp,posLV)
    tabBD = extractTab(spp,posBD)
    tabBDS = extractTab(spp,posBDS)
    tab3C = extractTab(spp,pos3C)
    tab3Cn = extractTab(spp,pos3Cn)
    tab3Cs = extractTab(spp,pos3Cs)
    tabCRF = parseNonBond(tabCR)    #Parsing using regular expression as shown on top
    tabLPF = parseNonBond(tabLP)
    tabLVF = parseNonBond(tabLV)
    tabBDF = parseBond(tabBD)
    tabBDSF = parseBond(tabBDS)
    tab3CF = parseNonBond(tab3C)
    tab3CnF = parseNonBond(tab3Cn)
    tab3CsF = parseNonBond(tab3Cs)
    tableCRF = pd.DataFrame(tabCRF) 
    print(tabCRF) 


def parseNonBond(file):
    reFloat = r"-?\d+\.\d+"
    atomicType = r"\d*[A-Z][A-Z]?[a-z]?"
    typelower = r"\d*[a-z][a-z]?"
    atom = r"[A-Z][a-z]?"
    cmonLine = namedRe('Energy',reFloat,before='allow',after='none')
    cmonLine += r"\*\[" + namedRe('NBO',r"\d+",after='none') + r"\]\:"
    cmonLine += namedRe('Type',atomicType,before='allow',after='allow')
    cmonLine += r"\(" + namedRe('QN',r"\d+",after='none')+r"\)"
    cmonLine += namedRe('Atom',atom,before='allow',after='allow')
    cmonLine += namedRe('Loc', r"\d+",before='none',after='allow')
    cmonLine += r"\(" + namedRe('type',typelower,before='none',after='none') + r"\)"
    cmonLineRe = re.compile(cmonLine)
    print(cmonLineRe)

    result = []
    for line in file:
        correct = cmonLineRe.search(line)
        if correct:
            result.append(correct.groupdict())
        else:
            print('IGNORE: ', line)
    return result

def parseBond(file):
    reFloat = r"-?\d+\.\d+"
    atomicType = r"\d?[A-Z][A-Z]?[a-z]?\*?"
    atom = r"[A-Z][a-z]?"
    cmonLine = namedRe('Energy',reFloat,before='allow',after='none')
    cmonLine += r"\*\[" + namedRe('NBO',r"\d+",after='none') + r"\]\:"
    cmonLine += namedRe('Type',atomicType,before='allow',after='allow')
    cmonLine += r"\(" + namedRe('QN',r"\d+",after='none')+r"\)"
    cmonLine += namedRe('Atom1',atom,before='allow',after='none')
    cmonLine += namedRe('Loc1', r"\d+",after='none') + r"\-"
    cmonLine += namedRe('Atom2',atom,before='allow',after='none')
    cmonLine += namedRe('Loc2', r"\d+",after='allow')
    cmonLineRe = re.compile(cmonLine)

    result = []
    for line in file:
        correct = cmonLineRe.search(line)
        if correct:
            result.append(correct.groupdict())
        else:
            print('IGNORE: ', line)
    return result