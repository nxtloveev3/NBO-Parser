import re
from .basicReadingFunctions import namedRe, find, findExact, extractTab, replacingDigit
import pandas as pd

def parseNonBond(file):
    reFloat = r"-?\d+\.\d+"
    atomicType = r"\d*[A-Z][A-Z]?[a-z]?\*?"
    atom = r"[A-Z][a-z]?\s*\d+"
    nboSumLine = namedRe('NBO',r"\d+",before='allow')
    nboSumLine += namedRe('Type',atomicType)
    nboSumLine += r"\(" + namedRe('AO',r"\d+",before='allow',after='none') + r"\)"
    nboSumLine += namedRe('Atom',atom,before='allow')
    nboSumLine += namedRe('Occupancy', reFloat)
    nboSumLine += namedRe('Energy', reFloat,after='allow')
    nboLineRe = re.compile(nboSumLine)

    result = []
    for line in file:
        correct = nboLineRe.search(line)
        if correct:
            result.append(correct.groupdict())
        else:
            print('IGNORE: ', line)
    return result

def parseBond(file):
    reFloat = r"-?\d+\.\d+"
    atomicType = r"\d*[A-Z][A-Z]?[a-z]?\*?"
    atom = r"[A-Z][a-z]?\s*\d+\-\s*[A-Z][a-z]?\s*\d+"
    nboSumLine = namedRe('NBO',r"\d+\.",before='allow')
    nboSumLine += namedRe('Type',atomicType,after='allow')
    nboSumLine += r"\(" + namedRe('AO',r"\d+",before='allow',after='none') + r"\)"
    nboSumLine += namedRe('Atom',atom,before='allow')
    nboSumLine += namedRe('Occupancy', reFloat)
    nboSumLine += namedRe('Energy', reFloat,after='allow')
    nboLineRe = re.compile(nboSumLine)

    result = []
    for line in file:
        correct = nboLineRe.search(line)
        if correct:
            result.append(correct.groupdict())
        else:
            print('IGNORE: ', line)
    return result

def parseNboSum(file):
    tmpFile = []
    for line in file:
        if len(line) > 0: 
            line = replacingDigit(0,line)
            tmpFile.append(line)
    posCR = find("CR",tmpFile)
    posLP = find("LP",tmpFile)
    posLV = find("LV",tmpFile)
    posBD = findExact("BD",tmpFile)
    posBDS = find("BD*(",tmpFile)
    pos3C = findExact("3C",tmpFile)
    pos3Cn = find("3Cn(",tmpFile)
    pos3Cs = find("3C*(",tmpFile)
    tabCR = extractTab(tmpFile,posCR)
    tabLP = extractTab(tmpFile,posLP)
    tabLV = extractTab(tmpFile,posLV)
    tabBD = extractTab(tmpFile,posBD)
    tabBDS = extractTab(tmpFile,posBDS)
    tab3C = extractTab(tmpFile,pos3C)
    tab3Cn = extractTab(tmpFile,pos3Cn)
    tab3Cs = extractTab(tmpFile,pos3Cs)
    tabCRF = parseNonBond(tabCR)
    tabLPF = parseNonBond(tabLP)
    tabLVF = parseNonBond(tabLV)
    tabBDF = parseBond(tabBD)
    tabBDSF = parseBond(tabBDS)
    tab3CF = parseNonBond(tab3C)
    tab3CnF = parseNonBond(tab3Cn)
    tab3CsF = parseNonBond(tab3Cs)
    tableCRF = pd.DataFrame(tabCRF)
    print(tableCRF)
    tableLPF = pd.DataFrame(tabLPF)
    tableLVF = pd.DataFrame(tabLVF)
    tableBDF = pd.DataFrame(tabBDF)
    tableBDSF = pd.DataFrame(tabBDSF)
    table3CF = pd.DataFrame(tab3CF)
    table3CnF = pd.DataFrame(tab3CnF)
    table3CsF = pd.DataFrame(tab3CsF)
    convertDict = {'NBO': int, 
                   'Type': str,
                   'AO': int,
                   'Atom':str,
                   'Occupancy': float,
                   'Energy':float} 
    tableCRF = tableCRF.astype(convertDict)
    tableLPF = tableLPF.astype(convertDict)
    tableLVF = tableLVF.astype(convertDict)
    tableBDF = tableBDF.astype(convertDict)
    tableBDSF = tableBDSF.astype(convertDict)
    table3CF = table3CF.astype(convertDict)
    table3CnF = table3CnF.astype(convertDict)
    table3CsF = table3CsF.astype(convertDict)
    return (tableCRF,tableLPF,tableLVF,tableBDF,tableBDSF,table3CF,table3CnF,table3CsF)