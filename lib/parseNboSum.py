import re
from .basicReadingFunctions import namedRe, find, findExact, extractTab, replacingDigit
import pandas as pd

def parseNonBond(file):
    reFloat = r"-?\d+\.\d+"
    atomicType = r"\d*[A-Z][A-Z]?[a-z]?\*?"
    atom = r"[A-Z][a-z]?"
    nboSumLine = namedRe('NBO',r"\d+",before='allow')
    nboSumLine += namedRe('Type',atomicType)
    nboSumLine += r"\(" + namedRe('AO',r"\d+",before='allow',after='none') + r"\)"
    nboSumLine += namedRe('Atom',atom,before='allow')
    nboSumLine += namedRe('Loc',r"\d+",after='allow')
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
    atom = r"[A-Z][a-z]?"
    atom2 = r"[A-Z][a-z]?"
    nboSumLine = namedRe('NBO',r"\d+",before='allow')
    nboSumLine += namedRe('Type',atomicType,after='allow')
    nboSumLine += r"\(" + namedRe('AO',r"\d+",before='allow',after='none') + r"\)"
    nboSumLine += namedRe('Atom1',atom,before='allow',after='allow')
    nboSumLine += namedRe('Loc1',r"\d+",after='none')+r"\-"
    nboSumLine += namedRe('Atom2',atom2,before='allow',after='allow')
    nboSumLine += namedRe('Loc2',r"\d+")
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
    tmpFile = [] # This part cleans out the digit with '.' at end so it is easier to parse
    for line in file:
        if len(line) > 0: 
            line = replacingDigit(0,line)
            tmpFile.append(line)
    posCR = find("CR",tmpFile)      #Finding the location of different type
    posLP = find("LP",tmpFile)
    posLV = find("LV",tmpFile)
    posBD = findExact("BD",tmpFile)
    posBDS = find("BD*(",tmpFile)
    pos3C = findExact("3C",tmpFile)
    pos3Cn = find("3Cn(",tmpFile)
    pos3Cs = find("3C*(",tmpFile)
    tabCR = extractTab(tmpFile,posCR)   #Putting the same type in one location
    tabLP = extractTab(tmpFile,posLP)
    tabLV = extractTab(tmpFile,posLV)
    tabBD = extractTab(tmpFile,posBD)
    tabBDS = extractTab(tmpFile,posBDS)
    tab3C = extractTab(tmpFile,pos3C)
    tab3Cn = extractTab(tmpFile,pos3Cn)
    tab3Cs = extractTab(tmpFile,pos3Cs)
    tabCRF = parseNonBond(tabCR)    #Parsing using regular expression as shown on top
    tabLPF = parseNonBond(tabLP)
    tabLVF = parseNonBond(tabLV)
    tabBDF = parseBond(tabBD)
    tabBDSF = parseBond(tabBDS)
    tab3CF = parseNonBond(tab3C)
    tab3CnF = parseNonBond(tab3Cn)
    tab3CsF = parseNonBond(tab3Cs)
    tableCRF = pd.DataFrame(tabCRF)     #Change the data type so they are not all strings or objects
    try:
        tableCRF[['NBO', 'AO','Loc']] = tableCRF[['NBO', 'AO','Loc']].astype(int) #adopted from https://stackoverflow.com/questions/15891038/change-data-type-of-columns-in-pandas
        tableCRF[['Type','Atom']] = tableCRF[['Type','Atom']].astype(str)
        tableCRF[['Occupancy','Energy']] = tableCRF[['Occupancy','Energy']].astype(float)
    except:
        print("TableCRF is empty or in wrong format.")
    tableLPF = pd.DataFrame(tabLPF)
    try:
        tableLPF[['NBO', 'AO','Loc']] = tableLPF[['NBO', 'AO','Loc']].astype(int) 
        tableLPF[['Type','Atom']] = tableLPF[['Type','Atom']].astype(str)
        tableLPF[['Occupancy','Energy']] = tableLPF[['Occupancy','Energy']].astype(float)
    except:
        print("TableLPF is empty or in wrong format.")
    tableLVF = pd.DataFrame(tabLVF)
    try:
        tableLVF[['NBO', 'AO','Loc']] = tableLVF[['NBO', 'AO','Loc']].astype(int) 
        tableLVF[['Type','Atom']] = tableLVF[['Type','Atom']].astype(str)
        tableLVF[['Occupancy','Energy']] = tableLVF[['Occupancy','Energy']].astype(float)
    except:
        print("TableLVF is empty or in wrong format.")
    tableBDF = pd.DataFrame(tabBDF)
    try:
        tableBDF[['NBO', 'AO','Loc1','Loc2']] = tableBDF[['NBO', 'AO','Loc1','Loc2']].astype(int) 
        tableBDF[['Type','Atom1','Atom2']] = tableBDF[['Type','Atom1','Atom2']].astype(str)
        tableBDF[['Occupancy','Energy']] = tableBDF[['Occupancy','Energy']].astype(float)
    except:
        print("TableBDF is empty or in wrong format.")
    tableBDSF = pd.DataFrame(tabBDSF)
    try:
        tableBDSF[['NBO', 'AO','Loc1','Loc2']] = tableBDSF[['NBO', 'AO','Loc1','Loc2']].astype(int) 
        tableBDSF[['Type','Atom1','Atom2']] = tableBDSF[['Type','Atom1','Atom2']].astype(str)
        tableBDSF[['Occupancy','Energy']] = tableBDSF[['Occupancy','Energy']].astype(float)
    except:
        print("TableBDSF is empty or in wrong format.")
    table3CF = pd.DataFrame(tab3CF)
    try:
        table3CF[['NBO', 'AO','Loc']] = table3CF[['NBO', 'AO','Loc']].astype(int) 
        table3CF[['Type','Atom']] = table3CF[['Type','Atom']].astype(str)
        table3CF[['Occupancy','Energy']] = table3CF[['Occupancy','Energy']].astype(float)
    except:
        print("Table3CF is empty or in wrong format.")
    table3CnF = pd.DataFrame(tab3CnF)
    try:
        table3CnF[['NBO', 'AO','Loc']] = table3CnF[['NBO', 'AO','Loc']].astype(int) 
        table3CnF[['Type','Atom']] = table3CnF[['Type','Atom']].astype(str)
        table3CnF[['Occupancy','Energy']] = table3CnF[['Occupancy','Energy']].astype(float)
    except:
        print("Table3CnF is empty or in wrong format.")
    table3CsF = pd.DataFrame(tab3CsF)
    try:
        table3CsF[['NBO', 'AO','Loc']] = table3CsF[['NBO', 'AO','Loc']].astype(int) 
        table3CsF[['Type','Atom']] = table3CsF[['Type','Atom']].astype(str)
        table3CsF[['Occupancy','Energy']] = table3CsF[['Occupancy','Energy']].astype(float)
    except:
        print("Table3CsF is empty or in wrong format.")
    return (tableCRF,tableLPF,tableLVF,tableBDF,tableBDSF,table3CF,table3CnF,table3CsF)