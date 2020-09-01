import re
from .basicReadingFunctions import namedRe, find, findExact, extractTab, replacingDigit
import pandas as pd

def parseNonBond(file, verbose=False):
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
        elif verbose:
            print('IGNORE: ', line)
    return result

def parseBond(file, verbose=False):
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
        elif verbose:
            print('IGNORE: ', line)
    return result

def parseNboSum(file, verbose=False, selected_tabs=["CR", "LP", "LV", "BD", "ABB", "3C", "3Cn", "3Cs"]):
    tmpFile = [] # This part cleans out the digit with '.' at end so it is easier to parse
    tables = {}
    for line in file:
        if len(line) > 0: 
            tmpFile.append(replacingDigit(0,line))

    def helper(df, list_1, list_2, list_3):
        df[list_1] = df[list_1].astype(int)
        df[list_2] = df[list_2].astype(str)
        df[list_3] = df[list_3].astype(float)
        return df
    
    if "CR" in selected_tabs:
        posCR = find("CR",tmpFile)      #Finding the location of different types of nbo
        tabCR = extractTab(tmpFile,posCR)       #Putting the same type in one location
        tabCRF = parseNonBond(tabCR, verbose=verbose)    #Parsing using regular expression as shown on top
        tableCRF = pd.DataFrame(tabCRF)     #Change the data type so they are not all strings or objects
        del posCR, tabCR, tabCRF
        try:
            tableCRF = helper(tableCRF, ['NBO', 'AO','Loc'], ['Type','Atom'], ['Occupancy','Energy'])
            tables['CR'] = tableCRF
        except:
            if verbose:
                print("TableCRF is empty or in wrong format.")
    
    if "LP" in selected_tabs:
        posLP = find("LP",tmpFile)
        tabLP = extractTab(tmpFile,posLP)
        tabLPF = parseNonBond(tabLP, verbose=verbose)
        tableLPF = pd.DataFrame(tabLPF)
        del posLP, tabLP, tabLPF
        try:
            tableLPF = helper(tableLPF, ['NBO', 'AO','Loc'], ['Type','Atom'], ['Occupancy','Energy'])
            tables['LP'] = tableLPF
        except:
            if verbose:
                print("TableLPF is empty or in wrong format.")

    if "LV" in selected_tabs:
        posLV = find("LV",tmpFile)
        tabLV = extractTab(tmpFile,posLV)
        tabLVF = parseNonBond(tabLV, verbose=verbose)
        tableLVF = pd.DataFrame(tabLVF)
        del posLV, tabLV, tabLVF
        try:
            tableLVF = helper(tableLVF, ['NBO', 'AO','Loc'], ['Type','Atom'], ['Occupancy','Energy'])
            tables['LV'] = tableLVF
        except:
            if verbose:
                print("TableLVF is empty or in wrong format.")

    if "BD" in selected_tabs:
        posBD = findExact("BD",tmpFile)
        tabBD = extractTab(tmpFile,posBD)
        tabBDF = parseBond(tabBD, verbose=verbose)
        tableBDF = pd.DataFrame(tabBDF)
        del posBD, tabBD, tabBDF
        try:
            tableBDF = helper(tableBDF, ['NBO', 'AO','Loc1','Loc2'], ['Type','Atom1','Atom2'], ['Occupancy','Energy'])
            tables['BD'] = tableBDF
        except:
            if verbose:
                print("TableBDF is empty or in wrong format.")
    
    if "ABB" in selected_tabs:
        posABB = find("BD*(",tmpFile)
        tabABB = extractTab(tmpFile,posABB)
        tabABBF = parseBond(tabABB, verbose=verbose)
        tableABBF = pd.DataFrame(tabABBF)
        del posABB, tabABB, tabABBF
        try:
            tableABBF = helper(tableABBF, ['NBO', 'AO','Loc1','Loc2'], ['Type','Atom1','Atom2'], ['Occupancy','Energy'])
            tables['ABB'] = tableABBF
        except:
            if verbose:
                print("TableBDSF is empty or in wrong format.")

    if "3C" in selected_tabs:
        pos3C = findExact("3C",tmpFile)
        tab3C = extractTab(tmpFile,pos3C)
        tab3CF = parseNonBond(tab3C, verbose=verbose)
        table3CF = pd.DataFrame(tab3CF)
        del pos3C, tab3C, tab3CF
        try:
            table3CF = helper(table3CF, ['NBO', 'AO','Loc'], ['Type','Atom'], ['Occupancy','Energy'])
            tables['3C'] = table3CF
        except:
            if verbose:
                print("Table3CF is empty or in wrong format.")

    if "3Cn" in selected_tabs:
        pos3Cn = find("3Cn(",tmpFile)
        tab3Cn = extractTab(tmpFile,pos3Cn)
        tab3CnF = parseNonBond(tab3Cn, verbose=verbose)
        table3CnF = pd.DataFrame(tab3CnF)
        del pos3Cn, tab3Cn, tab3CnF
        try:
            table3CnF = helper(table3CnF, ['NBO', 'AO','Loc'], ['Type','Atom'], ['Occupancy','Energy'])
            tables['3Cn'] = table3CnF.copy()
            del table3CnF
        except:
            if verbose:
                print("Table3CnF is empty or in wrong format.")

    if "3Cs" in selected_tabs:
        pos3Cs = find("3C*(",tmpFile)
        tab3Cs = extractTab(tmpFile,pos3Cs)
        tab3CsF = parseNonBond(tab3Cs, verbose=verbose)
        table3CsF = pd.DataFrame(tab3CsF)
        del pos3Cs, tab3Cs, tab3CsF
        try:
            table3CsF = helper(table3CsF, ['NBO', 'AO','Loc'], ['Type','Atom'], ['Occupancy','Energy'])
            tables['3Cs'] = table3CsF
        except:
            if verbose:
                print("Table3CsF is empty or in wrong format.")
    return tables