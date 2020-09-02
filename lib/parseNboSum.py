import re
from .basicReadingFunctions import namedRe, find, findExact, extractTab, replacingDigit
import pandas as pd

def parseNonBond(file, verbose=False):
    reFloat = r"-?\d+\.\d+"
    bondType = r"\d*[A-Z][A-Z]?[a-z]?\*?"
    atom = r"[A-Z][a-z]?"
    nboSumLine = namedRe('NBO',r"\d+",before='allow')
    nboSumLine += namedRe('Type',bondType)
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
    bondType = r"\d*[A-Z][A-Z]?[a-z]?\*?"
    atom = r"[A-Z][a-z]?"
    nboSumLine = namedRe('NBO',r"\d+",before='allow')
    nboSumLine += namedRe('Type',bondType,after='allow')
    nboSumLine += r"\(" + namedRe('AO',r"\d+",before='allow',after='none') + r"\)"
    nboSumLine += namedRe('Atom1',atom,before='allow',after='allow')
    nboSumLine += namedRe('Loc1',r"\d+",after='none')+r"\-"
    nboSumLine += namedRe('Atom2',atom,before='allow',after='allow')
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

def parse3CBond(file, verbose=False):
    reFloat = r"-?\d+\.\d+"
    atom = r"[A-Z][a-z]?"
    nboSumLine = namedRe('NBO',r'\d+.',before='allow', after='allow')
    nboSumLine += namedRe('Type', r'3C[a-z]*\**', before='allow', after='allow')
    nboSumLine += r"\(" + namedRe('AO',r"\d+",before='allow',after='none') + r"\)"
    nboSumLine += namedRe('Atom1',atom,before='allow',after='allow')
    nboSumLine += namedRe('Loc1',r"\d+",after='none')+r"\-"
    nboSumLine += namedRe('Atom2',atom,before='allow',after='allow')
    nboSumLine += namedRe('Loc2',r"\d+",after='none')+r"\-"
    nboSumLine += namedRe('Atom3',atom,before='allow',after='allow')
    nboSumLine += namedRe('Loc3',r"\d+")
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

def parseNboSum(file, verbose=False, selected=["CR", "LP", "LV", "BD", "ABB", "3C", "3Cn", "3Cs"]):

    if type(file[0])==str:
        for i,line in enumerate(file):
            if len(line) > 0: 
                file[i] = replacingDigit(0,line)

    tables = {}

    def helper(df, list_1, list_2, list_3):
        df[list_1] = df[list_1].astype(int)
        df[list_2] = df[list_2].astype(str)
        df[list_3] = df[list_3].astype(float)
        return df
    
    for key in ['CR', 'LP', 'LV']:
        if key in selected:
            pos = find(key, file)
            table = pd.DataFrame(parseNonBond(extractTab(file, pos), verbose=verbose))
            try:
                table = helper(table, ['NBO', 'AO','Loc'], ['Type','Atom'], ['Occupancy','Energy'])
                tables[key] = table
            except:
                if verbose:
                    print("Table of "+key+" is empty or in wrong format.")
    
    for key in ['BD', 'ABB']:
        if key in selected:
            if key=='ABB':
                pos = find("BD*(", file)
            else:
                pos = findExact(key, file)
            tab = pd.DataFrame(parseBond(extractTab(file, pos), verbose=verbose))
            try:
                tab = helper(tab, ['NBO','AO','Loc1','Loc2'], ['Type','Atom1','Atom2'], ['Occupancy','Energy'])
                tables[key] = tab
            except:
                if verbose:
                    print("Table of "+key+" is empty or in wrong format.")

    for key in ['3C', '3Cn', '3Cs']:
        if key in selected:
            if key=='3Cn':
                pos = find("3Cn(", file)
            elif key=='3Cs':
                pos = find("3C*(", file)
            else:
                pos = findExact("3C", file)
            tab = pd.DataFrame(parse3CBond(extractTab(file, pos), verbose=verbose))
            try:
                # tab = helper(tab, ['NBO', 'AO','Loc'], ['Type','Atom'], ['Occupancy','Energy'])
                tables[key] = tab
            except:
                if verbose:
                    print("Table of "+key+" is empty or in wrong format.")

    return tables