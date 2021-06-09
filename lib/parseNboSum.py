import re
from .basicReadingFunctions import namedRe, find, extract
import pandas as pd

def fix_info(info):
    info['NBO'] = info['NBO'][:-1]
    return info

reFloat = r"-?\d+\.\d+"
bondType = r"\d*[A-Z][A-Z]?[a-z]?\*?"
atom = r"[A-Z][a-z]?"

def parseNonBond(file, verbose=False):

    nboSumLine = namedRe('NBO',r"\d+.",before='allow')
    nboSumLine += namedRe('Type',bondType)
    nboSumLine += r"\(" + namedRe('AO',r"\d+",before='allow',after='none') + r"\)"
    nboSumLine += namedRe('Atom',atom,before='allow', after='allow')
    nboSumLine += namedRe('Loc',r"\d+",after='allow')
    nboSumLine += namedRe('Occupancy', reFloat)
    nboSumLine += namedRe('Energy', reFloat,after='allow')
    nboLineRe = re.compile(nboSumLine)
    
    result = []
    for line in file:
        correct = nboLineRe.search(line)
        if correct:
            result.append(fix_info(correct.groupdict()))
        elif verbose:
            print('IGNORE: ', line)
    return result

def parseBond(file, verbose=False):

    nboSumLine = namedRe('NBO',r"\d+.",before='allow')
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
            result.append(fix_info(correct.groupdict()))
        elif verbose:
            print('IGNORE: ', line)
    return result

def parse3CBond(file, verbose=False):
    
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
            result.append(fix_info(correct.groupdict()))
        elif verbose:
            print('IGNORE: ', line)
    return result

def parseNboSum(file, verbose=False, options=["CR", "LP", "LV", "BD", "ABB", "3C", "3Cn", "3Cs"]):

    tables = {}
        
    def helper(df, list_1, list_2, list_3):
        df[list_1] = df[list_1].astype(int)
        df[list_2] = df[list_2].astype(str)
        df[list_3] = df[list_3].astype(float)
        return df
    
    for key in ['CR', 'LP', 'LV']:
        if key in options:
            pos = find(key, file)
            table = pd.DataFrame(parseNonBond(extract(file, pos), verbose=verbose))
            try:
                table = helper(table, ['NBO', 'AO','Loc'], ['Type','Atom'], ['Occupancy','Energy'])
                tables[key] = table
            except:
                if verbose:
                    print("Table of "+key+" is empty or in wrong format.")
    
    for key in ['BD', 'ABB']:
        if key in options:
            if key=='ABB':
                pos = find("BD*(", file)
            else:
                pos = find('BD ', file)
            tab = pd.DataFrame(parseBond(extract(file, pos), verbose=verbose))
            try:
                tab = helper(tab, ['NBO','AO','Loc1','Loc2'], ['Type','Atom1','Atom2'], ['Occupancy','Energy'])
                tables[key] = tab
            except:
                if verbose:
                    print("Table of "+key+" is empty or in wrong format.")

    for key in ['3C', '3Cn', '3Cs']:
        if key in options:
            if key=='3Cn':
                pos = find("3Cn", file)
            elif key=='3Cs':
                pos = find("3C*", file)
            else:
                pos = find("3C ", file)
            tab = pd.DataFrame(parse3CBond(extract(file, pos), verbose=verbose))
            try:
                tab = helper(tab, ['NBO','AO','Loc1','Loc2', 'Loc3'], ['Type','Atom1','Atom2','Atom3'], ['Occupancy','Energy'])
                tables[key] = tab
            except:
                if verbose:
                    print("Table of "+key+" is empty or in wrong format.")

    return tables