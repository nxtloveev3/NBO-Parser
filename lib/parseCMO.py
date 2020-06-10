import re
from .basicReadingFunctions import find, findExact, extractTab
import pandas as pd

'''
Data structure for:
    MO   1 (occ): orbital energy =  -14.43595 a.u.
                0.707*[ 28]: CR ( 1) N41(cr)
                0.707*[ 40]: CR ( 1) N61(cr)
    MO  41 (occ): orbital energy =   -1.10114 a.u.
                0.426*[ 80]: BD ( 1) C22- N41
                0.426*[106]: BD ( 1) C42- N61
                0.413*[ 85]: BD ( 1) C23- N41
                0.413*[111]: BD ( 1) C43- N61 

    MO 129 (occ): orbital energy =   -0.31378 a.u.
                0.460*[ 56]: BD ( 2) C 3- C 4
                0.460*[ 69]: BD ( 2) C 8- C 9
               -0.387*[ 67]: BD ( 2) C 7- N20
               -0.387*[ 54]: BD ( 2) C 2- N21
               -0.257*[160]: BD*( 2) C10- C11*
               -0.257*[147]: BD*( 2) C 5- C 6*    
  129 = MO number
  -0.31378 = energy
   0.460 = coef   56 = nbo_num


res = dict()
res[MO_number] = {'energy': -0.31378
                  'wf' : [(coef, nbo_num), ...]}
res[1] = {'energy': -14.43595
          'wf' : [(0.707, 28), (0.707, 40)]}

re_energy_line
re_wf_line
res = {}
current_MO = None
unparsed = []
for line in line_range:
    if re_energy_line.match(line):
        # extract current_MO and energy
        current_MO =
        energy = 
       res[current_MO] = dict()
       res[current_MO]['energy'] = energy
       res[current_MO]['wf'] = []
    elif re_wf_line.match(line):
        nbo = 
        coef = 
       res[current_MO]['wf'].append( (nbo, coef) )
    else:
        unparsed.append(line)
        print("parseCMON: WARNING line not recognized "+ line)
        
if len(unparsed) > 0:
    raise Exception('parseCMON: unrecognized '+ str(len(unparsed))+ ' lines')

could add basic validation (later):
    MO_numbers are continuous from 1..N
    could parse the BD C 3- C 4 and check for agreement with info from
    parsing NBO table

'''


def parseCMON(file):
    loc = find("Molecular Orbital Atom-Atom Bonding Character",file)
    cmon = file[:loc[0]]
    tmpFile = []
    for line in cmon:
        newline = line
        if "BD*" in line:
            newline = line.replace("BD*","ABB") #ABB stand for antibonding bonds
        tmpFile.append(newline)
    posMO = find("MO",tmpFile)
    try:
        newFile = parseNonBond(tmpFile)
    except:
        newFile = parseBond(tmpFile)
    
    print(newFile)

def parseNonBond(file):
    reFloat = r"-?\d+\.\d+"
    atomicType = r"\d*[A-Z][A-Z]?[a-z]?"
    atom = r"[A-Z][a-z]?"
    cmonLine = r"\s*" + reFloat
    cmonLine += r"\*\[" + r"\s*\d+" + r"\]\:"
    cmonLine += r"\s*" + atomicType
    cmonLine += r"\s*" + r"\(" + r"\s*\d+" + r"\)"
    cmonLine += r"\s*" + atom
    cmonLine += r"\s*" + r"\d+"
    cmonLineRe = re.compile(cmonLine)

    for line in file:
        correct = cmonLineRe.search(line)
        if correct:
            line = correct.group()
        else:
            print('IGNORE: ', line)
    return file

def parseBond(file):
    reFloat = r"-?\d+\.\d+"
    atomicType = r"\d*[A-Z][A-Z][A-Z]?"
    atom = r"[A-Z][a-z]?"
    cmonLine = r"\s*" + reFloat
    cmonLine += r"\*\[" + r"\s*\d+" + r"\]\:"
    cmonLine += r"\s*" + atomicType
    cmonLine += r"\s*" + r"\(" + r"\s*\d+" + r"\)"
    cmonLine += r"\s*" + atom
    cmonLine += r"\s*" + r"\d+"
    cmonLine += r"\-" + r"\s*" + atom + r"\s*\d+" + r"\*"
    cmonLineRe = re.compile(cmonLine)

    for line in file:
        correct = cmonLineRe.search(line)
        if correct:
            line = correct.group()
        else:
            print('IGNORE: ', line)
    return file