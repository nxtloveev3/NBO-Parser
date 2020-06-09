import re
from .basicReadingFunctions import find, findExact, extractTab
import pandas as pd

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