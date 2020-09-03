import string
import pandas as pd
import numpy as np
import re

def readlines(file):
    f = open(file,'rb')
    lines = f.readlines()
    lines = [line.decode("utf-8").strip() for line in lines]
    return lines

# This function returns the index of specific line based on key phrase. It's primarily used 
# for finding different section of the file     
def find(text,file):
    result = []
    count = 0
    for line in file:
        if text in line:
            result.append(count)
        count += 1
    return result

def extract(file, pos):
    return [file[i] for i in pos]

#Following two equation takes in information and generates new table that takes unnecessary information out.
def fixEnding(posFile,tmpFile):
    newFile = []
    for num in posFile:
        line = tmpFile[num].split()
        newline = shortingLine(line)
        newFile.append(newline)
    return newFile

def shortingLine(file):
    count = 0
    for elem in file[-1]:
        if elem.isalpha():
            count += 1
    if count == 0:
        return file
    else:
        result = shortingLine(file[:-1])
        return result

def finalClean(file): # This part takes out all the nonessential component in the data table
    result = []
    for line in file:
        newstring = ""
        for elem in line:
            newstring += elem + " "
        newstring = newstring.replace(")", " ")
        newstring = newstring.replace("(", " ")
        newstring = newstring.replace("?", " ")
        newstring = newstring.replace("-", " ", 1)
        newstring = newstring.replace("*[", " ")
        newstring = newstring.replace("]", " ")
        newstring = newstring.replace("=", " ")
        newstring = newstring.replace(":", " ")
        newLine = newstring.split()
        result.append(newLine)
    return result

def extractTab(file,pos):
    result = []
    for i in pos:
        line = ""
        if isinstance(file[i],list):
            for elem in file[i]:
                line += elem + " "
        else:
            for elem in file[i].split():
                line += elem + " "
        result.append(line)
    return result

def toTable(file,titles):
    columnTitle = titles.split(",")
    totalCol = len(columnTitle) 
    dataInColForm = [[]for i in range(totalCol)]
    for line in file:
        for i in range(totalCol):
            dataInColForm[i].append(line[i])
    dataTable = pd.DataFrame(np.array(dataInColForm),columns=columnTitle)
    return dataTable


def namedRe(name, respec, before='none', after='require'): #function written by David Yaron
    '''
      wraps a regular expression in a block that gives it a name:
          (?P<name> respec)
      before and after are whitespace requirements before and after the match
         'none'   : do not allow any whitespace
         'allow' : accept but do not require whitespace
         'require' : require whitespace
    '''
    ws = {'none'    : r'',
          'allow'   : r'\s*',
          'require' : r'\s+'}
    res = ws[before] + "(?P<" + name + ">" + respec + ")" + ws[after]
    return res

def fix_badatom(atom):
    if ' ' in atom:
        return atom.split()
    a = ''
    i = ''
    for chr_ in atom:
        if chr_.isalpha():
            a += chr_
        else:
            i += chr_
    return [a, i]