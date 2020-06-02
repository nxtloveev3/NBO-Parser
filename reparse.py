#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 15:56:54 2020

@author: yaron
"""
from lib import basicReadingFunctions as brf
import re

singlet = "CN01_NN13_S_nbo.log"
lines = brf.readlines(singlet)

startNAO = brf.find("NATURAL POPULATIONS:  Natural atomic orbital occupancies",lines)
endNAO =  endNAO = brf.find("Summary of Natural Population Analysis:",lines)

sec = lines[startNAO[0] : endNAO[0]]

#%%
def named_re(name, respec, before = 'none', after='require'):
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

# format of 999.999  (forces digits before and after the period)
re_float = r"-?\d+\.\d+"
re_atomic_element = r"[A-Z][a-z]?"

# Naming convention for dictionary is that at the top of this section in 
# the gaussian output file.
nao_line  = named_re('NAO'   , r"\d+", before = 'allow')  
nao_line += named_re('Atom' , re_atomic_element)      
nao_line += named_re('No'   , r"\d+")
nao_line += named_re('lang' , r"\w+")
nao_line += named_re('Type' , r'[A-Z][a-z]+', after='none')
nao_line += r"\("  + named_re('AO', r"\d+[spdfg]", 
                              before = 'allow', after = 'none')  + r"\)"
nao_line += named_re('Occupancy', re_float, before='require')
nao_line += named_re('Energy', re_float,after='allow')

nao_line_re = re.compile(nao_line)
res = []
for line1 in sec:
    t1 = nao_line_re.search(line1)
    if t1:
        print('PARSING: ',line1)
        res.append( t1.groupdict())
    else:
        print('IGNORE: ', line1)
#%%
#for t1 in res:
#    print(t1)