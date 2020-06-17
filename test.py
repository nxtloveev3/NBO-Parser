from rootNbo import nbo
from lib.parseNboSum import *
from lib.parseNBO import*
from lib.parseCMO import *
from lib.basicReadingFunctions import *
from lib.parsePert import*
import pandas as pd
import numpy as np

singlet = "CN01_NN13_S_nbo.log"
triplet = "CN01_NN01_T_nbo.log"

#table = nbo(singlet)
#cmo = table.cmo
#nboSum = table.nboSum
table = nbo(triplet)
#cmoAlpha = table.cmoA
#cmoBeta = table.cmoB
#nboAlpha = table.nboA
#nboBeta = table.nboB
nboPertA = table.pertAlpha

#tabCRF,tabLPF,tabLVF,tabBDF,tabBDSF,tab3CF,tab3CnF,tab3CsF = parseNboSum(nboSum)
#tabCRF,tabBDF,tabBDSF,tab3CF,tab3CnF,tab3CsF = parseNBO(nboSum)
#nboA = parseNBO(nboAlpha)
#tableNBOA = pd.DataFrame(nboA)

pertA = parsePert(nboPertA)
print(pertA)
