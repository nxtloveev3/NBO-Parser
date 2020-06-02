from rootNbo import nbo
from lib.parseNBO import *
from lib.parseCMO import *

singlet = "CN01_NN13_S_nbo.log"
triplet = "CN01_NN01_T_nbo.log"

table = nbo(singlet)
cmo = table.cmo
nbo = table.nbo

#table = nbo(triplet)
#cmoAlpha = table.cmoA
#cmoBeta = table.cmoB
#nboAlpha = table.nboA
#nboBeta = table.nboB

cmoTab = parseCMO(cmo)
#nboTab = parseNBO(nbo)
print(cmoTab)


