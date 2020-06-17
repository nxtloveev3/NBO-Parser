from .basicReadingFunctions import namedRe, find
import re

def parsePert(file):
    tmpFile = []
    for line in file:
        newline = line
        if "BD*" in line:
            newline = line.replace("BD*","ABB") #ABB stand for antibonding bonds
        tmpFile.append(newline)
    
    ####pert Regex ####
    reFloat = reFloat = r"-?\d+\.\d+"
    atomicType = r"[A-Z][A-Z][A-Z]?"
    atom = r"[A-Z][a-z]?\s*\d+(\-\s*[A-Z][a-z]?\s*\d+)?"
    pertLine = namedRe("NBO1", r"\d+",before='allow',after='none') + r"\."
    pertLine += namedRe("Type1", atomicType,before='allow',after='allow')
    pertLine += r"\(" + namedRe("Index1", r"\d+",before='allow',after='none') + r"\)"
    pertLine += namedRe("AtomInfo1", atom,before='allow')
    pertLine += namedRe("NBO2", r"\d+",after='none') + r"\."
    pertLine += namedRe("Type2", atomicType,before='allow',after='allow')
    pertLine += r"\(" + namedRe("Index2", r"\d+",before='allow',after='none') + r"\)"
    pertLine += namedRe("AtomInfo2", atom,before='allow')
    pertLine += namedRe("E", reFloat)
    pertLine += namedRe("EE", reFloat)
    pertLine += namedRe("F", reFloat,after='allow')
    pertLineRe = re.compile(pertLine)

    result={}
    number = 0
    for line in tmpFile:
        if pertLineRe.search(line):
            number += 1
            count = str(number)
            result[count] = dict()
            info = pertLineRe.search(line).groupdict()
            donnor = dict()
            nbo1 = info["NBO1"]
            atomicType1 = info["Type1"]
            index1 = info["Index1"]
            atom1 = info["AtomInfo1"]
            donnor["NBO"] = int(nbo1)
            donnor["Type"] = atomicType1
            donnor["Index"] = int(index1)
            donnor["Atom"] = atom1
            acceptor = dict()
            nbo2 = info["NBO2"]
            atomicType2 = info["Type2"]
            index2 = info["Index2"]
            atom2 = info["AtomInfo2"]
            acceptor["NBO"] = int(nbo2)
            acceptor["Type"] = atomicType2
            acceptor["Index"] = int(index2)
            acceptor["Atom"] = atom2
            e = info["E"]
            ee = info["EE"]
            f = info["F"]
            result[count]["DonorNBO"] = donnor
            result[count]["AcceptorNBO"] = acceptor
            result[count]["E"] = float(e)
            result[count]["EE"] = float(ee)
            result[count]["F"] = float(f)
        else:
            print("Ignoring line:", line)
    return result


