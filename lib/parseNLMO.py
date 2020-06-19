from .basicReadingFunctions import namedRe
import re
import string

def parseNLMO(file):
    tmpFile = []
    for line in file:
        newline = line
        if "BD*" in line:
            newline = line.replace("BD*","ABB") #ABB stand for antibonding bonds
        tmpFile.append(newline)
    
    ####Nlmo Regex####
    reFloat = r"-?\d+\.\d+"
    bondType = r"[A-Z][A-Z][A-Z]?"
    nonBondType = r"[A-Z][A-Z]"
    atom = r"[A-Z][a-z]?\s*\d+"
    nlmoNonBondLine = namedRe("NLMO", r"\d+",before='allow',after='none') + r"\." 
    nlmoNonBondLine += r"\s+\(" + namedRe("Occupancy", reFloat,before='allow',after='none') + r"\)"
    nlmoNonBondLine += namedRe("PercentFromNBO", reFloat,before='allow',after='none') + r"\%"
    nlmoNonBondLine += namedRe("Type", nonBondType,before='allow')
    nlmoNonBondLine += r"\(" + namedRe("Index", r"\d+",before="allow",after="none") + r"\)"
    nlmoNonBondLine += namedRe("Atom", atom,before='allow',after='allow')
    nlmoNonBondLineRe = re.compile(nlmoNonBondLine)

    nlmoBondLine = namedRe("NLMO", r"\d+",before='allow',after='none') + r"\."
    nlmoBondLine += r"\s+\(" + namedRe("Occupancy", reFloat,before='allow',after='none') + r"\)"
    nlmoBondLine += namedRe("PercentFromNBO", reFloat,before='allow',after='none') + r"\%"
    nlmoBondLine += namedRe("Type", bondType,before='allow')
    nlmoBondLine += r"\(" + namedRe("Index", r"\d+",before="allow",after="none") + r"\)"
    nlmoBondLine += namedRe("Atom", atom,before='allow',after='none') + r"\-"
    nlmoBondLine += namedRe("Atom2", atom,before='allow',after='allow')
    nlmoBondLineRe = re.compile(nlmoBondLine)

    nlmoHybridsLine = namedRe("Percent", reFloat,before='allow',after='none') + r"\%"
    nlmoHybridsLine += namedRe("Atom", atom,before='allow')
    nlmoHybridsLine += r"[s]\s*\(" + namedRe("sOrbit", reFloat,before='allow',after='none') + r"\%\)"
    nlmoHybridsLine += r"\s*" + r"([p]\s*\-?\d+\.\d+\s*\()?" + namedRe("pOrbit", r"(d+\.\d+)?",before='allow',after='none') + r"(\))?"
    nlmoHybridsLine += r"\s*" + r"([d]\s*\-?\d+\.\d+\s*\()?" + namedRe("dOrbit", r"(d+\.\d+)?",before='allow',after='none') + r"(\))?"
    nlmoHybridsLineRe = re.compile(nlmoHybridsLine)

    result = {}
    currentNLMO = None
    for line in tmpFile:
        if nlmoNonBondLineRe.search(line):
            info = nlmoNonBondLineRe.match(line).groupdict()
            currentNLMO = info["NLMO"]
            result[currentNLMO] = dict()
            result[currentNLMO]["Occupancy"] = float(info["Occupancy"])
            result[currentNLMO]["PercentFromNBO"] = float(info["PercentFromNBO"])
            result[currentNLMO]["Type"] = info["Type"]
            result[currentNLMO]["Index"] = int(info["Index"])
            result[currentNLMO]["Atom"] = info["Atom"]
            result[currentNLMO]["HybridContribution"] = []
        elif nlmoBondLineRe.search(line):
            info = nlmoBondLineRe.match(line).groupdict()
            currentNLMO = info["NLMO"]
            result[currentNLMO] = dict()
            result[currentNLMO]["Occupancy"] = float(info["Occupancy"])
            result[currentNLMO]["PercentFromNBO"] = float(info["PercentFromNBO"])
            result[currentNLMO]["Type"] = info["Type"]
            result[currentNLMO]["Index"] = int(info["Index"])
            result[currentNLMO]["Atom"] = info["Atom"]
            result[currentNLMO]["Atom2"] = info["Atom2"]
            result[currentNLMO]["HybridContribution"] = []
        elif nlmoHybridsLineRe.search(line):
            info = nlmoHybridsLineRe.match(line).groupdict()
            info["Percent"] = float(info["Percent"])
            info["sOrbit"] = float(info["sOrbit"])
            if not info["pOrbit"].isspace:
                info["pOrbit"] = float(info["pOrbit"])
            if not info["dOrbit"].isspace:
                info["dOrbit"] = float(["dOrbit"])
            result[currentNLMO]["HybridContribution"].append(info)  
        else:
            print("Ignoring line", line)

    return result