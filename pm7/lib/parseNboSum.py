from .basicReadingFunctions import namedRe, find, fix_badatom
import re

reFloat = r"-?\d+\.\d+"
atom = r"([A-Z][a-z]?)?"
types = r"[A-Z]+\*?"
pLoc = r"(\d+)?"

nboSumLine = namedRe("NO", r"\d+", after = "allow") + r"\:\s+"
nboSumLine += namedRe("Occupancy", reFloat)
nboSumLine += namedRe("Energy",reFloat)
nboSumLine += namedRe("Type", types, after="allow")
nboSumLine += r"\(" + namedRe("Bondorbital", r"\d+", before ="allow", after="allow") + r"\)\s?"
nboSumLine += r"\[" + namedRe("Loc1", r"\d+", before = "allow", after="allow")
nboSumLine += namedRe("Atom1", atom, after = "allow") + r"\]\s?"
nboSumLine += r"\[?" + namedRe("Loc2", pLoc, before = "allow", after="allow")
nboSumLine += namedRe("Atom2", atom, after = "allow") + r"\]?"
nboSumRe = re.compile(nboSumLine)


def parseNboSum(file, verbose=False):
    result = []
    for line in file:
        newLine = {}
        if nboSumRe.match(line):
            info = nboSumRe.match(line).groupdict()
            newLine["NO"] = int(info["NO"])
            newLine["Atom1"] = info["Atom1"]
            newLine["Loc1"] = info["Loc1"]
            newLine["Atom2"] = info["Atom2"]
            newLine["Loc2"] = info["Loc2"]
            newLine["Type"] = info["Type"]
            newLine["Occupancy"] = float(info["Occupancy"])
            newLine["Energy"] = float(info["Energy"])
            result.append(newLine)
        else:
            if verbose:
                print("Ignoring this line:", line)
    return result