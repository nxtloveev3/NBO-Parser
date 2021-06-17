from .basicReadingFunctions import namedRe, find, fix_badatom
import re

reFloat = r"-?\d+\.\d+"
atom = r"\d+\s?[A-Z][a-z]?"
types = r"[A-Z]+\*?"

nboSumLine = namedRe("NO", r"\d+", after = "allow") + r"\:\s+"
nboSumLine += namedRe("Occupancy", reFloat)
nboSumLine += namedRe("Energy",reFloat)
nboSumLine += namedRe("Type", types, after="allow")
nboSumLine += r"\(" + namedRe("Bondorbital", r"\d+", before ="allow", after="allow") + r"\)"
nboSumLine += r"\[" + namedRe("Atom1", atom, before = "allow", after = "allow") + r"\]"
nboSumLine += r"\[" + namedRe("Atom2", atom, before = "allow", after = "allow") + r"\]"
nboSumRe = re.compile(nboSumLine)


def parseNboSum(file, verbose=False):
    result = []
    for line in file:
        newLine = {}
        if nboSumRe.match(line):
            info = nboSumRe.match(line).groupdict()
            newLine["NO"] = int(info["NO"])
            newLine["Atom1"] = info["Atom1"]
            newLine["Atom2"] = info["Atom2"]
            newLine["Type"] = info["Type"]
            newLine["Occupancy"] = float(info["Occupancy"])
            newLine["Energy"] = float(info["Energy"])
            result.append(newLine)
        else:
            if verbose:
                print("Ignoring this line:", line)
    return result