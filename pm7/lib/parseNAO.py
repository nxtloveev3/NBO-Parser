import re
from .basicReadingFunctions import namedRe

reFloat = r"-?\d+\.\d+"
energy = r"(-?\d+\.\d+)?"
lang = r"[A-Z]\d?[A-z]?[a-z]?"
atom = r"[A-Z][a-z]?"
types = r"[A-Z][a-z][a-z]"
naoLine = namedRe("NAO", r"\d+")
naoLine += namedRe("NO", r"\d+")
naoLine += namedRe("Atom", atom, after="allow")
naoLine += namedRe("lang", lang, after="allow")
naoLine += namedRe("Type", types, after="allow")
naoLine += r"\(" + namedRe("AO", r"\d+[a-z]", before="allow", after="none") + r"\)\s+"
naoLine += namedRe("Occupancy", reFloat, after="allowMore")
naoLine += namedRe("Energy",energy, after="allow")
naoLineRe = re.compile(naoLine)

def parseNAO(file, verbose=False):
    result = []
    for line in file:
        newLine = {}
        if naoLineRe.match(line):
            info = naoLineRe.match(line).groupdict()
            newLine["NAO"] = int(info["NAO"])
            newLine["NO"] = int(info["NO"])
            newLine["Atom"] = info["Atom"]
            newLine["lang"] = info["lang"]
            newLine["Type"] = info["Type"]
            newLine["AO"] = info["AO"]
            newLine["Occupancy"] = float(info["Occupancy"])
            if info["Energy"]:
                newLine["Energy"] = float(info["Energy"])
            result.append(newLine)
        else:
            if verbose:
                print("Ignoring this line:", line)
    return result