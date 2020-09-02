import re
from .basicReadingFunctions import namedRe, find, findExact, extractTab, replacingDigit

def parseNAO(file, verbose=True):
    reFloat = r"-?\d+\.\d+"
    lang = r"[a-z][a-z]?\d*[a-z]?\d*"
    atom = r"[A-Z][a-z]?"
    types = r"[A-Z][a-z][a-z]"
    naoLine = namedRe("NAO", r"\d+")
    naoLine += namedRe("Atom", atom, after="allow")
    naoLine += namedRe("NO", r"\d+")
    naoLine += namedRe("lang", lang, after="allow")
    naoLine += namedRe("Type", types, after="allow")
    naoLine += r"\(" + namedRe("AO", r"\d+[a-z]", before="allow", after="none") + r"\)\s+"
    naoLine += namedRe("Occupancy", reFloat)
    naoLine += namedRe("Energy",reFloat, after="allow")
    naoLineRe = re.compile(naoLine)
    
    result = []
    for line in file:
        newLine = {}
        if naoLineRe.match(line):
            info = naoLineRe.match(line).groupdict()
            newLine["NAO"] = int(info["NAO"])
            newLine["Atom"] = info["Atom"]
            newLine["NO"] = int(info["NO"])
            newLine["lang"] = info["lang"]
            newLine["Type"] = info["Type"]
            newLine["AO"] = info["AO"]
            newLine["Occupancy"] = float(info["Occupancy"])
            newLine["Energy"] = float(info["Energy"])
            result.append(newLine)
        else:
            if verbose:
                print("Ignoring this line:", line)
    return result