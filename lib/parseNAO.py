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
            nao = info["NAO"]
            newLine["NAO"] = int(nao)
            atom = info["Atom"]
            newLine["Atom"] = atom
            number = info["NO"]
            newLine["NO"] = int(number)
            lang = info["lang"]
            newLine["lang"] = lang
            types = info["Type"]
            newLine["Type"] = types
            ao = info["AO"]
            newLine["AO"] = ao
            occupancy = info["Occupancy"]
            newLine["Occupancy"] = float(occupancy)
            energy = info["Energy"]
            newLine["Energy"] = float(energy)
            result.append(newLine)
        else:
            if verbose:
                print("Ignoring this line:", line)
    return result