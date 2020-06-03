import re
from .basicReadingFunctions import namedRe

def parseNboSum(file):
    reFloat = r"-?\d+\.\d+"
    atomicType = r"\d*[A-Z][A-Z]?[a-z]?\*?"
    atom = r"[A-Z][a-z]?\s*\d+\-?\s*[A-Z]?[a-z]?\s*\d+?"
    nboSumLine = namedRe('NBO',r"\d+\.",after='allow')
    nboSumLine += namedRe('Type',atomicType,after='allow')
    nboSumLine += r"\(" + namedRe('AO',r"\d+",before='allow',after='none') + r"\)"
    nboSumLine += namedRe('Atom',atom)
    nboSumLine += namedRe('Occupancy', reFloat,after='require')
    nboSumLine += namedRe('Occupancy', reFloat,after='allow')
    nboLineRe = re.compile(nboSumLine)

    for line in file:
        correct = naoLineRe.search(line)
        if correct:
            print('PARSING: ',line)
            result.append(correct.groupdict())
        else:
            print('IGNORE: ', line)
    return result