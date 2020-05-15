def readlines(file):
    f = open(file,'rb')
    lines = f.readlines()
    lines = [line.decode("utf-8").strip() for line in lines]
    return lines

class nbo(object):
    def __init__(self,file):
        lines = readlines(file)
        nboSumstart = lines.index("NATURAL BOND ORBITALS (Summary):")
        startNAO = lines.index("NATURAL POPULATIONS:  Natural atomic orbital occupancies")
        endNAO = lines.index("Summary of Natural Population Analysis:")
        startNLMO = lines.index("NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:")
        self.nboAll = lines[startNAO:endNAO]
        print(len(self.nboAll))


nbo("CN01_NN01_opt_opt.log_r_nbo.log")