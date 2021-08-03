import pandas as pd
import openpyxl

def discriptorParserNao(file, option = None, verbose = False):
    fd = openpyxl.load_workbook(file)
    sheet = fd['Sheet1']

    if option == 'naoS':
        nao_S = [ 'Valence s occupancy Singlet',
            'Valence s energy Singlet',
            'Valence px occupancy Singlet',
            'Valence px energy Singlet',
            'Valence py occupancy Singlet',
            'Valence py energy Singlet',
            'Valence pz occupancy Singlet',
            'Valence pz energy Singlet']
        for num in range(len(nao_S)):
            nao_S[num] = []
    elif option == 'naoT':
        nao_T = [ 'Valence s occupancy Triplet',
            'Valence s energy Triplet',
            'Valence px occupancy Triplet',
            'Valence px energy Triplet',
            'Valence py occupancy Triplet',
            'Valence py energy Triplet',
            'Valence pz occupancy Triplet',
            'Valence pz energy Triplet',
            ]
        for num in range(len(nao_T)):
            nao_T[num] = []
    else: print("Require option naoS, or naoT")

    for row in range(2, sheet.max_row + 1):
        type = sheet['F' + str(row)].value
        if type == 'Val':
            lang = sheet['E' + str(row)].value
            occupancy = sheet['H' + str(row)].value
            energy = sheet['I' + str(row)].value
            if lang == 'S' and option == 'naoS': 
                nao_S[0].append(occupancy)
                nao_S[1].append(energy)
            elif lang == 'Px' and option == 'naoS': 
                nao_S[2].append(occupancy)
                nao_S[3].append(energy)
            elif lang == 'Py' and option == 'naoS': 
                nao_S[4].append(occupancy)
                nao_S[5].append(energy) 
            elif lang == 'Pz' and option == 'naoS': 
                nao_S[6].append(occupancy)
                nao_S[7].append(energy)
            elif lang == 'S' and option == 'naoT': 
                nao_T[0].append(occupancy)
                nao_T[1].append(energy)
            elif lang == 'Px' and option == 'naoT': 
                nao_T[2].append(occupancy)
                nao_T[3].append(energy)
            elif lang == 'Py' and option == 'naoT': 
                nao_T[4].append(occupancy)
                nao_T[5].append(energy) 
            elif lang == 'Pz' and option == 'naoT': 
                nao_T[6].append(occupancy)
                nao_T[7].append(energy)
    if option == 'naoS':
        result = pd.DataFrame(nao_S)
        result.fillna(0)
        if verbose: print(result)
        return result
    else:
        result = pd.DataFrame(nao_T)
        result.fillna(0)
        if verbose: print(result)
        return result