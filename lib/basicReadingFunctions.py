def readlines(file):
    f = open(file,'rb')
    lines = f.readlines()
    lines = [line.decode("utf-8").strip() for line in lines]
    return lines

# This function returns the index of specific line based on key phrase. It's primarily used 
# for finding different section of the file     
def find(text,file):
    result = []
    count = 0
    for elem in file:
        if text in elem:
            result.append(count)
        count += 1
    return result

# Similar to the last function but locates the exact text.
def findExact(text,file):
    result = []
    count = 0
    for elem in file:
        if text == elem:
            result.append(count)
        count += 1
    return result

def replacingDigit(index,file):
    file = file.split()
    if file[index].endswith("."): #Treat the number so it is easier to recognize as digit (Ex: 120. -> 120)
        file[index] = file[index][:-1]
    return file

#Following two equation takes in information and generates new table that takes unnecessary information out.
def fixEnding(posFile,tmpFile):
    newFile = []
    for num in posFile:
        line = tmpFile[num].split()
        newline = shortingLine(line)
        newFile.append(newline)
    return newFile

def shortingLine(file):
    count = 0
    for elem in file[-1]:
        if elem.isalpha():
            count += 1
    if count == 0:
        return file
    else:
        result = shortingLine(file[:-1])
        return result
    