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