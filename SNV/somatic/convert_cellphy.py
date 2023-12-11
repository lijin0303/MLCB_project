#Convert a mutation matrix to CellPhy PHYLIP format

import sys

heteroDict = {'AC':'M', 'CA':'M', 'AG':'R', 'GA':'R', 'AT':'W', 'TA':'W', 'CG':'S', 'GC':'S', 'CT':'Y', 'TC':'Y', 'GT':'K', 'TG':'K'}

inFile = open(sys.argv[1], "r")
lines = inFile.readlines()
nCells = len(lines) - 1
header = lines[0].strip().split()
snps = ["".join(snp.split(":")[2:]) for snp in header]
nSnps = len(snps)
print(nCells, nSnps)
for line in lines[1:]:
    tokens = line.strip().split()
    print(tokens[0], end="\t")
    for i in range(1, len(tokens)):
        gt = tokens[i]
        match gt:
            case '1/0':
                print(snps[i-1][0], end="")
            case '0/1':
                print(snps[i-1][1], end="")
            case '1/1':
                print(heteroDict[snps[i-1]], end="")
            case _:
                print('N', end="")
    print()

