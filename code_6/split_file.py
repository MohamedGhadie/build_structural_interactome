#----------------------------------------------------------------------------------------
# Split a text file into multiple files.
#----------------------------------------------------------------------------------------

import io
import numpy as np
from pathlib import Path

def main():
    
    # number of lines wanted in each split file
    outLen = 193
    
    # directory of input file
    inDir = '../data/processed/HuRI/model_based_4'
    
    # input file name
    inFile = 'model_annotated_interactome.txt'
    
    # path to input file
    inPath = Path(inDir) / inFile
    
    # count number of lines to split, not including header line
    with io.open(inPath, "r", encoding="utf-8") as f:
        for inLen, _ in enumerate(f):
            pass
    
    numSplits = int(np.ceil(inLen / outLen))
    filename = inFile.split('.')[0] if inFile.endswith('.txt') else inFile
    outPaths = [Path(inDir) / (filename + '_' + str(i) + '.txt') for i in np.arange(1, numSplits + 1)]
    if outPaths:
        with io.open(inPath, "r", encoding="utf-8") as f:
            headers = f.readline()
            fout = open(outPaths[0], "w")
            fout.write(headers)
            for i, line in enumerate(f):
                ind = int(i/outLen)
                if (ind > 0) and (i % outLen == 0):
                    fout.close()
                    fout = open(outPaths[ind], "w")
                    fout.write(headers)
                fout.write(line)
            fout.close()

if __name__ == "__main__":
    main()
