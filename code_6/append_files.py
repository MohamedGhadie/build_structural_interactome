#----------------------------------------------------------------------------------------
# Append multiple text files together.
#----------------------------------------------------------------------------------------

import io
import numpy as np
from pathlib import Path

def main():
    
    # directory of input files
    inDir = '../data/processed/HuRI/model_based'
    
    # input file names
    inFiles = ['structural_interactome_withDuplicates_%d.txt' % i for i in np.arange(1,11)]
    
    # path to output file
    outPath = Path(inDir) / 'structural_interactome_withDuplicates.txt'
    
    with io.open(outPath, 'w') as fout:
        for i, inFile in enumerate(inFiles):
            inPath = Path(inDir) / inFile
            with io.open(inPath, "r", encoding="utf-8") as f:
                headers = f.readline()
                if i == 0:
                    fout.write(headers)
                for line in f:
                    fout.write(line)

if __name__ == "__main__":
    main()