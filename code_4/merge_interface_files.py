#----------------------------------------------------------------------------------------
# Merge multiple protein chain interface files into one file.
#----------------------------------------------------------------------------------------

import io
import numpy as np
from pathlib import Path
from structural_annotation import clear_interfaces, read_chain_interfaces, write_chain_interfaces

def main():
    
    # directory of input files
    inDir = '../data/processed/HuRI/model_based'
    
    # input file names
    inFiles = ['model_interfaces_%d.txt' % i for i in np.arange(1,11)]
    
    # path to output file
    outPath = Path(inDir) / 'model_interfaces.txt'
    
    with io.open(outPath, 'w') as fout:
        for i, inFile in enumerate(inFiles):
            inPath = Path(inDir) / inFile
            with io.open(inPath, "r", encoding="utf-8") as f:
                headers = f.readline()
                if i == 0:
                    fout.write(headers)
                for line in f:
                    fout.write(line)
    
    clear_interfaces()
    read_chain_interfaces (outPath)
    write_chain_interfaces (outPath)
    
if __name__ == "__main__":
    main()
