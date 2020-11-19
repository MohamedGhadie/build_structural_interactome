#----------------------------------------------------------------------------------------
# Split files in one folder into multiple folders.
#----------------------------------------------------------------------------------------

import sys
from os import listdir
from subprocess import run

def main():
    
    # number of files wanted in each split folder
    numFiles = 111
    
    # path to folder to split
    inDir = './jobs'
    
    infiles = [f for f in listdir(inDir) if not f.startswith('.')]
    numDir, n = 0, len(infiles)
    
    for i, f in enumerate(infiles):
        sys.stdout.write('  File %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        if i % numFiles == 0:
            numDir += 1
            outDir = str(inDir) + '_split_' + str(numDir) + '/'
            run(['mkdir', outDir])
        run(['cp', '/'.join([inDir, f]), outDir])
    print()

if __name__ == "__main__":
    main()