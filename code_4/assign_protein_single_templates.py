#----------------------------------------------------------------------------------------
# Select one template chain per protein from the protein-chain alignment table.
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
from structural_annotation import single_chain_per_protein
from pdb_tools import download_structures

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # download missing template structures from PDB
    download_template_structures = True
    
    # parent directory of all data files
    #dataDir = Path('/Volumes/MG_Samsung/edgotype_fitness_effect_full_model/data')
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
        
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    templateBasedDir = interactomeDir / 'template_based'
    
    # directory of PDB structure files
    pdbDir = Path('../../pdb_files')
    
    # input data files
    chainMapFile = templateBasedDir / 'struc_interactome_chain_map.txt'
    chainSeqFile = templateBasedDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = templateBasedDir / 'protein_chain_strucRes.pkl'
    
    # output data files
    singleChainMapFile = templateBasedDir / 'single_chain_map_per_protein.txt'
    chainIDFile = templateBasedDir / 'single_model_chainIDs.txt'
    pdbIDFile = templateBasedDir / 'single_model_complexIDs.txt'
    
    # create output directories if not existing
    if not templateBasedDir.exists():
        os.makedirs(templateBasedDir)
    if not pdbDir.exists():
        os.makedirs(pdbDir)

    #------------------------------------------------------------------------------------
    # select one chain model per protein
    #------------------------------------------------------------------------------------
    
    print('selecting one chain model per protein')
    single_chain_per_protein (chainMapFile,
                              singleChainMapFile,
                              chainSeqFile = chainSeqFile,
                              chainStrucResFile = chainStrucResFile,
                              pdbDir = pdbDir)
    
    proteinModels = pd.read_table (singleChainMapFile, sep='\t')
    uniqueChains = set(proteinModels["Subject"].values)
    uniquePDBs = {id.split('_')[0] for id in uniqueChains}
    print('%d unique chains in %d unique PDB structures' % (len(uniqueChains), len(uniquePDBs)))
    
    with open(chainIDFile, 'w') as f:
        for i in sorted(uniqueChains):
            f.write("%s\n" % i)
    with open(pdbIDFile, 'w') as f:
        for i in sorted(uniquePDBs):
            f.write("%s\n" % i)
    
    if download_template_structures:
        print('downloading structures for selected protein models')
        download_structures (pdbIDFile, pdbDir)

if __name__ == "__main__":
    main()
