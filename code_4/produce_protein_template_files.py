#----------------------------------------------------------------------------------------
# Produce protein template files required for single protein structural modelling.
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
from id_mapping import produce_chainSeq_dict
from modelling_tools import (set_pdb_dir,
                             set_template_dir,
                             enable_pdb_downloads,
                             disable_pdb_warnings,
                             write_protein_template_sequences,
                             produce_fullmodel_chain_strucRes_dict,
                             produce_protein_template_files)

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # allow downloading of PDB structures while constructing the structural interactome
    allow_pdb_downloads = False
    
    # suppress PDB warnings when constructing the structural interactome
    suppress_pdb_warnings = True
    
    # parent directory of all data files
    dataDir = Path('../data')
        
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed template-related data files specific to interactome
    templateBasedDir = interactomeDir / 'template_based'
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # directory for PDB structure files
    pdbDir = Path('../../pdb_files')
    
    # directory for template structure files
    templateDir = modelBasedDir / 'protein_templates'
    
    # input data files
    templateMapFile = templateBasedDir / 'single_chain_map_per_protein.txt'
    chainSeqFile = templateBasedDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = templateBasedDir / 'protein_chain_strucRes.pkl'
    
    # output data files
    templateSeqFastaFile = modelBasedDir / 'protein_template_sequences.fasta'
    templateSeqFile = modelBasedDir / 'protein_template_sequences.pkl'
    templateStrucResFile = modelBasedDir / 'protein_template_strucRes.pkl'
    
    # create output directories if not existing
    if not modelBasedDir.exists():
        os.makedirs(modelBasedDir)
    if not pdbDir.exists():
        os.makedirs(pdbDir)
    if not templateDir.exists():
        os.makedirs(templateDir)
    
    # set directory of raw PDB coordinate files for modelling tools
    set_pdb_dir (pdbDir)
    
    # set directory of template coordinate files for modelling tools
    set_template_dir (templateDir)
    
    # enable or disable PDB downloads
    enable_pdb_downloads (allow_pdb_downloads)
    
    # suppress or allow PDB warnings
    disable_pdb_warnings (suppress_pdb_warnings)
    
    print('Extracting coordinate files for protein templates')
    produce_protein_template_files (templateMapFile, chainSeqFile, chainStrucResFile)
    
    print('writing protein template sequences to Fasta file')
    templateMap = pd.read_table(templateMapFile, sep='\t')
    templateIDs = []
    for template in templateMap["Subject"].values:
        pdbid, chainID = template.split('_')
        templateID = '-'.join([pdbid, chainID])
        templateIDs.append((templateID, chainID))
    write_protein_template_sequences (templateIDs,
                                      chainSeqFile,
                                      chainStrucResFile,
                                      templateSeqFastaFile)
    
    print('producing protein template sequence dictionary')
    produce_chainSeq_dict (templateSeqFastaFile, templateSeqFile)
    
    print('producing protein template structured residue label file')
    produce_fullmodel_chain_strucRes_dict (templateSeqFile, templateStrucResFile)

if __name__ == "__main__":
    main()
