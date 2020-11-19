#----------------------------------------------------------------------------------------
# Produce protein full model mapping data required for further structural calculations.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from text_tools import produce_item_list
from id_mapping import produce_chain_dict
from modelling_tools import (produce_protein_fullmodel_chainSeq_dict,
                             produce_fullmodel_chain_strucRes_dict,
                             produce_protein_fullmodel_pos_mapping)

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # input data files
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    templateMapFile = modelBasedDir / 'single_template_map_per_protein.txt'
    
    # output data files
    chainSeqFile = modelBasedDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = modelBasedDir / 'protein_chain_strucRes.pkl'
    chainMapFile = modelBasedDir / 'single_chain_map_per_protein.txt'
    chainListFile = modelBasedDir / 'protein_model_chains.list'
    modelChainsFile = modelBasedDir / 'protein_model_chains.pkl'
    
    if not modelBasedDir.exists():
        os.makedirs(modelBasedDir)
    
    print('producing protein model chain sequence dictionary')
    produce_protein_fullmodel_chainSeq_dict (templateMapFile, proteinSeqFile, chainSeqFile)
    
    print('producing protein model chain structured residue label file')
    produce_fullmodel_chain_strucRes_dict (chainSeqFile, chainStrucResFile)
    
    print('producing protein model chain position mapping file')
    produce_protein_fullmodel_pos_mapping (templateMapFile, chainSeqFile, chainMapFile)
    
    print('producing model chain ID list')
    produce_item_list (chainMapFile, "Subject", chainListFile)
    
    print('producing model chain dictionary from chain list file')
    produce_chain_dict (chainListFile, modelChainsFile)

if __name__ == "__main__":
    main()
