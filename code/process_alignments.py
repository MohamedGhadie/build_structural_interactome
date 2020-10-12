import os
from pathlib import Path
from text_tools import parse_blast_file
from structural_annotation import filter_chain_annotations
from alignment_tools import extend_alignments

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # Maximum e-value cutoff to filter out protein-chain annotations
    evalue = 1e-10
    
    # Minimum protein coverage fraction required for protein-chain annotation
    proteinCov = 0
    
    # Minimum chain coverage fraction required for protein-chain annotation
    chainCov = 0
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    blastFile = interactomeDir / 'interactome_template_blast_alignments_e-5'
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    chainStrucSeqFile = interactomeDir / 'chain_struc_sequences.pkl'
    
    # output data files
    alignmentFile1 = interactomeDir / 'interactome_template_alignments.txt'
    alignmentFile2 = interactomeDir / 'interactome_template_filtered_alignments.txt'
    alignmentFile3 = interactomeDir / 'interactome_template_extended_alignments.txt'
    
    # create output directories if not existing
    if not dataDir.exists():
        os.makedirs(str(dataDir))
    if not procDir.exists():
        os.makedirs(str(procDir))
    if not interactomeDir.exists():
        os.makedirs(str(interactomeDir))
    
    if not alignmentFile1.is_file():
        print( 'Parsing BLAST protein-chain alignment file' )
        parse_blast_file (blastFile, alignmentFile1)
    
    if not alignmentFile2.is_file():
        print('Filtering alignments')
        filter_chain_annotations (alignmentFile1,
                                  alignmentFile2,
                                  evalue = evalue,
                                  prCov = proteinCov,
                                  chCov = chainCov)
    
    if not alignmentFile3.is_file():
        print('Extending alignments to full sequences')
        extend_alignments (alignmentFile2,
                           proteinSeqFile,
                           chainStrucSeqFile,
                           alignmentFile3)
    
if __name__ == "__main__":
    main()
