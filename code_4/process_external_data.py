#----------------------------------------------------------------------------------------
# Process raw data files from external sources.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from text_tools import parse_fasta_file, reduce_fasta_headers
from id_mapping import (produce_geneName_dict,
                        produce_uniqueGene_swissProtIDs,
                        produce_uniqueGene_sequences,
                        produce_proteinSeq_dict,
                        produce_uniprotID_dict,
                        produce_rnaToProtein_refseqID_dict)

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # input data files
    uniprotRefSeqFile = extDir / 'UP000005640_9606.fasta'
    idmapFile = extDir / 'HUMAN_9606_idmapping.dat'
    proteomeListFile = extDir / 'uniprot_reviewed_human_proteome.list'
    refseqIDFile = extDir / 'LRG_RefSeqGene'
    pdbSeqresFile = extDir / 'pdb_seqres.txt'
    
    # output data files
    refSeqFile = procDir / 'human_reference_sequences.fasta'
    sequenceFile = procDir / 'human_reference_sequences.txt'
    geneMapFile = procDir / 'to_human_geneName_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    uniqueGeneSequenceFile = procDir / 'human_unique_gene_reference_sequences.txt'
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    rnaToProteinRefseqIDMapFile = procDir / 'human_rnaToProtein_refseqID_map.pkl'
    seqresFile = procDir / 'pdb_seqres_reduced.fasta'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    
    if not refSeqFile.is_file():
        print('Reducing headers in protein sequence fasta file')
        reduce_fasta_headers (uniprotRefSeqFile, '|', 2, 2, refSeqFile)
    
    if not sequenceFile.is_file():
        print('reading protein sequence fasta file')
        parse_fasta_file (refSeqFile, sequenceFile)
    
    if not geneMapFile.is_file():
        print('producing UniProtID-to-geneName dictionary')
        produce_geneName_dict (idmapFile, proteomeListFile, geneMapFile)
    
    if not uniqueGeneSwissProtIDFile.is_file():
        print('producing list of unique-gene UniProt IDs')
        produce_uniqueGene_swissProtIDs (proteomeListFile, geneMapFile, uniqueGeneSwissProtIDFile)
    
    if not uniqueGeneSequenceFile.is_file():
        print('producing sequence file for unique-gene UniProt IDs')
        produce_uniqueGene_sequences (sequenceFile, uniqueGeneSwissProtIDFile, geneMapFile, uniqueGeneSequenceFile)
    
    if not proteinSeqFile.is_file():
        print('producing protein sequence dictionary')
        produce_proteinSeq_dict (refSeqFile, proteinSeqFile)
    
    if not uniprotIDmapFile.is_file():
        print('producing to-UniProt-ID dictionary')
        produce_uniprotID_dict (idmapFile, uniqueGeneSwissProtIDFile, uniprotIDmapFile)
    
    if not rnaToProteinRefseqIDMapFile.is_file():
        print('producing rna to protein RefSeq ID dictionary')
        produce_rnaToProtein_refseqID_dict (refseqIDFile, rnaToProteinRefseqIDMapFile)
    
    if not seqresFile.is_file():
        print('reducing headers in PDB chain sequence file')
        reduce_fasta_headers (pdbSeqresFile, ' ', 1, 1, seqresFile)

if __name__ == "__main__":
    main()
