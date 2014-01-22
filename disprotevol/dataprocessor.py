# separate Disprot proteins into different
# fasta files
from __future__ import division

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord 
from csb.bio.io.hhpred  import HHOutputParser
import logging

import os
import re

#TODO: move all the paths to a config file
#configs
FASTfilename = "data/raw/disprot.fasta"   
MOBIDBseqfile = "data/raw/mobidb/sequences.fasta"
MOBIDBannotfile = "data/raw/mobidb/annotations.fasta"
HHblitsdir = "data/hhsuite/hhblits/1e-4/"
tblProteinFilename = "data/tables/proteins.txt"
fastaDisorderFilename = "data/tables/disorder.fasta"

"""
Takes disprot.fasta and splits it into 
individual fasta files for each protein.

Throws away all the annotation information
and keeps only the disprot id. 

"""
def separateDisprot():
    f = open(FASTfilename,"rU")
    seqs = SeqIO.parse(f,"fasta")

    for record in seqs: 
        annotations = record.id.split("|")
        disprot_id = annotations[1] 
        protein = SeqRecord(record.seq,id=disprot_id,description="")
        f_result = open("data/disprot/"+disprot_id+".fasta", "w+")
        SeqIO.write(protein, f_result, "fasta")
        f_result.close()
        
    f.close()

"""
Input: sequences.fasta from MOBIDB

Format of the file

><Sequence ID>|refseq|<Protein name>|<Organism>

Takes only reference sequences from this file and saves to 
a new fasta file.

"""
def getRefSeqsMOBIDB(file = MOBIDBseqfile):
    f = open(MOBIDBseqfile,"rU")
    mobidb =  SeqIO.parse(f,"fasta")
    fout = open(file + "-refseqs", "w+") #need to remove if exists
    
    for record in mobidb:
        idsplit = record.description.split('|')
        seqtype = idsplit[1]
        if seqtype == 'refseq':
            SeqIO.write(record, fout, "fasta")
    fout.close()
    f.close()
    return
    
"""
Input: sequences.fasta from MOBIDB

Format of the file

><Sequence ID>|refseq|<Protein name>|<Organism>

Takes only reference sequences from this file 
and returns their annotations as a dictionary.

"""
def getRefSeqsInfoMOBIDB(file = MOBIDBseqfile):
    f = open(MOBIDBseqfile,"rU")
    mobidb =  SeqIO.parse(f,"fasta")
    refseqs = {}
    for record in mobidb:
        idsplit = record.description.split('|')
        seqtype = idsplit[1]
        if seqtype == 'refseq':
            refseq = {}
            refseq['name'] = idsplit[2]
            refseq['species'] = idsplit[3]
            uniprot_id = idsplit[0]
            refseqs[uniprot_id] = refseq
    f.close()
    return refseqs


"""
Takes a disorder annotation sequence from Disprot,
transforms to a protein sequence 10100001... 
with 0 standing for order, 1 for disorder.
Function returns the count of 1s in this sequence.

Disordered regions are denoted by symbol "#" in format:

#<starting residue>-<ending residue>

Ordered (structurally determined) parts of proteins are denoted by symbol "&" in format:

<starting residue>-<ending residue>

The functional class(es) and subclass(es) (if known) of each structurally determined region follow
the starting-ending residue. Functional classes are denoted by the symbol "*", functional subclasses
are denoted by the symbol ":". For example:

<starting residue>-<ending residue>*Molecular recognition effectors:Protein-protein binding

"""
def getDisorderSequence(length, s = '&78-90 #216-261 #1-7'):
    disorder = [1] * length # no disorder by default
    regions = s.split()
    for region in regions:
        if region[0]=='#':
            region = region[1:].split('-')
            for i in xrange(int(region[0])-1, int(region[1])): 
                disorder[i] = 0             
    return disorder

"""
Input: *.hhr file for a protein in DisProt
Output: number of homologs and NEFF(alignment diversity score)

"""
def getHomologsInfo(dp_id, disseq):
    homologs = {}
    try:
        alignment = list(SeqIO.parse(open(HHblitsdir + dp_id + ".oa3m"), "fasta"))
    except Exception as exc:
         print exc
         return homologs

    homologs["num"]=len(alignment)

    sumgaps = 0
    numdis = disseq.count(0)
    for i in range(0,len(disseq)): 
        numgaps = 0
        if disseq[i] == 0:
            for seq in alignment:
                if (i<len(seq)) & (seq[i] == '-'): numgaps+=1
        if numdis != 0:
            sumgaps += numgaps/disseq.count(0)

    homologs["gaps"] = sumgaps/len(alignment)
    return homologs

def createDisorderFasta():
    f = open(FASTfilename,"rU")
    fout = open(fastaDisorderFilename, "w+") #need to remove if exists
    seqs = SeqIO.parse(f,"fasta")
    for record in seqs: 
        annotations = record.description.split("|")
        disprot_id = annotations[1] 
        disorder_residues = annotations[len(annotations)-1]
        seqlength = len(record.seq)
        disseq = getDisorderSequence(seqlength, disorder_residues)
        disseq = ''.join(str(x) for x in disseq) #make a string
        record.id =  disprot_id 
        record.seq = Seq(disseq)
        SeqIO.write(record, fout, "fasta")
    fout.close()

#uniprotid / sourceID /  name / organism / length / disorder content 
def createProteinsTable():
    f = open(FASTfilename,"rU")
    seqs = SeqIO.parse(f,"fasta")

    mobidb_seqs =  getRefSeqsInfoMOBIDB(MOBIDBseqfile)
   
    f = open(MOBIDBannotfile,"rU")
    mobidb_annot = SeqIO.to_dict(SeqIO.parse(f,"fasta"))

    #Remove the table file before starting
    try:
        os.remove(tblProteinFilename)
    except OSError:
        pass

    table_file = open(tblProteinFilename, "w+")
    table_file.write("disprotID,uniprotID,name,species,"+
                     "seqlength,numdisorder,discontent," +
                     "numhomologs,gaps\n")

    # Annotations string looks like this: 
    # DisProt|DP00001|uniprot|Q9HFQ6|sp|RLA3_CANAL #1-108
    for record in seqs: 
        annotations = record.description.split("|")
        disprot_id = annotations[1] 
        uniprot_id = annotations[3]
        seqlength = len(record.seq)

        disorder_residues = annotations[len(annotations)-1]
        disseq = getDisorderSequence(seqlength, disorder_residues)
        disnum = disseq.count(0) #0 - for disorder
        discontent = round((disnum/float(seqlength)) * 100,1)

        # Integrate with MOBIDB information where possible
        try:
            name = mobidb_seqs[uniprot_id]['name']
            species = mobidb_seqs[uniprot_id]['species']
            species = species.split()[0] + ' ' + species.split()[1]
        except KeyError, e:
            logging.warning("MobiDB can't find UniprotId " + uniprot_id + ", DisprotID " + disprot_id)
            species = ''
            name = ''

        try:
            # Get homologs information from hhr files
            hominfo = getHomologsInfo(disprot_id, disseq) 
            homnum = hominfo['num']
            homgaps = hominfo['gaps']
        except KeyError, e:
             logging.warning("Could not retrieve homologs info, DisprotID " + disprot_id)
             homnum = ''
             homgaps = ''

        table_file.write(disprot_id + ',' 
            + uniprot_id + ',' 
            + '"' + name + '",' 
            + species + ','  
            + str(seqlength) + ',' 
            + str(disnum) + ',' 
            + str(discontent) + ',' 
            + str(homnum) + ',' + str(homgaps) +'\n')  
        
    table_file.close() 

#TODO: make a separate table for homnums with different cuttoffs
