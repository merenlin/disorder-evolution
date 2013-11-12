# separate Disprot proteins into different
# fasta files
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
from csb.bio.io.hhpred  import HHOutputParser
import os
import re

#TODO: move all the paths to a config file
#configs
FASTfilename = "data/raw/disprot.fasta"   
MOBIDBseqfile = "data/raw/mobidb/sequences.fasta"
MOBIDBannotfile = "data/raw/mobidb/annotations.fasta"
HHsearchdir = "data/hhsuite/hhsearch/"
tblProteinFilename = "data/tables/proteins.txt"

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

Takes only reference sequences from this file 
and returns their annotations as a dictionary.

"""
def getRefSeqsMOBIDB(file = MOBIDBseqfile):
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
def getDisorderDisprot(length, s = '&78-90 #216-261 #1-7'):
    disorder = [0] * length # no disorder by default
    regions = s.split()
    for region in regions:
        print region
        if region[0]=='#':
            region = region[1:].split('-')
            for i in xrange(int(region[0]), int(region[1])): 
                disorder[i] = 1             
    return disorder.count(1)

"""
Input: *.hhr file for a protein in DisProt
Output: number of homologs and NEFF(alignment diversity score)

"""
def getHomologsInfo(dp_id):
    homologs = {}
    try:
        hhr = HHOutputParser().parse_file(HHsearchdir + dp_id + ".fa.hhr")
        homologs["neff"]=hhr.neff
        homologs["num"]=hhr.no_of_seqs.split()[0]
    except Exception as exc:
        print dp_id
        print str(exc)
        pass
    return homologs

#uniprotid / sourceID /  name / organism / length / disorder content 
def createProteinsTable():
    f = open(FASTfilename,"rU")
    seqs = SeqIO.parse(f,"fasta")

    mobidb_seqs =  getRefSeqsMOBIDB(MOBIDBseqfile)
   
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
                     "numhomologs,alignmentdiv\n")

    # Annotations string looks like this: 
    # DisProt|DP00001|uniprot|Q9HFQ6|sp|RLA3_CANAL #1-108
    for record in seqs: 
        annotations = record.description.split("|")
        disprot_id = annotations[1] 
        uniprot_id = annotations[3]
        seqlength = len(record.seq)

        disorder_residues = annotations[len(annotations)-1]
        disnum = getDisorderDisprot(seqlength, disorder_residues)
        discontent = round((disnum/float(seqlength)) * 100,1)

        # Integrate with MOBIDB information where possible
        try:
            name = mobidb_seqs[uniprot_id]['name']
            species = mobidb_seqs[uniprot_id]['species']
            species = species.split()[0] + ' ' + species.split()[1]
        except KeyError, e:
            species = ''
            name = ''

        try:
            # Get homologs information from hhr files
            hominfo = getHomologsInfo(disprot_id) 
            print hominfo
            neff = hominfo['neff']
            homnum = hominfo['num']
        except KeyError, e:
             neff = ''
             homnum = ''

        table_file.write(disprot_id + ',' 
            + uniprot_id + ',' 
            + '"' + name + '",' 
            + species + ','  
            + str(seqlength) + ',' 
            + str(disnum) + ',' 
            + str(discontent) + ',' 
            + homnum + ',' 
            + str(neff) + '\n')  

    table_file.close() 

#TODO: make a separate table for homnums with different cuttoffs
