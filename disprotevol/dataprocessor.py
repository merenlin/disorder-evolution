# separate Disprot proteins into different
# fasta files
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
import re

#TODO: move all the paths to a config file
#configs
FASTfilename = "data/raw/disprot.fasta"   
MOBIDBseqfile = "data/raw/mobidb/sequences.fasta"
MOBIDBannotfile = "data/raw/mobidb/annotations.fasta"

#TODO: switch to pyfasta for this instead
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
Format of the file

><Sequence ID>|refseq|<Protein name>|<Organism>

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
Disordered regions are denoted by symbol "#" in format:

#<starting residue>-<ending residue>

In the example above residues 1 to 195, 217 to 222, 260-276, and 268 to 271 are disordered.

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

#uniprotid / sourceID /  name / organism / length / disorder content 
def createProteinsTable():
    #TODO: remove table file before starting
    f = open(FASTfilename,"rU")
    seqs = SeqIO.parse(f,"fasta")

    mobidb_seqs =  getRefSeqsMOBIDB(MOBIDBseqfile)
   
    f = open(MOBIDBannotfile,"rU")
    mobidb_annot = SeqIO.to_dict(SeqIO.parse(f,"fasta"))

    table_file = open("data/tables/proteins.txt", "w+")

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

        table_file.write(disprot_id + ',' 
            + uniprot_id + ',' 
            + '"' + name + '",' 
            + species + ','  
            + str(seqlength) + ',' 
            + str(disnum) + ',' 
            + str(discontent) + '\n')  

        # TODO: add NEFF and num_orthologs
    table_file.close() 
