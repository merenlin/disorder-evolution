# separate Disprot proteins into different
# fasta files
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
import re

#TODO: move all the paths to a config file
def separateDisprot():
    FASTfilename = "data/raw/disprot.fasta"   
    f = open(FASTfilename,"rU")
    seqs = SeqIO.parse(f,"fasta")

    for record in seqs: 
        disprot_id = record.id.split("|")[1]  #parse out the individual id
        protein = SeqRecord(record.seq,id=disprot_id,description="")
        f_result = open("data/disprot/"+disprot_id+".fasta", "w+")
        SeqIO.write(protein, f_result, "fasta")
        f_result.close()
        
    f.close()
