# BioServices 1: obtain the PDB ID from a given uniprot ID (P43403 i.e. ZAP70)
from bioservices import *
import pandas as pn 

workdir = 'data/pconsc/'

def get_uniprot_from_disprot(disprotID = 'DP00001'):
    ## Get info from proteins table
    df = pn.read_csv('data/tables/proteins.txt',dtype=str)
    row = df[df.disprotID==disprotID]
    return str(row['uniprotID'].values[0])
    
def get_pdb_from_uniprot():
    dirs = [(x[0],x[0].split('/')[-1]) for x in os.walk(workdir)]

    for protein_dir, protein_id in dirs:
        if protein_id != '':    
            #fasta_filename = protein_dir + '/' + protein_id + '.fa'
            uniprotID = get_uniprot_from_disprot(protein_id)

            print("Retrieving Uniprot ID " + uniprotID)
            u = UniProt(verbose=False)
            res = u.mapping(fr="ID", to="PDB_ID", query=uniprotID, format="tab")
            
            if (res!=[]) & (res!={}):
                print res[uniprotID] 
                pdb_id = res[uniprotID][0]   

                # TODO: switch this urllib method to biopython.pdb.fetch
                import urllib
                url = 'http://www.rcsb.org/pdb/files/'+pdb_id+'.pdb'
                print("Fetching PDB file "  + url )
                out = protein_dir + '/' + protein_id + '.pdb'
                urllib.urlretrieve(url, out)
                
               