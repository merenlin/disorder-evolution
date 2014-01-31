import disprotevol as disevol 

#
#TODO:
#process arguments
#flag 0 - fata source : disprot
#flag 1 - prepare data
#flag 2 - disprot summary 
#flag 3 

def prepare_disprot(source = "disprot"):
    #disevol.dataprocessor.separateDisprot() 
    #disevol.dataprocessor.createProteinsTable()
    #disevol.dataprocessor.createDisorderFasta()
    #disevol.dataprocessor.getRefSeqsMOBIDB("data/raw/mobidb/sequences_longdisorder.fasta")
    pass

def plot_figures():
    disevol.visualizing.contactmaps.plotmaps()

def prepare_pdb():
    disevol.pdb.get_pdb_from_uniprot()
    
if __name__ == '__main__':
    plot_figures()
    