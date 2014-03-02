import disprotevol as disevol 

#
#TODO:
#process arguments
#flag 0 - fata source : disprot
#flag 1 - prepare data
#flag 2 - disprot summary 
#flag 3 

def prepare_disprot(source = "disprot"):
    disevol.dataprocessor.separateDisprot() 
    disevol.dataprocessor.getRefSeqsMOBIDB("data/raw/mobidb/sequences_longdisorder.fasta")
    disevol.dataprocessor.createDisorderFasta()
    disevol.dataprocessor.createProteinsTable()
    
def plot_figures():
    disevol.visualizing.contactmaps.plotmaps()

def prepare_pdb():
    disevol.pdb.get_pdb_from_uniprot()

def prepare_mobidb():
    disevol.mobidb.get_proteins()
    disevol.mobidb.getHomologsInfo()

if __name__ == '__main__':
    disevol.mobidb.classifyDisorder()
    disevol.mobidb.createRegionsTable()
    