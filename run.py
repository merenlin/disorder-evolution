import disprotevol as disevol 

#TODO:
#process arguments
#flag 0 - fata source : disprot
#flag 1 - prepare data
#flag 2 - disprot summary 
#flag 3 - disprot alignments summary

def prepare_data(source = "disprot"):
    disevol.dataprocessor.separateDisprot() 

if __name__ == '__main__':
    prepare_data()