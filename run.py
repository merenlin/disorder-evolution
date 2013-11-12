import disprotevol as disevol 

#TODO:
#process arguments
#flag 0 - fata source : disprot
#flag 1 - prepare data
#flag 2 - disprot summary 
#flag 3 

def prepare_data(source = "disprot"):
    #disevol.dataprocessor.separateDisprot() 
    disevol.dataprocessor.createProteinsTable()

if __name__ == '__main__':
	#if args[2] == "data":
    prepare_data()
    