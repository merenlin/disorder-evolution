import os
from contactvis import plot_contact_map

workdir = 'data/pconsc/'
resultdir = 'results/contacts'
c_filename = 'pconsc.slurm.out'

def plotmaps():
    
    # Get all the subdirectories and a protein id from the 
    # directory id 'DP00006/DP00006.fa'
    dirs = [(x[0],x[0].split('/')[-1]) for x in os.walk(workdir)]

    for protein_dir, protein_id in dirs:
        if protein_id != '':    
            fasta_filename = protein_dir + '/' + protein_id + '.fa'
            contact_filename = protein_dir + '/' + c_filename
            pdb_filename = protein_dir + '/' + protein_id + '.pdb'
            out_filename = resultdir + '/'+ protein_id + '.pdf'

            if os.path.exists(pdb_filename): 
                plot_contact_map.plot_map(fasta_filename, contact_filename,
                    pdb_filename = pdb_filename,
                    factor=float(2.0), outfilename=out_filename)
