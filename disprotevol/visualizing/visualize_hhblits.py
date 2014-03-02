# run hhblits
import os
import pylab 
from Bio import SeqIO
import matplotlib.pyplot as plt
from csb.bio.io.hhpred  import HHOutputParser

# get all protein names from disprot dir 
inputdir = 'data/hhblits/'
hhrdir = os.getcwd() + '/data/disprot/'
alignments = (os.listdir(inputdir))  #all alignments in the folder

a_lengths = []
neffs = []    #diversity of the alignment
for file_name in alignments: 
    alignment = list(SeqIO.parse(open(inputdir + file_name), "fasta"))
    a_len = len(alignment)
    a_lengths.append(a_len)
    if a_len > 2000:
        try:
            hhr = HHOutputParser().parse_file(hhrdir + file_name[:-5] + ".hhr")
            neffs.append(hhr.neff)
        except Exception as exc:
            print file_name
            print str(exc)

plt.figure()
plt.hist(a_lengths, color="purple", alpha= 0.3)
plt.suptitle("Disprot proteins homologs counts")
plt.xlabel('Number of homologs')
plt.ylabel('Frequency')
plt.savefig('hhblits_num_orthologs_hist.png')

print neffs
plt.figure()
plt.hist(neffs, color="purple", alpha= 0.3)
plt.suptitle("Alignment diversity for > 2000 orthologs")
plt.xlabel('NEFF')
plt.ylabel('Frequency')
plt.savefig('hhblits_neffs_hist.png')

print ">3000 " + str(len([x for x in a_lengths if x > 3000]))
print ">2000 " + str(len([x for x in a_lengths if x > 2000]))
print ">1000 " + str(len([x for x in a_lengths if x > 1000]))