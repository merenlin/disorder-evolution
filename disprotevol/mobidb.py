import urllib2
import json
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord 

species = 'Saccharomyces+cerevisiae'
workingdir = "data/mobidb/Saccharomyces+cerevisiae/"

#u'start': 849, u'ann': u'd', u'end': 850}, {u'start': 851, u'ann': u's', u'end': 851}
def get_proteins(query = 'organism:'+species):
    url = 'http://mobidb.bio.unipd.it/ws/search?q=' + query
    response = json.load(urllib2.urlopen(url))

    table_file = open('data/tables/mobidb/all_' + species + 'human.csv', "w+")
        
    table_file.write("uniprotID,name,species,numdisorder,seqlength\n")

    for s in response: 
        
        url2 = 'http://mobidb.bio.unipd.it/ws/entries/'+ s['acc']+'/disorder'
        disorder = json.load(urllib2.urlopen(url2))
    
        disnum = 0
        full_consensus = disorder['consensus']['full']
        for entry in full_consensus: 
            if entry['ann'].lower() == 'd':
                disnum += (entry['end']-entry['start'])+1
        
        if (disnum > 30):
            seqlength = get_fasta(s['acc'], disorder)  
            table_file.write(s['acc'] + ',"' + s['pname'] + '",' + s['org'] +
                ',' + str(disnum) + ',' + str(seqlength) +'\n')
        
    table_file.close()  

"""
Takes uniprot id, save a fasta file with a sequence
and returns sequence length
"""
def get_fasta(proteinID, disorder):
    url = "http://mobidb.bio.unipd.it/ws/entries/"+ proteinID +"/uniprot"
    data = json.load(urllib2.urlopen(url))

    protein = SeqRecord(Seq(data["sequence"]),id=proteinID,description="")
    
    proteindir = "data/mobidb/"+ species +"/" + proteinID+"/"
    if not os.path.exists(proteindir):
        os.makedirs(proteindir)
    else:
        return len(data["sequence"])

    f_result = open(proteindir+proteinID+".fa", "w+")
    SeqIO.write(protein, f_result, "fasta")
    f_result.close()

    with open(proteindir+proteinID+"_disorder.json", "w+") as f_dis:
        json.dump(disorder, f_dis) 

    return len(data["sequence"])


def getHomologsInfo():
    homologs = {}
    for proteindir in os.listdir(workingdir):
        try:
            alignment = list(SeqIO.parse(open(proteindir + "/"+ proteindir +".oa3m"), "fasta"))
        except Exception as exc:
             print exc
             return homologs

    homologs = len(alignment)

    table_file = open('data/tables/mobidb/' + species + '_homologs.csv', "w+")
    table_file.write("uniprotID,numhomologs\n")
    for i in homologs: 
        table_file.write(i)
    table_file.close() 

    """ sumgaps = 0
    numdis = disseq.count(0)
    for i in range(0,len(disseq)): 
        numgaps = 0
        if disseq[i] == 0:
            for seq in alignment:
                if (i<len(seq)) & (seq[i] == '-'): numgaps+=1
        if numdis != 0:
            sumgaps += numgaps/disseq.count(0)

    homologs["gaps"] = sumgaps/len(alignment)
    """
    return homologs    

def annotationToString(entries, flag = 'd'):
    ids = []
    for entry in entries: 
        if entry['ann'].lower() == flag:
            for i in xrange(int(entry['start']), int(entry['end'])+1): 
                ids.append(i)
    return ids

def classifyDisorder():
     for proteindir in os.listdir(workingdir):
        with open(workingdir + proteindir+"/"+proteindir+"_disorder.json") as f_dis:
            disorder = json.load(f_dis) 
            
            disprot = annotationToString(disorder['consensus']['disprot'])
            pdb = annotationToString(disorder['consensus']['pdb']) #just disorder
            pdbs = annotationToString(disorder['consensus']['pdb'],'s') #just order
            
            pred = annotationToString(disorder['consensus']['predictors'])
            
            disseq = ""
     
            full_consensus = disorder['consensus']['full']
            for entry in full_consensus: 
                if entry['ann'].lower() == 'd':
                    for i in xrange(int(entry['start'])-1, int(entry['end'])): 
                        if ((i+1) in disprot): 
                            if ((i+1) in pdb):
                                disseq += '0'
                            elif ((i+1) in pdbs):
                                disseq += '1'
                            else:
                                disseq += '0'
                        else: 
                            if ((i+1) in pdb):
                                disseq += '4'
                            elif ((i+1) in pred):
                                disseq += '5'
                else:
                    for i in xrange(int(entry['start'])-1, int(entry['end'])): 
                        disseq += '1'
           
            fout = open(workingdir + proteindir+"/"+proteindir+"_disorder.fasta", "w+")       
            record = SeqRecord(id = proteindir,seq = Seq(disseq))

            SeqIO.write(record, fout, "fasta")
            fout.close()
    

def createRegionsTable():
    table_file = open('data/tables/mobidb/' + species + '_disorder.csv', "w+")
    table_file.write("uniprotID,disnum,pdbnum,prednum\n")
    for proteindir in os.listdir(workingdir):
        f_dis = open(workingdir + proteindir+"/"+proteindir+"_disorder.fasta")
        seqs = SeqIO.parse(f_dis,"fasta")
        for record in seqs:
            disseq = [int(i) for i in record.seq] 
            disnum = disseq.count(0)
            pdbnum = disseq.count(4)
            prednum = disseq.count(5)
        disannot = proteindir + ',' + str(disnum) + ',' + str(pdbnum) + ',' + str(prednum) +"\n"
        table_file.write(disannot)
    table_file.close()
   
