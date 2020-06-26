#!/usr/bin/env python
# coding: utf-8

# # Rosalind Problems

# ### Rosalind_Problems are problems are found at rosalind.info, a platform for learning bioinformatics through problem solving. The problems below deal with analyzing sequencing data using python, while additional problems dealing with bioinformatic problems can be found at https://github.com/agalianese/Pytnon_Code.

# In[160]:


import Bio
from Bio import Entrez
from Bio import SwissProt
from Bio import ExPASy
from Bio import SeqIO


# In[113]:


#Genbank Introduction

#read in file and set variable names
with open('C:/Users/mc/Downloads/rosalind_gbk.txt') as data_in:
    genus, begin_date, end_date = [line.strip() for line in data_in.readlines()]
        
#contact email
Entrez.email = 'amgalianese@gmail.com'

#input search variables of the genus, the start date, and the end date
handle = Entrez.esearch(db='nucleotide', term=genus+'[Organism]', mindate=begin_date, maxdate = end_date, datetype = 'pdat')

#returns which organisms fit this profile
record = Entrez.read(handle)
    
#return the amount of organisms 
print(record["Count"])


# In[ ]:





# In[157]:


#Data Formats


def short_seq(ids):
    """Function that takes in a list of ids, find the one with the shortest sequence length and returns its description
    and sequence"""
    #contact email
    Entrez.email = "amgalianese@gmail.com"

    #ids being searched using fasta format
    handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta")

    #create a list of the the SeqIO objects in Fasta format
    records = list(SeqIO.parse(handle, "fasta")) 

    #begin with the shortest length and id to be compared
    shortest = 1000000000000
    short_id = ''

    #loop over the sequences to find the shortest, then log those results
    for record in records:
        if len(record) < shortest:
            shortest = len(record)
            short_id = record
    #prints the description of the shortest id and its fasta format sequence
    print('>' + short_id.description)
    print(short_id.seq)


short_seq(["NM_001251956, JX445144, NM_001159758, JX475048, NM_131329, JX317622, JQ712982, NM_001133698, JX491654, NM_001197168"])


# In[ ]:




