#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 09:05:04 2017

@author: jonasfunke
"""
#%%
import os
import sys
import numpy as np
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
#import nanodesign converter module
base_path = '/Users/jonasfunke/NANODESIGN/nanodesign'
sys.path.append( base_path )
import nanodesign
from nanodesign.converters import Converter
sys.path = sys.path[:-1]

#%%

#create path to file to load
#current_working_directory = os.getcwd()
#file_name=sys.argv[1]
#file_name=current_working_directory+"/"+file_name
#print('file_name',file_name)
file_name = '/Users/jonasfunke/Dropbox/FRET_STAGE/Designs/FS-v6_spectrometer/FS-v6_009_build/FS-v6_009_02_recode_04.json'
seq_name = 'p7704'

#create output path
#output_name = file_name + "_gillespieData.txt"
#print('output_name',output_name)

#%% Define functions
def read_file(file_name, seq_name):
    """ Read in a cadnano file. """
    converter = Converter()
    seq_file = None
    converter.read_cadnano_file(file_name, seq_file, seq_name)
    return converter

#%%
# Read cadnano file and create dna structure.
converter = read_file(file_name, seq_name)
dna_structure = converter.dna_structure

#determine scaffold strand id
for strand in dna_structure.strands:
    if strand.is_scaffold:
        scaffold_id = strand.id
        
#%% Count thymines
reportcolor = [0, 0, 0]

# compile a list of scaffold indices coressponding to staples that were colored with reportcolor
scaffold_indices =[]
scaffold_bases_init =[]
for strand in dna_structure.strands:
    #print str(strand.color)+' '+str(num_bases)
    if strand.color==reportcolor:
        #seq = [base.seq for base in strand.tour]
        strand_scaffold_base = []
        strand_scaffold_index = []
        for base in strand.tour:
            if base.across is not None:
                strand_scaffold_base.append(base.across.seq)
                strand_scaffold_index.append(dna_structure.strands[scaffold_id].get_base_index(base.across))
                
        strand_scaffold_base.reverse() # reverse the squence to obtain 5 - 3 p 
        strand_scaffold_index.reverse() # reverse the sequence to obtain 5' to 3'        
        scaffold_indices.append(strand_scaffold_index)
        scaffold_bases_init.append(strand_scaffold_base)
        
print scaffold_indices
print len(scaffold_indices)

scaffold_length = len(dna_structure.strands[scaffold_id].tour)


# loop through scaffold permutations

permutations = []
print scaffold_bases_init
for i in range(0,scaffold_length):
    scaffold_sequences = []
    for strand in scaffold_indices:
        tmp = []
        for baseindex in strand:
            tmp.append(dna_structure.strands[scaffold_id].tour[(baseindex+i)%7704].seq)
        scaffold_sequences.append(tmp)
    print 'Iteration '+str(i)
    
    #print scaffold_sequences
    
    tmp = [i]
    # compute numerb of thymines
    for strand in scaffold_sequences:
        tmp.append(strand.count('T'))
    permutations.append(tmp)
    
file_out =  file_name[:-5]+'_permutations_thymine-count.txt'
a = np.matrix(permutations)
np.savetxt(file_out, a, fmt='%i', delimiter='\t', newline='\n') 


#%%
# REPORT 
i = 0
scaffold_sequences = []
for strand in scaffold_indices:
    tmp = []
    for baseindex in strand:
        tmp.append(dna_structure.strands[scaffold_id].tour[(baseindex+i)%7704].seq)
        
    tmp.reverse()
    print ''.join(tmp)



#%% Measure melting temperature 
reportcolor = [1, 0, 0]

# compile a list of scaffold indices coressponding to staples that were colored with reportcolor
scaffold_indices =[]
scaffold_bases_init =[]
for strand in dna_structure.strands:
    #print str(strand.color)+' '+str(num_bases)
    if strand.color==reportcolor:
        #seq = [base.seq for base in strand.tour]
        strand_scaffold_base = []
        strand_scaffold_index = []
        for base in strand.tour:
            if base.across is not None:
                strand_scaffold_base.append(base.across.seq)
                strand_scaffold_index.append(dna_structure.strands[scaffold_id].get_base_index(base.across))
                
        strand_scaffold_base.reverse() # reverse the squence to obtain 5 - 3 p 
        strand_scaffold_index.reverse() # reverse the sequence to obtain 5' to 3'        
        scaffold_indices.append(strand_scaffold_index)
        scaffold_bases_init.append(strand_scaffold_base)
        
print scaffold_indices
print len(scaffold_indices)

scaffold_length = len(dna_structure.strands[scaffold_id].tour)


# loop through scaffold permutations

permutations = []
print scaffold_bases_init
for i in range(0,scaffold_length):
    scaffold_sequences = []
    for strand in scaffold_indices:
        tmp = []
        for baseindex in strand:
            tmp.append(dna_structure.strands[scaffold_id].tour[(baseindex+i)%7704].seq)
        scaffold_sequences.append(tmp)
    if i%1000 is 0:
        print 'Iteration '+str(i)
    
    tmp = [i]
    # compute melting temp of each strand
    for strand in scaffold_sequences:
        tmp.append(mt.Tm_NN(Seq(''.join(strand))))
    permutations.append(tmp)
    
file_out =  file_name[:-5]+'_permutations_melting-temp.txt'
a = np.matrix(permutations)
np.savetxt(file_out, a, fmt='%i', delimiter='\t', newline='\n') 




#%% Measure melting temperature 
reportcolor = [1.0, 0.5019607843137255, 0.0]

# compile a list of scaffold indices coressponding to staples that were colored with reportcolor
scaffold_indices =[]
scaffold_bases_init =[]
for strand in dna_structure.strands:
    #print str(strand.color)+' '+str(num_bases)
    if strand.color==reportcolor:
        #seq = [base.seq for base in strand.tour]
        strand_scaffold_base = []
        strand_scaffold_index = []
        for base in strand.tour:
            if base.across is not None:
                strand_scaffold_base.append(base.across.seq)
                strand_scaffold_index.append(dna_structure.strands[scaffold_id].get_base_index(base.across))
                
        strand_scaffold_base.reverse() # reverse the squence to obtain 5 - 3 p 
        strand_scaffold_index.reverse() # reverse the sequence to obtain 5' to 3'        
        scaffold_indices.append(strand_scaffold_index)
        scaffold_bases_init.append(strand_scaffold_base)
        
print scaffold_indices
print len(scaffold_indices)

scaffold_length = len(dna_structure.strands[scaffold_id].tour)


# loop through scaffold permutations

permutations = []
print scaffold_bases_init
for i in range(0,scaffold_length):
    scaffold_sequences = []
    for strand in scaffold_indices:
        tmp = []
        for baseindex in strand:
            tmp.append(dna_structure.strands[scaffold_id].tour[(baseindex+i)%7704].seq)
        scaffold_sequences.append(tmp)
    if i%1000 is 0:
        print 'Iteration '+str(i)
    
    tmp = [i]
    # compute melting temp of each strand
    for strand in scaffold_sequences:
        tmp.append(mt.Tm_NN(Seq(''.join(strand))))
    permutations.append(tmp)
    
file_out =  file_name[:-5]+'_permutations_melting-temp_orange.txt'
a = np.matrix(permutations)
np.savetxt(file_out, a, fmt='%i', delimiter='\t', newline='\n') 




#%% determine scaffold sequence of best permutation
i = 3128
tmp= []
for j in range(0,10):
    tmp.append(dna_structure.strands[scaffold_id].tour[-i+j].seq)
    
print ''.join(tmp)


#%%
scaffold_sequences = []
for strand in scaffold_indices:
    tmp = []
    for baseindex in strand:
        tmp.append(dna_structure.strands[scaffold_id].tour[(baseindex+i)%7704].seq)
    scaffold_sequences.append(tmp)
print 'Iteration '+str(i)




