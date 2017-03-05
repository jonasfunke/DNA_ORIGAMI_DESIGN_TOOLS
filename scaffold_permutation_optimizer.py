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


#create path to file to load
#current_working_directory = os.getcwd()
#file_name=sys.argv[1]
#file_name=current_working_directory+"/"+file_name
#print('file_name',file_name)
file_name = '/Users/jonasfunke/Dropbox/FRET_STAGE/Designs/FS-v6_spectrometer/FS-v6_009_build/FS-v6_009_03_recode_03.json'
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

def get_scaffold_indices(dna_structure, reportcolor):
    # compile a list of scaffold indices coressponding to staples that were colored with reportcolor
    scaffold_indices =[]
    scaffold_bases_init =[]
    for strand in dna_structure.strands:
        #print str(strand.color)
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
    return scaffold_indices

def get_scaffold_sequences(dna_structure, scaffold_indices, scaffold_rotation):
    scaffold_length = len(dna_structure.strands[scaffold_id].tour)
    scaffold_sequences = []
    for strand in scaffold_indices:
        tmp = []
        for baseindex in strand:
            tmp.append(dna_structure.strands[scaffold_id].tour[(baseindex+scaffold_rotation)%scaffold_length].seq)
        scaffold_sequences.append(tmp)
    return scaffold_sequences
#%%
# Read cadnano file and create dna structure.
converter = read_file(file_name, seq_name)
dna_structure = converter.dna_structure

#determine scaffold strand id
for strand in dna_structure.strands:
    if strand.is_scaffold:
        scaffold_id = strand.id

scaffold_length = len(dna_structure.strands[scaffold_id].tour)

#%% Count thymines
scaffold_indices_black = get_scaffold_indices(dna_structure, [0, 0, 0])
        
print scaffold_indices_black
print 'Number of black oligos: '+str(len(scaffold_indices_black))


# loop through scaffold permutations
permutations_black = []
for i in range(0,scaffold_length):
    scaffold_sequences = get_scaffold_sequences(dna_structure, scaffold_indices_black, i)
    if i%1000 is 0: print 'Permutation '+str(i)    
    
    tmp = [i]
    # compute numerb of thymines
    for strand in scaffold_sequences:
        tmp.append(strand.count('T'))
    permutations_black.append(tmp)
    
file_out =  file_name[:-5]+'_permutations_thymine-count.txt'
a = np.matrix(permutations_black)
np.savetxt(file_out, a, fmt='%i', delimiter='\t', newline='\n') 



#%% Measure melting temperature of red staples
scaffold_indices_red = get_scaffold_indices(dna_structure, [1, 0, 0])
print 'Number of red oligos: '+str(len(scaffold_indices_red))

# loop through scaffold permutations
permutations_red = []
for i in range(0,scaffold_length):
    scaffold_sequences = get_scaffold_sequences(dna_structure, scaffold_indices_red, i)
    if i%1000 is 0: print 'Permutation '+str(i)
    
    tmp = [i]
    # compute melting temp of each strand
    for strand in scaffold_sequences:
        tmp.append(mt.Tm_NN(Seq(''.join(strand)), nn_table=mt.DNA_NN4))
    permutations_red.append(tmp)
    
    
file_out =  file_name[:-5]+'_permutations_melting-temp.txt'
a = np.matrix(permutations_red)
np.savetxt(file_out, a, fmt='%i', delimiter='\t', newline='\n') 


#%% Measure melting temperature of orange staples
scaffold_indices_orange = get_scaffold_indices(dna_structure, [1.0, 0.5019607843137255, 0.0])
print 'Number of orange oligos: '+str(len(scaffold_indices_orange))

# loop through scaffold permutations
permutations_orange = []
for i in range(0,scaffold_length):
    scaffold_sequences = get_scaffold_sequences(dna_structure, scaffold_indices_orange, i)
    if i%1000 is 0: print 'Permutation '+str(i)
    
    tmp = [i]
    # compute melting temp of each strand
    for strand in scaffold_sequences:
        tmp.append(mt.Tm_NN(Seq(''.join(strand)), nn_table=mt.DNA_NN4))
    permutations_orange.append(tmp)

file_out =  file_name[:-5]+'_permutations_melting-temp_orange.txt'
a = np.matrix(permutations_orange)
np.savetxt(file_out, a, fmt='%i', delimiter='\t', newline='\n') 


#%% REPORT CURRENT PERMUTATION
print 'STATS OF PERMUTATION '+str(0)
print 'Thymine content of scaffold loops: {}'.format(sum(permutations_black[0][1:]))+' '+str(permutations_black[0][1:])
print 'Minimum melt. temp of red staples: {}'.format(round(min(permutations_red[0][1:])))
print 'Minimum melt. temp of orange staples: {}'.format(round(min(permutations_orange[0][1:])))
print 'Red melting temperatures ' +str([round(melt,0) for melt in permutations_red[0][1:]])
print 'Orange melting temperatures ' +str([round(melt,0) for melt in permutations_orange[0][1:]])

#%% REPORT ith PERMUTATION
i = 6728
print 'STATS OF PERMUTATION '+str(i)
print 'Thymine content of scaffold loops: {}'.format(sum(permutations_black[i][1:]))+' '+str(permutations_black[i][1:])
print 'Minimum melt. temp of red staples: {}'.format(round(min(permutations_red[i][1:])))
print 'Minimum melt. temp of orange staples: {}'.format(round(min(permutations_orange[i][1:])))
print 'Red melting temperatures ' +str([round(melt,0) for melt in permutations_red[i][1:]])
print 'Orange melting temperatures ' +str([round(melt,0) for melt in permutations_orange[i][1:]])

#%% REPORT BLACK sequences
i = 0
scaffold_sequences = get_scaffold_sequences(dna_structure, scaffold_indices_black, i)
for strand in scaffold_sequences:
    print ''.join(strand)+str(' ')+str(strand.count('T'))


#%% determine scaffold sequence of best permutation
i = 6728
print 'Helix: ' + str(dna_structure.strands[scaffold_id].tour[(-i)%scaffold_length].h)
print 'Index: ' + str(dna_structure.strands[scaffold_id].tour[(-i)%scaffold_length].p)
tmp= []
for j in range(0,20):
    tmp.append(dna_structure.strands[scaffold_id].tour[(-i+j)%scaffold_length].seq)
    
print ''.join(tmp)



#%%
scaffold_sequences = get_scaffold_sequences(dna_structure, scaffold_indices_red, 0)
for strand in scaffold_sequences:
    #tmp.append(mt.Tm_NN(Seq(''.join(strand))))
    print str( round(mt.Tm_NN(Seq(''.join(strand)))) )  + ' ' + str(round(mt.Tm_NN(Seq(''.join(strand)), nn_table=mt.DNA_NN4) ))   + ' ' + ''.join(strand)




