#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:42:58 2017

@author: jonasfunke
    
    This script computes the sequences of the black staples (colored in cadnano version 1) for each scaffold permutation
    usage: python scaffold_permutations.py path/to/json_file.json p7560
"""

# Imports
import os
import csv
import __future__

# try to impot nanodesign package. I assume the script is in x/somename/design_statistics.py and the nanodesign package is in x/nanodesign
try:
    import nanodesign
except ImportError:
    import sys
    #base_path = '/Users/jonasfunke/NANODESIGN/nanodesign'
    base_path = os.path.abspath( os.path.join( os.path.dirname(os.path.abspath( __file__)), '../nanodesign'))
    sys.path.append(base_path)
    import nanodesign
    # If the import fails now, we let the exception go all the way up to halt execution.
    sys.path = sys.path[:-1]    
    from nanodesign.converters import Converter

# Functions
def read_file(file_name, seq_name): # read file
    converter = Converter()
    seq_file = None
    converter.read_cadnano_file(file_name, seq_file, seq_name)
    return converter


def get_scaffold_indices(dna_structure, reportcolor, scaffold_id, physical_index):
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
                if base.seq is not 'N':
                    strand_scaffold_base.append(base.across.seq)
                    strand_scaffold_index.append(physical_index[dna_structure.strands[scaffold_id].get_base_index(base.across)])
                    
            strand_scaffold_base.reverse() # reverse the squence to obtain 5 - 3 p 
            strand_scaffold_index.reverse() # reverse the sequence to obtain 5' to 3'        
            scaffold_indices.append(strand_scaffold_index)
            scaffold_bases_init.append(strand_scaffold_base)
    #for strand in scaffold_bases_init:
    #    print(''.join(strand))
    return scaffold_indices

def get_sequence(strand):
    cur_seq = []
    for base in strand.tour:
        cur_seq.append(base.seq)
    return ''.join(cur_seq)

def get_scaffold_sequences(dna_structure, scaffold_indices, scaffold_rotation, scaffold_id):
    # physical scaffold sequence
    scaffold_sequence = get_sequence(dna_structure.strands[scaffold_id]).replace('N', '')
    # physical scaffold length
    scaffold_length = len(scaffold_sequence)
    
    #loop through strands
    scaffold_sequences = []
    for strand in scaffold_indices:
        # loop through bases in strand
        tmp = []
        for baseindex in strand:
            # physical index in scaffold
            i_physical = (baseindex+scaffold_rotation)%scaffold_length
                         
            #dna_structure.strands[scaffold_id].tour[i_physical+offset].seq
            tmp.append(scaffold_sequence[i_physical])
        scaffold_sequences.append(''.join(tmp))

    return scaffold_sequences

def get_staple_sequences_from_scaffold_indices(dna_structure, scaffold_indices, scaffold_rotation, scaffold_id):
    scaffold_sequences = get_scaffold_sequences(dna_structure, scaffold_indices, scaffold_rotation, scaffold_id)
    
    staple_sequences = []
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for strand in scaffold_sequences:
        staple_sequences.append("".join(complement.get(base, base) for base in reversed(strand)))
        #print strand
        #print "".join(complement.get(base, base) for base in reversed(strand))
    return staple_sequences


def get_index_lists(dna_structure, scaffold_id):
    scaffold_sequence = get_sequence(dna_structure.strands[scaffold_id])
    
    physical_to_design = []
    design_to_physical =[]
    
    physical_index = 0
    for design_index in range(len(scaffold_sequence)):
        #physical_index.append(i-scaffold_sequence[0:design_index+1].count('N'))
        if scaffold_sequence[design_index] is 'N':
            design_to_physical.append(physical_index)
        else:
            physical_to_design.append(design_index)
            design_to_physical.append(physical_index)
            physical_index = physical_index+1
    return physical_to_design, design_to_physical
    
    
#%%
def main():
    # Parse arguments, TODO: use parser object
    if (len(sys.argv) != 3):
        sys.stderr.write("**** ERROR: Wrong number of arguments.\n") 
        sys.stderr.write("Usage: design_statistics.py <filename> scaffold_sequence\n")
        sys.stderr.write("Output: If <filename> is path/name.json, output will be placed in path/name_analysis/\n")
        sys.exit(1)
        
    file_full_path_and_name = os.path.abspath( os.path.expanduser( sys.argv[1] ))
    seq_name = sys.argv[2]

    #%%
    #file_full_path_and_name = '/Users/jonasfunke/Dropbox/FRET_STAGE/Designs/FS-v6_spectrometer/twist_screen/FS-v6_019_deacivated.json'
    #seq_name = 'p7704'
    
    # parse filename and create output directory
    file_name = os.path.basename( file_full_path_and_name )
    output_path = os.path.dirname(file_full_path_and_name) + os.sep + file_name[:-5] + '_scaffold-permutations' + os.sep
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = output_path+file_name[:-5]
    
    # Read cadnano file and create dna structure.
    converter = read_file(file_full_path_and_name, seq_name)
    dna_structure = converter.dna_structure     
    dna_structure.compute_aux_data() # compute domain data
        
   #determine scaffold strand id
    for strand in dna_structure.strands:
        if strand.is_scaffold:
            scaffold_id = strand.id
    
    #%% make the maps, that map design to physical indeces
    design_index, physical_index = get_index_lists(dna_structure, scaffold_id)
     
    #%% Loop through scaffold permutations and write sequencs of black oligos 
    scaffold_indices_black = get_scaffold_indices(dna_structure, [0.2, 0.2, 0.2], scaffold_id, physical_index)       
    #print scaffold_indices_black
    print('Number of black oligos: '+str(len(scaffold_indices_black)))
    
    # print sequences of oligos
    print('Complementary sequences of black oligos for current scaffold permutation:')
    tmp = get_scaffold_sequences(dna_structure, scaffold_indices_black, 0, scaffold_id)
    for strand in tmp:
        print(strand)
        
    #%% Loop through scaffold permutations and write sequencs of black oligos 
    #sequence of scaffold in design, this includes skips as 'N' 
    design_scaffold_sequence = get_sequence(dna_structure.strands[scaffold_id])
    #physical sequence of scaffold, this is the true scaffold sequence 
    physical_scaffold_length = len(design_scaffold_sequence)-design_scaffold_sequence.count('N')
    design_scaffold_length = len(design_scaffold_sequence)
    
    file_out = output_path+'_permutations_black-sequences.csv'
    with open(file_out, 'wb') as csvfile:
        outputwriter = csv.writer(csvfile, delimiter=';')
        tmp = ['Permutation', 'Start helix', 'Start position']
        for j in range(len(scaffold_indices_black)):
            tmp.append('Black oligo ' + str(j))
        outputwriter.writerow(tmp)
        
        # loop through scaffold permutations in physical space
        for i in range(0,physical_scaffold_length):
            if i%1000 is 0: print('Computing permutation '+str(i) + ' of ' + str(physical_scaffold_length))
            
            # determine start of scaffold in design
            i_d = (-design_index[i])%design_scaffold_length

            # sequence of black oligos of rotation i
            staple_sequences = get_staple_sequences_from_scaffold_indices(dna_structure, scaffold_indices_black, i, scaffold_id) 
            
            # write data for this rotation
            tmp = [str(i), str(dna_structure.strands[scaffold_id].tour[i_d].h), str(dna_structure.strands[scaffold_id].tour[i_d].p)]
            for strand in staple_sequences:
                tmp.append(strand)
            outputwriter.writerow(tmp)
        print('Output written to: ' + file_out)
#%%

if __name__ == '__main__':
    main()