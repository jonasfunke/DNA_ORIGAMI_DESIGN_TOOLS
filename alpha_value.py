#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 14:51:51 2017

@author: jonasfunke
    This script will return the alpha value for each scaffold permutation
"""

# Imports
import os
import csv
import __future__
from Bio.SeqUtils import MeltingTemp # compute melting temperatures
from Bio.Seq import Seq #create biopython sequences

# try to impot nanodesign package. I assume the script is in x/somename/design_statistics.py and the nanodesign package is in x/nanodesign
try:
    import nanodesign
except ImportError:
    import sys
    #base_path = '/Users/jonasfunke/NANODESIGN/nanodesign'
    #base_path = os.path.abspath( os.path.join( os.path.dirname(os.path.abspath( __file__)), '../nanodesign'))
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
    
def get_alpha_value(dna_structure, staple_indices, scaffold_rotation, scaffold_id, T_crit):
    # physical scaffold sequence
    scaffold_sequence = get_sequence(dna_structure.strands[scaffold_id]).replace('N', '')
    # physical scaffold length
    scaffold_length = len(scaffold_sequence)
    
    #loop through strands
    staple_domain_melt = []
    for strand in staple_indices:
        cur_domain_melt = []
        # loop through domain
        for domain in strand:
            # loop through bases in DOMAIN
            tmp = []
            for baseindex in domain:
                # physical index in scaffold
                i_physical = (baseindex+scaffold_rotation)%scaffold_length
                             
                #dna_structure.strands[scaffold_id].tour[i_physical+offset].seq
                tmp.append(scaffold_sequence[i_physical])
            if len(tmp)>1:
                # compute melting temperature of domain
                cur_domain_melt.append(MeltingTemp.Tm_NN(Seq(''.join(tmp))))
            else:
                cur_domain_melt.append(0.)    
            
        staple_domain_melt.append(max(cur_domain_melt))
    # calculate alpha value
    N_good = 0
    for T in staple_domain_melt:
        if T >= T_crit:
            N_good = N_good + 1
    return N_good / float(len(staple_domain_melt))
    
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
    #file_full_path_and_name = '/Users/jonasfunke/Dropbox/FRET_STAGE/Designs/FS-v7_spectrometer/FS-v7_order_2017-03-18/FS-v7.json'
    #seq_name = 'p8064'
    
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
     
    staple_indices = []
    for strand in dna_structure.strands:
        if not strand.is_scaffold:
            cur_strand = []
            for domain in strand.domain_list:
                cur_domain = []
                for base in domain.base_list:
                    if base.seq is not 'N':
                        # index of base on the physical scaffold
                        i = physical_index[dna_structure.strands[scaffold_id].get_base_index(base.across)]
                        cur_domain.append(i)
                cur_strand.append(cur_domain)
        staple_indices.append(cur_strand)
                    
                    
        
    #%% Loop through scaffold permutations and write sequencs of black oligos 
    #sequence of scaffold in design, this includes skips as 'N' 
    design_scaffold_sequence = get_sequence(dna_structure.strands[scaffold_id])
    #physical sequence of scaffold, this is the true scaffold sequence 
    physical_scaffold_length = len(design_scaffold_sequence)-design_scaffold_sequence.count('N')
    design_scaffold_length = len(design_scaffold_sequence)
    
    T_crit = 45.
    alpha_value = []
    # loop through scaffold permutations in physical space
    for i in range(0,physical_scaffold_length):
        if i%100 is 0: print('Computing permutation '+str(i) + ' of ' + str(physical_scaffold_length))
        
        # determine start of scaffold in design
        i_d = (-design_index[i])%design_scaffold_length

        # calculate alpha value
        alpha_value.append(get_alpha_value(dna_structure, staple_indices, i, scaffold_id, T_crit))
        #print('Rotation '+str(i) + ' alpha = '+str(round(alpha_value[-1],2)) )
        
       
        
        
        
    #%%
        
    file_out = output_path+'_permutations_black-sequences.csv'
    with open(file_out, 'wb') as csvfile:
        outputwriter = csv.writer(csvfile, delimiter=';')
        tmp = ['Permutation', 'Start helix', 'Start position', 'Alpha_value']
        outputwriter.writerow(tmp)
        
        # loop through scaffold permutations in physical space
        for i in range(0,physical_scaffold_length):
            
            # determine start of scaffold in design
            i_d = (-design_index[i])%design_scaffold_length


            # write data for this rotation
            tmp = [str(i), str(dna_structure.strands[scaffold_id].tour[i_d].h), str(dna_structure.strands[scaffold_id].tour[i_d].p), alpha_value[i]]            
            outputwriter.writerow(tmp)
        print('Output written to: ' + file_out)   
        
        
        
        
#%%

if __name__ == '__main__':
    main()