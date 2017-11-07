#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 12:41:12 2017

@author: jonasfunke
"""


# Imports
from __future__ import print_function

import os
#import csv
import argparse
import numpy
#from Bio.SeqUtils import MeltingTemp # compute melting temperatures
#from Bio.Seq import Seq #create biopython sequences
#from Bio.Alphabet import generic_dna
import matplotlib.pyplot as plt #plotting 


# try to impot nanodesign package. I assume the script is in x/somename/design_statistics.py and the nanodesign package is in x/nanodesign
try:
    import nanodesign
except ImportError:
    import sys
    base_path = '/Users/jonasfunke/NANODESIGN/nanodesign'
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

def get_sequence(strand):
    cur_seq = []
    for base in strand.tour:
        cur_seq.append(base.seq)
    return ''.join(cur_seq)

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
    
def get_staple_loop_lengths(path_to_file, scaffold_name):
    # parse filename and create output directory
    file_name = os.path.basename( path_to_file )
    output_path = os.path.dirname(path_to_file) + os.sep + file_name[:-5] + '_scaffold-permutations' + os.sep
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = output_path+file_name[:-5]
    
    # Read cadnano file and create dna structure.
    converter = read_file(path_to_file, scaffold_name)
    dna_structure = converter.dna_structure     
    dna_structure.compute_aux_data() # compute domain data
        
    #determine scaffold strand id
    scaffold_id = -1
    for strand in dna_structure.strands:
        if strand.is_scaffold:
            print('Scaffold strand has index '+ str(strand.id))
            if scaffold_id is -1 :
                scaffold_id = strand.id
            else:
                print('WARNING: Multiple scaffolds detected')
    
    # create the maps, that map design to physical indeces
    design_index, physical_index = get_index_lists(dna_structure, scaffold_id)

    #sequence of scaffold in design, this includes skips as 'N' 
    design_scaffold_sequence = get_sequence(dna_structure.strands[scaffold_id])
    physical_scaffold_length = len(design_scaffold_sequence)-design_scaffold_sequence.count('N')

    # get the indices of the staples on the scaffold strand and scaffold loop length
    staple_indices = []
    loop_lengths = []
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
            print(cur_strand)
            
            cur_strand_loop_length = []
            for i in range(1,len(cur_strand)):
                if len(cur_strand[i])>0 and len(cur_strand[i-1])>0:
                    d = abs(cur_strand[i][0]-cur_strand[i-1][-1])
                    cur_strand_loop_length.append(min(d, physical_scaffold_length-d))
                
            print(cur_strand_loop_length)
            loop_lengths.append(cur_strand_loop_length)
                
    return loop_lengths

#%%
def main():
    #%%      
    # Parse arguments, TODO: use parser object
    #parser = argparse.ArgumentParser(description='Do scaffold permutations.', prog='scaffold_permutations.py')
    #parser.add_argument('path_to_json', type=str, nargs=1, help='Path to json file')
    #parser.add_argument('scaffold', type=str, nargs=1, help='Scaffold: ...p7704, p8064', choices=['M13mp18', 'p7308', 'p7560', 'p7704', 'p8064', 'p8100', 'p8634', 'M13KO7'])
    #parser.add_argument('--positions', type=str, nargs='+', help='Base to report on (scaffold base) as HelixID,Position. Example: --positions 0,100 27,212')
    #parser.add_argument('--alpha_value', action='store_const', const=True , help='Compute alpha value for each rotation')
    #parser.add_argument('--alpha_value', type=str, nargs='?', const='45', help='Threshold temperatures, for which the alpha-value should be computed. Example: --alpha_value 40,45,55 ')
    #parser.add_argument('--black_oligos', action='store_const', const=True , help='Compute sequences of black oligos for each rotation')
    
    #args = parser.parse_args()
        
    #file_full_path_and_name = os.path.abspath( os.path.expanduser( args.path_to_json[0] ))
    #seq_name = args.scaffold[0]

 
        
    #print(alpha_temperatures)
    #%%
    path_to_file = '/Users/jonasfunke/Dropbox/FRET_STAGE/test/FS-v6.json'
    seq_name = 'p8064'
   
    loop_lengths = get_staple_loop_lengths(path_to_file, seq_name)
    
    ll_v6 = []
    for i in range(len(loop_lengths)):
        for j in range(len(loop_lengths[i])):
            ll_v6.append(loop_lengths[i][j])
        
        
    path_to_file = '/Users/jonasfunke/Dropbox/FRET_STAGE/test/FS-v7.json'
    seq_name = 'p8064'
   
    loop_lengths = get_staple_loop_lengths(path_to_file, seq_name)
    
    ll_v7 = []
    for i in range(len(loop_lengths)):
        for j in range(len(loop_lengths[i])):
            ll_v7.append(loop_lengths[i][j])
            
    path_to_file = '/Users/jonasfunke/Dropbox/FRET_STAGE/test/FS-v8.json'
    seq_name = 'p8064'
   
    loop_lengths = get_staple_loop_lengths(path_to_file, seq_name)
    
    ll_v8 = []
    for i in range(len(loop_lengths)):
        for j in range(len(loop_lengths[i])):
            ll_v8.append(loop_lengths[i][j])
            
    
    path_to_file = '/Users/jonasfunke/Dropbox/FRET_STAGE/test/FS-v9.json'
    seq_name = 'p8064'
   
    loop_lengths = get_staple_loop_lengths(path_to_file, seq_name)
    
    ll_v9 = []
    for i in range(len(loop_lengths)):
        for j in range(len(loop_lengths[i])):
            ll_v9.append(loop_lengths[i][j])        
    #%%
    output_path = '/Users/jonasfunke/Dropbox/FRET_STAGE/test/'
    bin_width = 500.
    x_min = 0
    x_max = 5000
    bins= numpy.arange(x_min-bin_width/2, x_max+bin_width/2+bin_width, bin_width)
    
    fig = plt.figure()
    plt.hist(ll_v6, bins, alpha=0.5,  histtype='bar', label='v6') 
    plt.hist(ll_v7, bins, alpha=0.5, histtype='bar', label='v7') 
    plt.hist(ll_v8, bins, alpha=0.5, histtype='bar', label='v8') 
    plt.hist(ll_v9, bins, alpha=0.5, histtype='bar', label='v9') 
    #plt.title("Average  " + name + " = " +str(round(numpy.mean(data),1)) + unit)
    plt.xlabel('loop length [bases]')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    #plt.xticks(numpy.arange(x_min, x_max+1, bin_width))
    fig.savefig(output_path+'ScaffoldLoopLengthDistribution_log.pdf')
    plt.show()
    #plt.close() 

    
    bin_width = 0.5
    x_min = 1
    x_max = numpy.log(5000)
    bins= numpy.arange(x_min-bin_width/2, x_max+bin_width/2+bin_width, bin_width)
    
    fig = plt.figure()
    plt.hist(numpy.log(ll_v6), bins, alpha=0.5,  histtype='step', label='v6') 
    plt.hist(numpy.log(ll_v7), bins, alpha=0.5, histtype='step', label='v7') 
    plt.hist(numpy.log(ll_v8), bins, alpha=0.5, histtype='step', label='v8') 
    plt.hist(numpy.log(ll_v9), bins, alpha=0.5, histtype='step', label='v9') 
    #plt.title("Average  " + name + " = " +str(round(numpy.mean(data),1)) + unit)
    plt.xlabel('log(loop length)')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    #plt.xticks(numpy.arange(x_min, x_max+1, bin_width))
    fig.savefig(output_path+'ScaffoldLoopLengthDistribution_log.pdf')
    plt.show()
    
    print(numpy.mean(numpy.log(ll_v6)))
    print(numpy.mean(numpy.log(ll_v7)))
    print(numpy.mean(numpy.log(ll_v8)))
    print(numpy.mean(numpy.log(ll_v9)))        
#%%
    path_to_file = '/Users/jonasfunke/Dropbox/FRET_STAGE/test/T1_v1.2.json'
    seq_name = 'p8064'
   
    loop_lengths = get_staple_loop_lengths(path_to_file, seq_name)
    
    ll_v1 = []
    for i in range(len(loop_lengths)):
        for j in range(len(loop_lengths[i])):
            ll_v1.append(loop_lengths[i][j])
            
    path_to_file = '/Users/jonasfunke/Dropbox/FRET_STAGE/test/T1_v2.json'
    seq_name = 'p8064'
   
    loop_lengths = get_staple_loop_lengths(path_to_file, seq_name)
    
    ll_v2 = []
    for i in range(len(loop_lengths)):
        for j in range(len(loop_lengths[i])):
            ll_v2.append(loop_lengths[i][j])
            
    
    path_to_file = '/Users/jonasfunke/Dropbox/FRET_STAGE/test/T3D_v1.json'
    seq_name = 'p8064'
   
    loop_lengths = get_staple_loop_lengths(path_to_file, seq_name)
    
    ll_v3 = []
    for i in range(len(loop_lengths)):
        for j in range(len(loop_lengths[i])):
            ll_v3.append(loop_lengths[i][j])        
    #%%
    output_path = '/Users/jonasfunke/Dropbox/FRET_STAGE/test/'
    bin_width = 500.
    x_min = 0
    x_max = 5000
    bins= numpy.arange(x_min-bin_width/2, x_max+bin_width/2+bin_width, bin_width)
    
    fig = plt.figure()
    plt.hist(ll_v1, bins, alpha=0.5,  histtype='bar', label='T1_v1.2') 
    plt.hist(ll_v2, bins, alpha=0.5, histtype='bar', label='T1_v2') 
    plt.hist(ll_v3, bins, alpha=0.5, histtype='bar', label='T3D') 
    #plt.title("Average  " + name + " = " +str(round(numpy.mean(data),1)) + unit)
    plt.xlabel('loop length [bases]')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    #plt.xticks(numpy.arange(x_min, x_max+1, bin_width))
    fig.savefig(output_path+'T_ScaffoldLoopLengthDistribution.pdf')
    plt.show()
    #plt.close() 

    
    bin_width = 0.1
    x_min = 1
    x_max = numpy.log(5000)
    bins= numpy.arange(x_min-bin_width/2, x_max+bin_width/2+bin_width, bin_width)
    
    fig = plt.figure()
    plt.hist(numpy.log(ll_v1), bins, alpha=0.5,  histtype='step', label='T1_v1.2') 
    plt.hist(numpy.log(ll_v2), bins, alpha=0.5, histtype='step', label='T1_v2') 
    plt.hist(numpy.log(ll_v3), bins,  alpha=0.5, histtype='step', label='T3D') 
    #plt.title("Average  " + name + " = " +str(round(numpy.mean(data),1)) + unit)
    plt.xlabel('log(loop length)')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    #plt.xticks(numpy.arange(x_min, x_max+1, bin_width))
    fig.savefig(output_path+'T_ScaffoldLoopLengthDistribution_log.pdf')
    plt.show()
    
    
    #%%
    bin_width = 0.1
    x_min = 1
    x_max = numpy.log(5000)
    bins= numpy.arange(x_min-bin_width/2, x_max+bin_width/2+bin_width, bin_width)
    
    fig = plt.figure()
    plt.hist(numpy.log(ll_v1), bins, alpha=0.5,  histtype='bar', label='T1_v1.2') 
    #plt.hist(numpy.log(ll_v2), bins, alpha=0.5, histtype='bar', label='T1_v2') 
    #plt.hist(numpy.log(ll_v3), bins,  alpha=0.5, histtype='bar', label='T3D') 
    
    plt.xlabel('log(loop length)')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    #plt.xticks(numpy.arange(x_min, x_max+1, bin_width))
    fig.savefig(output_path+'T_ScaffoldLoopLengthDistribution_log2.pdf')
    plt.show()
    
    #%%
    print(numpy.mean(numpy.log(ll_v6)))
    print(numpy.mean(numpy.log(ll_v7)))
    print(numpy.mean(numpy.log(ll_v8)))
    print(numpy.mean(numpy.log(ll_v9)))  
    
    #%%

if __name__ == '__main__':
    main()
