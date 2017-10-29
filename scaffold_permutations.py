#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 14:51:51 2017

@author: jonasfunke
    This script will compute various metrics/sequences for each scaffold permutation
"""

# Imports
from __future__ import print_function

import os
import csv
import argparse
import numpy
from Bio.SeqUtils import MeltingTemp # compute melting temperatures
from Bio.Seq import Seq #create biopython sequences
#from Bio.Alphabet import generic_dna


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
    

# find the domain with the highest melting temperature of each staple 
def get_max_domain_melt(dna_structure, staple_indices, scaffold_rotation, scaffold_id, print_staples):
    # physical scaffold sequence
    scaffold_sequence = get_sequence(dna_structure.strands[scaffold_id]).replace('N', '')
    # physical scaffold length
    scaffold_length = len(scaffold_sequence)
    #print(staple_indices)
    #loop through strands
    staple_domain_melt = []
    for strand in staple_indices:
        #cur_strand= []
        cur_domain_melt = []
        # loop through domain
        for domain in strand:
            # loop through bases in DOMAIN
            cur_domain = []
            for baseindex in domain:
                # physical index in scaffold
                i_physical = (baseindex+scaffold_rotation)%scaffold_length
                             
                #dna_structure.strands[scaffold_id].tour[i_physical+offset].seq
                cur_domain.append(scaffold_sequence[i_physical])
            if len(cur_domain)>1:
                # compute melting temperature of domain; reverse sequence of cur_domain, since it is on the scaffold and the indices follow staples
                cur_domain_melt.append(MeltingTemp.Tm_NN(Seq(''.join(cur_domain[::-1]))))
            else:
                cur_domain_melt.append(0.)
            #domain_seq_on_scaffold = Seq(''.join(cur_domain[::-1]), generic_dna)
            #cur_strand.append(str(domain_seq_on_scaffold.reverse_complement()))
           
        staple_domain_melt.append(max(cur_domain_melt))
        #if print_staples:
        #    print(str(cur_strand))
    return staple_domain_melt


# calulate alpha values for each tempereture given in T_crit
def get_alpha_values(staple_domain_melt, T_crit):
    
    alpha_values = []
    for T_cur in T_crit:
        N_good = 0
        for T in staple_domain_melt:
            if T >= T_cur:
                N_good = N_good + 1
        alpha_values.append( N_good / float(len(staple_domain_melt)) )
    
    
    #print('Number of oligos with a domain with T_m >= '+str(T_crit)+'°C: '+str(N_good)+' of '+str(len(staple_domain_melt)))
    return alpha_values
    
#%%
def main():
    #%%      
    # Parse arguments, TODO: use parser object
    parser = argparse.ArgumentParser(description='Do scaffold permutations.', prog='scaffold_permutations.py')
    parser.add_argument('path_to_json', type=str, nargs=1, help='Path to json file')
    parser.add_argument('scaffold', type=str, nargs=1, help='Scaffold: ...p7704, p8064', choices=['M13mp18', 'p7308', 'p7560', 'p7704', 'p8064', 'p8100', 'p8634', 'M13KO7'])
    parser.add_argument('--positions', type=str, nargs='+', help='Base to report on (scaffold base) as HelixID,Position. Example: --positions 0,100 27,212')
    #parser.add_argument('--alpha_value', action='store_const', const=True , help='Compute alpha value for each rotation')
    parser.add_argument('--alpha_value', type=str, nargs='?', help='Threshold temperatures, for which the alpha-value should be computed. Example: --alpha_value 40,45,55 ')
    parser.add_argument('--black_oligos', action='store_const', const=True , help='Compute sequences of black oligos for each rotation')
    
    #args = parser.parse_args(['/Users/jonasfunke/Documents/FRET_STAGE/test', 'p8064', '--positions', '2,3', '--alpha_value'])
    args = parser.parse_args()
        
    file_full_path_and_name = os.path.abspath( os.path.expanduser( args.path_to_json[0] ))
    seq_name = args.scaffold[0]

    
    base_coords = []
    if args.positions > 0:
        for pos in args.positions:
            base_coords.append([int(a) for a in pos.split(',')]) # [[H,p], [H,p], ...]
        
    #print(args.alpha_value)
    if len(args.alpha_value) > 0:
        #args.alpha_value = tmp;
        alpha_temperatures = [int(a) for a in args.alpha_value.split(',')] # [[H,p], [H,p], ...]
    #print(alpha_temperatures)
    #%%
    #file_full_path_and_name = '/Users/jonasfunke/Dropbox/FRET_STAGE/test/small.json'
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
    scaffold_id = -1
    for strand in dna_structure.strands:
        if strand.is_scaffold:
            print('Scaffold strand has index '+ str(strand.id))
            if scaffold_id is -1 :
                scaffold_id = strand.id
            else:
                print('WARNING: Multiple scaffolds detected')
    
    #%% create the maps, that map design to physical indeces
    design_index, physical_index = get_index_lists(dna_structure, scaffold_id)
     
    # get the indices of the staples on the scaffold strand
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
                        
    # get indices of black colored oligos
    reportcolor = 3355443 # this is the decimal value corresponding to 33333
    reportoligo_scaffold_indices =[]
    scaffold_bases_init =[]
    if args.black_oligos:
        for strand in dna_structure.strands:
            if strand.icolor==reportcolor:
                strand_scaffold_indeces = []
                for base in strand.tour:
                    if base.seq is not 'N':
                        strand_scaffold_indeces.append(physical_index[dna_structure.strands[scaffold_id].get_base_index(base.across)])          
                #strand_scaffold_indeces.reverse() # reverse the sequence to obtain 5' to 3'        
                reportoligo_scaffold_indices.append(strand_scaffold_indeces)
        print('Number of black oligos: '+str(len(reportoligo_scaffold_indices)))
        
    # get indeces on scaffold from candnano lattice positions     
    base_coords_on_scaffold = []
    for cur_coord in base_coords:
       # get base from helix and position coordinates
       cur_base = dna_structure.structure_helices_map[cur_coord[0]].scaffold_pos[cur_coord[1]]
       i = physical_index[dna_structure.strands[scaffold_id].get_base_index(cur_base)] 
       base_coords_on_scaffold.append(i)

    #%%          
    #sequence of scaffold in design, this includes skips as 'N' 
    design_scaffold_sequence = get_sequence(dna_structure.strands[scaffold_id])
    physical_scaffold_sequence = design_scaffold_sequence.replace('N', '')
    
    #physical sequence of scaffold, this is the true scaffold sequence 
    physical_scaffold_length = len(design_scaffold_sequence)-design_scaffold_sequence.count('N')
    design_scaffold_length = len(design_scaffold_sequence)
    
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


    #T_crit = 45.
    alpha_values = []
    base_list = []
    reportoligo_sequences = []
    mean_var = []
    # loop through scaffold permutations in physical space
    for i in range(0, physical_scaffold_length):
    
        if i%100 is 0: 
            #print('Computing permutation '+str(i) + ' of ' + str(physical_scaffold_length), end='')
        
            done = 100*float(i)/physical_scaffold_length
            # the exact output you're looking for:
            sys.stdout.write("[%-100s] %d%% done, iteration %i of %i" % ('='*int(done/1), done, i, physical_scaffold_length))
            sys.stdout.write('\r')            
            sys.stdout.flush()
        
        # calculate alpha value
        if len(alpha_temperatures)>0:
            staple_domain_max_melt = get_max_domain_melt(dna_structure, staple_indices, i, scaffold_id, 0)
            alpha_values.append(get_alpha_values(staple_domain_max_melt, alpha_temperatures))
            #if i==0: print("Current alpha value is " + str(alpha_value[-1]))
            mean_var.append([numpy.mean(staple_domain_max_melt), numpy.var(staple_domain_max_melt)])
        
        # get the bases at the current rotation
        cur_base_list = []
        for cur_base_index in base_coords_on_scaffold:
            cur_base_list.append(physical_scaffold_sequence[(cur_base_index+i)%physical_scaffold_length])
    
        base_list.append(cur_base_list)
        
        # black oligos
        if args.black_oligos:
            tmp = []
            for strand in reportoligo_scaffold_indices:
                cur_strand = []
                for baseindex in strand:
                    cur_strand.append(complement[physical_scaffold_sequence[(baseindex+i)%physical_scaffold_length]])
                tmp.append(''.join(cur_strand))
            reportoligo_sequences.append(tmp)

    sys.stdout.write('\n')            

    #print(alpha_values)
    #%%
    #if args.alpha_value:
    #    print('Maximal alpha values:')
    #    i_sorted = sorted(range(len(alpha_value)), key=lambda k: alpha_value[k])
    #    for j in range(10):
    #        i_d = design_index[(-i_sorted[-1-j])%physical_scaffold_length]
    #        print('Alpha value ' + str(alpha_value[i_sorted[-1-j]]) + ' on Helix H' + str(dna_structure.strands[scaffold_id].tour[i_d].h) + ' at position ' +  str(dna_structure.strands[scaffold_id].tour[i_d].p) + '  (permutation ' + str(i_sorted[-1-j]) + ')' )
              
    #print(mean_var)

    #%%
    file_out = output_path+'_scaffold-permutations.csv'
    with open(file_out, 'wb') as csvfile:
        outputwriter = csv.writer(csvfile, delimiter=';')
        tmp = ['Permutation', 'Scaffold start helix', 'Scaffold start position']
        
        for i in range(0,len(alpha_temperatures)):
            tmp.append('Alpha_value_' + str(alpha_temperatures[i]))
            
        if len(alpha_temperatures)>0:
            tmp.append('Mean of largest domain Tm')
            tmp.append('Variance of largest domain Tm')
            
        for cur_coord in base_coords:
            tmp.append('H' + str(cur_coord[0]) + ',' + str(cur_coord[1]))
        
        for i in range(0,len(reportoligo_scaffold_indices)):
            tmp.append('Black oligo '+ str(i))
            
        for i in range(0,len(reportoligo_scaffold_indices)):
            tmp.append('Black oligo '+ str(i)+' T_m')
        
        outputwriter.writerow(tmp)
        
        # loop through scaffold permutations in physical space
        for i in range(0,physical_scaffold_length):
            
            # determine start of scaffold in design
            i_d = design_index[(-i)%physical_scaffold_length]

            # write data for this rotation
            tmp = [str(i), str(dna_structure.strands[scaffold_id].tour[i_d].h), str(dna_structure.strands[scaffold_id].tour[i_d].p)]            
            for j in range(0, len(alpha_temperatures)):
                tmp.append(alpha_values[i][j])
            
            if len(alpha_temperatures)>0:
                tmp.append(mean_var[i][0])
                tmp.append(mean_var[i][1])
            
            for cur_base in base_list[i]:
                tmp.append(cur_base)
                
            if args.black_oligos:
                for strand in reportoligo_sequences[i]:
                    tmp.append(strand)
                    
            if args.black_oligos:
                for strand in reportoligo_sequences[i]:
                    tmp.append(round(MeltingTemp.Tm_NN(Seq(''.join(strand))),1))
                    
            outputwriter.writerow(tmp)
        print('Output written to: ' + file_out)   

    if args.alpha_value:
        print('REMARK: alpha values of the modified design might differ slightly form the predicted alpha values. Altered staple segmentation caused by shifting the scaffold starting point can slightly change the final alpha value.')        
     
        
        
#%%

if __name__ == '__main__':
    main()
