#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:42:58 2017

@author: jonasfunke
    This scripted computes strand statistics for a cadnano design. It uses the nanodesign package from autodesk
    usage: python design_statistics.py path/to/json_file.json p7560
"""

# compatibility to python 3.5
#from __future__ import (absolute_import, division, print_function, unicode_literals)
import __future__
import os
import numpy #for writing matrices
import matplotlib.pyplot as plt #plotting 
import argparse # parsing arguments
from matplotlib import cm #colormaps for coloring staples
from Bio.SeqUtils import MeltingTemp # compute melting temperatures
from Bio.Seq import Seq #create biopython sequences
import csv

# try to impot nanodesign package. I assume the script is in x/somename/design_statistics.py and the nanodesign package is in x/nanodesign
try:
    import nanodesign
except ImportError:
    import sys
    base_path = '/Users/jonasfunke/NANODESIGN/nanodesign'
    #base_path = os.path.abspath( os.path.join( os.path.dirname(os.path.abspath( __file__)), '../nanodesign/'))
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

def get_strand_length(strand): # computes total length of strand in bases
    length = 0
    for base in strand.tour:
        length = length + 1 + base.num_deletions
    return length

def bases_bound_to_scaffold(strand): # computes number of bases that are bound to the scaffold
    sequence= []
    for base in strand.tour:
        sequence.append(base.seq)
    tmp = ''.join(sequence)
    return len(tmp.translate(None, 'N'))



def my_histogram(data_in, output_path, name, unit):
    data=[data_in[i][1] for i in range(len(data_in))]
    bin_width = 1.
    fig = plt.figure()
    x_min = numpy.ceil(min(data))
    x_max = numpy.floor(max(data))
    plt.hist(data, bins= numpy.arange(x_min-bin_width/2, x_max+bin_width/2+bin_width, bin_width), histtype='bar') 
    plt.title("Average  " + name + " = " +str(round(numpy.mean(data),1)) + unit)
    plt.xlabel(name + ' [' + unit + ']')
    plt.ylabel('Frequency')
    #plt.xticks(numpy.arange(x_min, x_max+1, bin_width))
    fig.savefig(output_path+'_stats_' + name.replace(' ','_')+ '.pdf')
    #plt.show()
    plt.close() 


    
    

def write_colored_json(data, dna_structure, colormap, output_file, limits): # write cadnano file colorcoded using a 2D list [[index, value], [index+1, value],...]
    if len(limits)>0: 
        x_min = limits[0]
        x_max = limits[1]
    else: # if no limits are given, use max and min
        x_min = min([data[i][1] for i in range(len(data))])
        x_max = max([data[i][1] for i in range(len(data))])
    dna_structure_out = dna_structure
    for i in range(len(data)):
        x = data[i][1]
        colorindex = int(254*(x-x_min)/(x_max-x_min))+1 # asumin 0-255 colors in colormap
        dna_structure_out.strands[data[i][0]].color = list(colormap(colorindex)[0:-1])                
       
    con = Converter()
    con.dna_structure = dna_structure_out
    con.write_cadnano_file(output_file)

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

#%%
def main():
    parser = argparse.ArgumentParser(description='Analyze cadnano design.', prog='design_statistics.py')
    parser.add_argument('path_to_json', type=str, nargs=1, help='Path to json file')
    parser.add_argument('scaffold', type=str, nargs=1, help='Scaffold: ...p7704, p8064', choices=['M13mp18', 'p7308', 'p7560', 'p7704', 'p8064', 'p8100', 'p8634', 'M13KO7'])
    parser.add_argument('--alpha_value', type=str, nargs='?', help='Threshold temperatures, for which the alpha-value should be computed. Example: --alpha_value 40,45,55 ')
    parser.add_argument('--T_crit', type=float, nargs='?', default=45.0 , help='Melting temperture to color json. If not used, 45deg will be used. Example: --T_crit 35')
    args = parser.parse_args()
    
    T_crit = args.T_crit
    file_full_path_and_name = os.path.abspath( os.path.expanduser( args.path_to_json[0] )) #os.path.abspath( os.path.expanduser( sys.argv[1] ))
    seq_name = args.scaffold[0] #sys.argv[2]

    #%%
    file_full_path_and_name = '/Users/jonasfunke/Desktop/LGV3_4st_9_final.json'
    seq_name = 'p8064'
    T_crit = 45.
    
    # parse filename and create output directory
    file_name = os.path.basename( file_full_path_and_name )
    output_path = os.path.dirname(file_full_path_and_name) + os.sep + file_name[:-5] + '_analysis' + os.sep
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = output_path+file_name[:-5]
    
    
    # Read cadnano file and create dna structure.
    converter = read_file(file_full_path_and_name, seq_name)
    dna_structure = converter.dna_structure     
    dna_structure.compute_aux_data() # compute domain data
    
    #%%        
    # compute staple statistics
    staple_length=[] # length of staple
    avg_domain_length = [] # average domain length of staple
    domain_max_melt = [] # largest domain melting-temperature of staple
    domain_max_length = [] #longest domain of staple
    number_of_domains = [] #number of domains of staple
    number_of_long_domains = [] #length of domains
    strand_melting_temp = []
    larger_than_T_crit = [] # bool whether max domain melt is larger than T_crit
    for strand in dna_structure.strands:
        if not strand.is_scaffold:
            cur_strand = []
            staple_length.append([strand.id, get_strand_length(strand)])
            cur_domain_len = []
            cur_domain_mt = []
            n_long = 0
            for domain in strand.domain_list:
                cur_domain_len.append(len(domain.sequence.replace('N', '')))
                cur_strand.append(domain.sequence.replace('N', ''))
                if len(domain.sequence.replace('N', ''))>0:
                    cur_domain_mt.append(MeltingTemp.Tm_NN(Seq(domain.sequence.replace('N', ''))))
                # check if staple has more than one long domain (>10 bp) (this can be optimized / made nicer)
                if len(domain.sequence.replace('N', '')) > 10:
                    n_long = n_long + 1
            domain_max_melt.append([strand.id, max(cur_domain_mt)])
            avg_domain_length.append([strand.id, numpy.mean(cur_domain_len)])
            domain_max_length.append([strand.id, max(cur_domain_len)])
            number_of_domains.append([strand.id, len(cur_domain_mt)]) # this will not incude polyT overhangs as separate domains
            number_of_long_domains.append([strand.id, n_long])        
            strand_melting_temp.append([strand.id, MeltingTemp.Tm_NN(Seq(''.join(cur_strand)))])
            larger_than_T_crit.append([strand.id, int(max(cur_domain_mt)>=T_crit)])
            
            #if max(cur_domain_mt) >= 45. :
            #print(cur_strand)

            
    
    
    #%%
        
    # compute histograms
    my_histogram(staple_length, output_path, 'Length of staples', 'bases')
    my_histogram(avg_domain_length, output_path, 'Average staple-domain length', 'bp')
    my_histogram(number_of_domains, output_path, 'Number of domains per staple', '')
    my_histogram(domain_max_melt, output_path, 'Melting temperature of domain with highest melting temperature', 'deg')
    my_histogram(domain_max_length, output_path, 'Length of longest domain', '')
    my_histogram(strand_melting_temp, output_path, 'Strand melting temperature', 'deg')
        
    # wirte json files
    write_colored_json(staple_length, dna_structure, cm.afmhot, output_path+'_colorcoded_StapleLength.json', [])
    write_colored_json(domain_max_melt, dna_structure, cm.afmhot, output_path+'_colorcoded_MeltTempOfLargestDomainMelt_0-60.json', [0, 60] )
    write_colored_json(domain_max_length, dna_structure, cm.afmhot, output_path+'_colorcoded_LengthOfLongestDomain.json', [] )
    write_colored_json(number_of_domains, dna_structure, cm.afmhot, output_path+'_colorcoded_NumberOfDomains.json', [] )
    write_colored_json(number_of_long_domains, dna_structure, cm.bwr, output_path+'_colorcoded_NumberOfLongDomains.json',[])
    write_colored_json(number_of_long_domains, dna_structure, cm.coolwarm, output_path+'_colorcoded_NumberOfLongDomains_binary.json', [0, 2])
    write_colored_json(strand_melting_temp, dna_structure, cm.afmhot, output_path+'_colorcoded_StrandMeltingTemperature_0-60.json', [0, 60])
    write_colored_json(larger_than_T_crit, dna_structure, cm.bwr, output_path+'_colorcoded_LargerEqualThan'+str(T_crit)+'deg.json',[])
    # compute domain length statitics
    domain_lengths = []
    for domain in dna_structure.domain_list:
        if len(domain.sequence.replace('N', ''))> 0:
            domain_lengths.append([domain.id, len(domain.sequence.replace('N', ''))]) #length of domain, removed N for skipped bases
    
    my_histogram(domain_lengths, output_path, 'Length of domain', 'bp')
    

    # compute alpha value
    T_crit = 45.
    N_good = 0
    for T in domain_max_melt:
        #print(T[1])
        if T[1] >= T_crit:
            N_good = N_good + 1
    print('Number of oligos with a domain with T_m >= '+str(T_crit)+'°C: '+str(N_good)+' of '+str(len(domain_max_melt)))
    
    #determine scaffold strand id
    for strand in dna_structure.strands:
        if strand.is_scaffold:
            scaffold_id = strand.id
    print('Alpha value is '+str(N_good / float(len(domain_max_melt))) + ' on H' + str(dna_structure.strands[scaffold_id].tour[0].h) + ' at position '  + str(dna_structure.strands[scaffold_id].tour[0].p) )
    
    n_long_tot = 0
    for s in number_of_long_domains:
        if s[1]>1:
            n_long_tot = n_long_tot + 1
    print('Number of oligos with more than one long (>10 bp) domains: ' + str(n_long_tot))
                 

    #%% scaffold loop lengths distribution
    #determine scaffold strand id
    scaffold_id = -1
    for strand in dna_structure.strands:
        if strand.is_scaffold:
            print('Scaffold strand has index '+ str(strand.id))
            if scaffold_id is -1 :
                scaffold_id = strand.id
            else:
                print('WARNING: Multiple scaffolds detected')

    # create maps, that map design to physical scaffold indeces
    design_index, physical_index = get_index_lists(dna_structure, scaffold_id)

    #sequence of scaffold in design, this includes skips as 'N' 
    design_scaffold_sequence = get_sequence(dna_structure.strands[scaffold_id])
    physical_scaffold_length = len(design_scaffold_sequence)-design_scaffold_sequence.count('N')
    
    #%%
    # get the indices of the staples on the scaffold strand and scaffold loop length

    f= open("/Users/jonasfunke/Desktop/test.bild","w+")
    staple_indices = []
    #staple_base_coords = []
    loop_lengths = []
    for strand in dna_structure.strands:
        if not strand.is_scaffold:  
            cur_strand = []
            cur_strand_base_coord = []
            for domain in strand.domain_list:
                cur_domain = []
                cur_domain_base_coord = []
                for base in domain.base_list:
                    if base.seq is not 'N':
                        # index of base on the physical scaffold
                        i = physical_index[dna_structure.strands[scaffold_id].get_base_index(base.across)]
                        cur_domain.append(i)
                        cur_domain_base_coord.append(base.nt_coords)
                        #f.write(".color 1 1 1\n" )
                        #f.write(".sphere %f %f %f 0.1\n" % (base.nt_coords[0],base.nt_coords[1],base.nt_coords[2])) 
                        
                for i in range(1,len(domain.base_list)):
                    if (domain.base_list[i-1].seq is not 'N') and (domain.base_list[i].seq is not 'N'):
                        cur_base = domain.base_list[i]
                        prev_base = domain.base_list[i-1]
                        f.write(".color 1 1 1\n" )
                        f.write(".v %f %f %f %f %f %f 0.1\n" % (prev_base.nt_coords[0],prev_base.nt_coords[1],prev_base.nt_coords[2], cur_base.nt_coords[0],cur_base.nt_coords[1],cur_base.nt_coords[2])) 
                cur_strand.append(cur_domain)
                cur_strand_base_coord.append(cur_domain_base_coord)
            staple_indices.append(cur_strand)
            #staple_base_coords.append(cur_strand_base_coord)
            cur_strand_loop_length = []
            for i in range(1,len(cur_strand)):
                if len(cur_strand[i])>0 and len(cur_strand[i-1])>0:
                    d = abs(cur_strand[i][0]-cur_strand[i-1][-1])
                    cur_strand_loop_length.append(min(d, physical_scaffold_length-d))
                    if min(d, physical_scaffold_length-d) > 10:
                        #draw connection
                        x = min(d, physical_scaffold_length-d)/float(physical_scaffold_length/2)
                        f.write(".color %f %f 0\n" % (2*x, 2*(1-x)) )
                        f.write(".cylinder %f %f %f %f %f %f 0.1\n" % (cur_strand_base_coord[i][0][0], cur_strand_base_coord[i][0][1], cur_strand_base_coord[i][0][2], cur_strand_base_coord[i-1][-1][0], cur_strand_base_coord[i-1][-1][1], cur_strand_base_coord[i-1][-1][2]))
                    loop_lengths.append(cur_strand_loop_length)
                    
        
    f.close()
    #%%
    # make historgram
    ll_tmp = []
    for i in range(len(loop_lengths)):
        for j in range(len(loop_lengths[i])):
            ll_tmp.append(loop_lengths[i][j])        
    
    fig = plt.figure(figsize = (12,6))
    
    plt.subplot(121)
    bin_width = 500.
    x_min = 0
    x_max = 5000
    bins= numpy.arange(x_min-bin_width/2, x_max+bin_width/2+bin_width, bin_width)
    plt.hist(ll_tmp, bins,  histtype='bar', label='loop length') 
    plt.title("Average loop length = " +str(round(numpy.mean(ll_tmp),1)) + ' bases')
    plt.xlabel('Loop length [bases]')
    plt.ylabel('Frequency')
    
    plt.subplot(122)
    bin_width = 0.5
    x_min = 1
    x_max = numpy.log(5000)
    bins= numpy.arange(x_min-bin_width/2, x_max+bin_width/2+bin_width, bin_width)
    plt.hist(numpy.log(ll_tmp), bins,  histtype='bar') 
    plt.title("Average log(loop length) = " +str(round(numpy.mean(numpy.log(ll_tmp)),1)) )
    plt.xlabel('Log(Loop length [bases])')
    plt.ylabel('Frequency')
    
    fig.savefig(output_path+'_stats_ScaffoldLoopLengthDistribution.pdf')
    #plt.show()
    plt.close()                
    print('Average loop length is ' + str(numpy.mean(ll_tmp)) + ' bases.')
    print('Average log(loop length) is ' + str(numpy.mean(numpy.log(ll_tmp))) + ' bases.')

    #%% Print max domain Tm of oligos
    file_out = output_path+'_MaxDomainMeltingTemp.txt'
    with open(file_out, 'wb') as csvfile:
        outputwriter = csv.writer(csvfile, delimiter='\t')
        
        for i in range(len(domain_max_melt)):
            outputwriter.writerow([domain_max_melt[i][1]])
            #print(max(loop_lengths[i]))

    file_out = output_path+'_LoopLengths.txt'
    with open(file_out, 'wb') as csvfile:
        outputwriter = csv.writer(csvfile, delimiter='\t')
        
        for i in range(len(loop_lengths)):
            for j in range(len(loop_lengths[i])):
                outputwriter.writerow([loop_lengths[i][j]])

    print('Output written to: ' + file_out)  
    
    
    
    
    #%%

if __name__ == '__main__':
    main()