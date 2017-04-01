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
from matplotlib import cm #colormaps for coloring staples
from Bio.SeqUtils import MeltingTemp # compute melting temperatures
from Bio.Seq import Seq #create biopython sequences

# try to impot nanodesign package. I assume the script is in x/somename/design_statistics.py and the nanodesign package is in x/nanodesign
try:
    import nanodesign
except ImportError:
    import sys
    #base_path = '/Users/jonasfunke/NANODESIGN/nanodesign'
    base_path = os.path.abspath( os.path.join( os.path.dirname(os.path.abspath( __file__)), '../nanodesign/'))
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

def write_colored_json(data, dna_structure, colormap, output_file ): # write cadnano file colorcoded using a 2D list [[index, value], [index+1, value],...]
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
    #seq_name = 'p7308'
    
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
    for strand in dna_structure.strands:
        if not strand.is_scaffold:
            staple_length.append([strand.id, get_strand_length(strand)])
            cur_domain_len = []
            cur_domain_mt = []
            for domain in strand.domain_list:
                cur_domain_len.append(len(domain.sequence.replace('N', '')))
                if len(domain.sequence.replace('N', ''))>0:
                    cur_domain_mt.append(MeltingTemp.Tm_NN(Seq(domain.sequence.replace('N', ''))))
            domain_max_melt.append([strand.id, max(cur_domain_mt)])
            avg_domain_length.append([strand.id, numpy.mean(cur_domain_len)])
            domain_max_length.append([strand.id, max(cur_domain_len)])
            number_of_domains.append([strand.id, len(cur_domain_mt)]) # this will not incude polyT overhangs as separate domains

    # compute histograms
    my_histogram(staple_length, output_path, 'Length of staples', 'bases')
    my_histogram(avg_domain_length, output_path, 'Average staple-domain length', 'bp')
    my_histogram(number_of_domains, output_path, 'Number of domains per staple', '')
    my_histogram(domain_max_melt, output_path, 'Melting temperature of domain with highest melting temperature', 'deg')
    my_histogram(domain_max_length, output_path, 'Length of longest domain', '')
        
    # wirte json files
    write_colored_json(staple_length, dna_structure, cm.afmhot, output_path+'_colorcoded_StapleLength.json' )
    write_colored_json(domain_max_melt, dna_structure, cm.afmhot, output_path+'_colorcoded_LargestDomainMelt.json' )
    write_colored_json(domain_max_length, dna_structure, cm.afmhot, output_path+'_colorcoded_LongestDomainMelt.json' )
    write_colored_json(number_of_domains, dna_structure, cm.afmhot, output_path+'_colorcoded_NumberOfDomains.json' )
    
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
    print('Alpha value is '+str(N_good / float(len(domain_max_melt))) )
                          
#%%

if __name__ == '__main__':
    main()