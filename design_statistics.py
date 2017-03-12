#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:42:58 2017

@author: jonasfunke
"""


#%% Imports
import os
import numpy #for writing matrices
import matplotlib.pyplot as plt #plotting 
from matplotlib import cm #colormaps for coloring staples
from Bio.SeqUtils import MeltingTemp # compute melting temperatures
from Bio.Seq import Seq #create biopython sequences

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

#%% Functions


def read_file(file_name, seq_name):
    """ Read in a cadnano file. """
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


def main():
    #%% Parse arguments
    if (len(sys.argv) != 3):
        sys.stderr.write("**** ERROR: Wrong number of arguments.\n") 
        sys.stderr.write("Usage: design_statistics.py <filename> scaffold_sequence\n")
        sys.stderr.write("Output: If <filename> is path/name.json, output will be placed in path/name_analysis/\n")
        sys.exit(1)
    #%% Parse arguments
    file_full_path_and_name = os.path.abspath( os.path.expanduser( sys.argv[1] ))
    #file_full_path_and_name = '/Users/jonasfunke/Dropbox/FRET_STAGE/Designs/FS-v6_spectrometer/twist_screen/FS-v6_019_deacivated.json'
    #seq_name = 'p7704'
    
    file_name = os.path.basename( file_full_path_and_name )

    output_path = os.path.dirname(file_full_path_and_name) + os.sep + file_name[:-5] + '_analysis' + os.sep
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = output_path+file_name[:-5]
    
    seq_name = sys.argv[2]
    
    # Read cadnano file and create dna structure.
    converter = read_file(file_full_path_and_name, seq_name)
    dna_structure = converter.dna_structure
    #%%
    tmp = 0
    for strand in dna_structure.strands:
        if strand.is_scaffold:
            scaffold_id = strand.id
            print 'Length of scaffold strand ' + str(strand.id) + ': ' + str(get_strand_length(strand))
            tmp = tmp+1
    print 'Number of scaffolds: ' + str(tmp)        
            
    #%% compute staple length statistics
    staple_length=[]
    for strand in dna_structure.strands:
        if not strand.is_scaffold:
            staple_length.append(get_strand_length(strand))
            #print 'Length: '+str(get_physical_length(strand))

    # plot histogram of staple length
    l_min = min(staple_length)
    l_max = max(staple_length)
    dl = 1.0 # staple length increment
    
    fig = plt.figure()
    plt.hist(staple_length, bins=numpy.linspace(l_min-dl/2, l_max+dl/2, l_max-l_min+2)) #, bins=numpy.linspace(l_min-0.5, l_max+0.5, l_max-l_min+2))  # plt.hist passes it's arguments to np.histogram
    plt.title("Average staple length = "+str(round(numpy.mean(staple_length),1)) + " bases")
    plt.xlabel('Length of staple [bp]')
    plt.ylabel('Frequency')
    fig.savefig(output_path+'_stats_StapleLength.pdf')
    #plt.show()
    plt.close()        
        
    # write json which is colorcoded according to staple length        
    dna_structure_length_colored = dna_structure
    for strand in dna_structure_length_colored.strands:
        if not strand.is_scaffold:
            l = get_strand_length(strand)
            colorindex = int(254*(l-l_min)/(l_max-l_min))+1
            strand.color = list(cm.afmhot(colorindex)[0:-1])
       
    con = Converter()
    con.dna_structure = dna_structure_length_colored
    con.write_cadnano_file(output_path+'_colorcoded_StapleLength.json')
        

    #%% compute domain length statitics
    dna_structure.compute_aux_data() # compute domains
    domain_lengths = []
    
    for domain in dna_structure.domain_list:
        if len(domain.sequence.translate(None, 'N'))> 0:
            domain_lengths.append(len(domain.sequence.translate(None, 'N'))) #length of domain, removed N for skipped bases
        
    # plot histogram of domain length
    l_min = min(domain_lengths)
    l_max = max(domain_lengths)
    dl = 1.0 # domain length increment
    
    fig = plt.figure()
    plt.hist(domain_lengths, bins=numpy.linspace(l_min-dl/2, l_max+dl/2, l_max-l_min+2) )# , bins=numpy.linspace(min(domain_lengths)-0.5, max(domain_lengths)+0.5, max(domain_lengths)-min(domain_lengths)+2))  # plt.hist passes it's arguments to np.histogram
    plt.title("Average domain length = "+str(round(numpy.mean(domain_lengths),4)) + " bases")
    plt.xlabel('Length of domain [bp]')
    plt.ylabel('Frequency')
    #plt.show()
    fig.savefig(output_path+'_stats_DomainLength.pdf')
    plt.close()
    
    # average domain length 
    avg_domain_length = []
    domain_max_melt = []
    domain_max_length = []
    number_of_domains = []
    for strand in dna_structure.strands:
        if not strand.is_scaffold:
            cur_len = []
            cur_mt = []
            for domain in strand.domain_list:
                cur_len.append(len(domain.sequence.translate(None, 'N')))
                if len(domain.sequence.translate(None, 'N'))>0:
                    cur_mt.append(MeltingTemp.Tm_NN(Seq(domain.sequence.translate(None, 'N'))))
                    
            domain_max_melt.append( max(cur_mt))
            avg_domain_length.append( numpy.mean(cur_len))
            domain_max_length.append( max(cur_len))
            number_of_domains.append(len(cur_mt)) # this will not incude polyT overhangs as separate domains
    
    fig = plt.figure()
    plt.hist(avg_domain_length, bins=numpy.linspace(l_min-dl/2, l_max+dl/2, l_max-l_min+2) )# , bins=numpy.linspace(min(domain_lengths)-0.5, max(domain_lengths)+0.5, max(domain_lengths)-min(domain_lengths)+2))  # plt.hist passes it's arguments to np.histogram
    plt.title("Average of average domain length = "+str(round(numpy.mean(avg_domain_length),4)) + " bases")
    plt.xlabel('Average staple-domain length [bp]')
    plt.ylabel('Frequency')
    #plt.show()
    fig.savefig(output_path+'_stats_AverageDomainLength.pdf')
    plt.close()
    
    fig = plt.figure()
    plt.hist(number_of_domains, 
             bins=numpy.linspace(min(number_of_domains)-0.5, max(number_of_domains)+0.5, max(number_of_domains)-min(number_of_domains)+2) )# , bins=numpy.linspace(min(domain_lengths)-0.5, max(domain_lengths)+0.5, max(domain_lengths)-min(domain_lengths)+2))  # plt.hist passes it's arguments to np.histogram
    plt.title("Average number of domains = "+str(round(numpy.mean(number_of_domains),4)) + " bases")
    plt.xlabel('Number of domains per staple')
    plt.ylabel('Frequency')
    #plt.show()
    fig.savefig(output_path+'_stats_NumberOfDomains.pdf')
    plt.close()
    #%%
    
    fig = plt.figure()
    plt.hist(domain_max_melt, bins=numpy.linspace(min(domain_max_melt)-0.5, max(domain_max_melt)+0.5, max(domain_max_melt)-min(domain_max_melt)+2))
    #plt.title("Domain melting-temperature distribution")
    plt.xlabel('Melting temp. of domain with highest melting temperature [deg]')
    plt.ylabel('Frequency')
    #plt.show()
    fig.savefig(output_path+'_stats_MaxDomainMeltingTemperature.pdf')
    plt.close()
    dmt_max = max(domain_max_melt)
    dmt_min = min(domain_max_melt)
    
    fig = plt.figure()
    plt.hist(domain_max_length, bins=numpy.linspace(min(domain_max_length)-0.5, max(domain_max_length)+0.5, max(domain_max_length)-min(domain_max_length)+2))
    #plt.hist(domain_max_length)
    #plt.title("Domain melting-temperature distribution")
    plt.xlabel('Length of longest domain')
    plt.ylabel('Frequency')
    #plt.show()
    fig.savefig(output_path+'_stats_MaxDomainLength.pdf')
    plt.close()
    
    #%%
    # write json colorcoded for melting temperature of domain with largest melting temperature  
    dmt_max = max(domain_max_melt)
    dmt_min = min(domain_max_melt)
    dna_structure_domainmelt_colored = dna_structure
    for strand in dna_structure_domainmelt_colored.strands:
        if not strand.is_scaffold:
            cur_mt = []
            for domain in strand.domain_list:
                if len(domain.sequence.translate(None, 'N'))>0:
                    cur_mt.append(MeltingTemp.Tm_NN(Seq(domain.sequence.translate(None, 'N'))))
            dmt = max(cur_mt)
            colorindex = int(255*(dmt-dmt_min)/(dmt_max-dmt_min))
            strand.color = list(cm.afmhot(colorindex)[0:-1])
        #print strand.color     
       
    con = Converter()
    con.dna_structure = dna_structure_domainmelt_colored
    con.write_cadnano_file(output_path+'_colorcoded_LargestDomainMelt.json')
  
    
    
    # write json colorcoded for length of longest domain
    dlen_min = min(domain_max_length)    
    dlen_max = max(domain_max_length)
    dna_structure_domainlength_colored = dna_structure
    for strand in dna_structure_domainlength_colored.strands:
        if not strand.is_scaffold:
            cur_len = []
            for domain in strand.domain_list:
                cur_len.append(len(domain.sequence.translate(None, 'N')))
            dlen = max(cur_len)
            colorindex = int(255*(dlen-dlen_min)/(dlen_max-dlen_min))
            strand.color = list(cm.afmhot(colorindex)[0:-1])
        #print strand.color     
       
    con = Converter()
    con.dna_structure = dna_structure_domainlength_colored
    con.write_cadnano_file(output_path+'_colorcoded_LongestDomainLength.json')


    # write json colorcoded for number of domains
    N_min = min(number_of_domains)    
    N_max = max(number_of_domains)
    dna_structure_out = dna_structure
    for strand in dna_structure_out.strands:
        if not strand.is_scaffold:
            cur_N = 0
            for domain in strand.domain_list:
                if len(domain.sequence.translate(None, 'N'))>0:
                    cur_N = cur_N +1
            colorindex = int(255*(cur_N-N_min)/(N_max-N_min))
            strand.color = list(cm.bwr(colorindex)[0:-1])
        #print strand.color     
       
    con = Converter()
    con.dna_structure = dna_structure_domainlength_colored
    con.write_cadnano_file(output_path+'_colorcoded_NumberOfDomains.json')
#%%

if __name__ == '__main__':
    main()