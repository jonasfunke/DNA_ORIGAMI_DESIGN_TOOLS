# -*- coding: utf-8 -*-
"""
@author: jean-philippe sobczak

Loads cadnano .json design, performs multiple quality checks and creates new
.json designs with rule-breaking staples coloured

first argument: cadnano json file path, relative to current working directory
second argument: scaffold type: M13mp18, p7308, p7560, p7704, p8064, p8100, p8634
"""

# compatibility to python 3.5
from __future__ import (absolute_import, division, print_function, unicode_literals)

import os
import sys

# copied from nanodesign/strand-statistics.py
# The following imports are designed to try and get the nanodesign package
# imported regardless of whether or not you have it installed in the
# site-packages. If you have it in site packages, the first block should just work.
# Otherwise, it will assume you can import it based on a relative path from this source
# file's directory, and try to do so by adjusting the system paths temporarily.
try:
    import nanodesign
except ImportError:
    base_path = os.path.abspath( os.path.join( os.path.dirname(os.path.abspath( __file__)), '../nanodesign/'))
    sys.path.append(base_path)
    import nanodesign
    # If the import fails now, we let the exception go all the way up to halt execution.
    sys.path = sys.path[:-1]
    
from nanodesign.converters import Converter

def read_file(file_name, seq_name):
    """ Read in a cadnano file. """
    converter = Converter()
    seq_file = None
    converter.read_cadnano_file(file_name, seq_file, seq_name)
    return converter


def main():

    # create path to caDNAno file to load
    current_working_directory = os.getcwd()
    file_name = sys.argv[1]
    file_name = os.path.normpath(current_working_directory + "/" + file_name)
    print('file_name to be opened' , file_name)
    
    # Set sequence to assign to scaffold.
    scaffold_seq_name = sys.argv[2]
    
    # Read cadnano file and create dna structure.
    converter = read_file(file_name, scaffold_seq_name)
    dna_structure = converter.dna_structure
   
    # determine domain information
    dna_structure.get_domains()
    
    # determine scaffold strand id, alert if multiple scaffold strands
    # prints starting location if not circular strand
    
    # strand id of last found scaffold strand
    scaffold_id = -1
    
    for strand in dna_structure.strands:
        if strand.is_scaffold:
            if scaffold_id != -1:
                print('multiple scaffold strands!')
            scaffold_id = strand.id
            if not strand.is_circular:
                print('scaffold strand starts in helix ', strand.tour[0].h, ' base number ', strand.tour[0].p)
            else:
                print('scaffold strand is circular, passes through helix ', strand.tour[0].h, ' base number ', strand.tour[0].p)
            
    # print double-stranded domain lengths of all staples and starting locations of staples
    for strand in dna_structure.strands:
        if not strand.is_scaffold:
            print('staple\t', strand.id, '\thelix\t', strand.tour[0].h, '\tbase\t', strand.tour[0].p, '\tdomains:\t', end = '')
            for domain in strand.domain_list:
                # check if domain is double-stranded
                if domain.base_list[0].across != None:
                    print( len( domain.base_list ), end = '\t')
            print('\n', end = '')


if __name__ == '__main__':
    main()