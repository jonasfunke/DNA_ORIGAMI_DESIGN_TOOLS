#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 17:01:56 2017

@author: jonasfunke
"""

#%%

#from __future__ import (absolute_import, division, print_function, unicode_literals)
import __future__
import numpy #for writing matrices
import matplotlib.pyplot as plt #plotting 


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
#compute contact matrix
cur_seq = 'ACTGCCTAGCATCAGT'#reporter_sequences[0];

contacts = numpy.zeros((len(cur_seq),len(cur_seq)))
for i in range(len(cur_seq)):
    for j in range(i+1, len(cur_seq)):
        if (cur_seq[i]==complement[cur_seq[j]]) and abs(i-j)>2:
            contacts[i][j] = 1
                    
for i in range(len(contacts)):
    for j in range(i+1, len(contacts)):
        if contacts[][]
        

plt.matshow(contacts)
plt.show()        


