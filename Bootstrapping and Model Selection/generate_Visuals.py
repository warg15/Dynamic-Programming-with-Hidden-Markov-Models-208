#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 21:37:24 2020

@author: femw90
"""

import pickle
import numpy as np
from matplotlib.pylab import plt
import statistics

if __name__ == "__main__":
    
    
    '''
    seqs = read_FASTA(infile); infile.close()
    pickle_out = open("dict.pickle","wb")
    pickle.dump(seqs, pickle_out)
    pickle_out.close()
    '''
    
    
    
    pickle_in = open("storeNormal.pickle","rb")
    storeNormal = pickle.load(pickle_in)
    #keys = list(example_dict.keys())
    #a = example_dict[keys[0]]
    
    pickle_in = open("storeParams.pickle","rb")
    storeParams = pickle.load(pickle_in)
    
    
    pickle_in = open("storeTrees.pickle","rb")
    storeTrees = pickle.load(pickle_in)
    
    
    #Filter through the dictionaries, recording any key 
    #that has a p-value anywhere in the list of less than 0.05
    allTrees = list(storeTrees.keys())
    goodTrees = list(storeTrees.keys())
    badtrees = set()
    for tree in storeTrees.keys():
        for i in range(10):
            if storeTrees[tree][i] <= 0.05:
                badtrees.add(tree)
                #print(storeTrees[tree][i])
    for tree in badtrees:
        goodTrees.remove(tree)
    badtreesShow = list(badtrees)
    
    
    allParams = list(storeParams.keys())
    goodParams = list(storeParams.keys())
    badparams = set()
    for param in storeParams.keys():
        for i in range(10):
            if storeParams[param][i] <= 0.05:
                badparams.add(param)
                #print(storeParams[param][i])
    for param in badparams:
        goodParams.remove(param)
    badparamShow = list(badparams)
    
    
    
    averageTrees = dict()
    averageTreesList = []
    for tree in storeTrees.keys():
        averageTrees[tree] = sum(storeTrees[tree]) / len(storeTrees[tree])
        averageTreesList.append(averageTrees[tree])
        
    averageParams = dict()
    averageParamsList = []
    for param in storeParams.keys():
        averageParams[param] = sum(storeParams[param]) / len(storeParams[param])
        averageParamsList.append(averageParams[param])
    
    paramRange = list(range(1,6))
    
    plt.scatter(paramRange, averageParamsList)
    plt.xticks(range(len(allParams)), allParams, size='small', rotation=70 )
    plt.title ('Average P-Value of parameters over all 10 FASTA Files')
    plt.ylabel ('P-Value')
    plt.grid()
    plt.yticks(np.arange(0, 0.5, .03))
    plt.tight_layout()
    plt.savefig('Parameters.png', dpi = 600)
    plt.show()
    
    treeRange = list(range(1,28))
    plt.scatter(treeRange, averageTreesList)
    plt.xticks(range(len(allTrees)), allTrees, size='small', rotation=90 )
    plt.title ('Average P-Value of trees over all 10 FASTA Files')
    plt.ylabel ('P-Value')
    plt.grid()
    plt.yticks(np.arange(0, 0.5, .05))
    plt.tight_layout()
    plt.savefig('Trees.png', dpi = 600)
    plt.show()


    