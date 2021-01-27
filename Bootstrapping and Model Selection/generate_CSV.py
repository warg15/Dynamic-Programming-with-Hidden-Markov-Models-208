#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 19:40:06 2020

@author: femw90
"""
import argparse
import numpy as np
import treeswift 
from treeswift import read_tree_newick
import scipy
from scipy.linalg import expm, sinm, cosm
import pickle
import random
import csv
import os
from os import listdir
from os.path import isfile, join

def likelihood(tree, seqs, gtr_probs, gtr_rates):
    
    def popLeafs(node):
        label = node.get_label()
        sequence = seqs[label]
        
        matrix = np.full((4, l), 0)
        i = 0
        for instance in sequence:
                if instance == 'A': matrix[0,i] = 1
                if instance == 'G': matrix[1,i] = 1
                if instance == 'C': matrix[2,i] = 1
                if instance == 'T': matrix[3,i] = 1
                i = i+1
        matrices[node] = matrix
    
    def updateMatrix(node):
        #matrices[node] = np.zeros((4,l))
        children = node.child_nodes()
        if len(children) == 2:
            child1 = children[0]
            child2 = children[1]
            branch1 = child1.get_edge_length()
            branch2 = child2.get_edge_length()
            Q1 = scipy.linalg.expm(R*branch1/d)
            Q2 = scipy.linalg.expm(R*branch2/d)
            
            mat1 = matrices[child1]
            mat2 = matrices[child2]
            matrix = np.zeros((4,l))
            
            for i in range(0,l): #for each item in the sequence
                for s in range(0,4): #for each row, aka for each A, G, C, T
                    #print(s)
                    L1 = np.transpose(Q1[s,:])@mat1[:,i]
                    L2 = np.transpose(Q2[s,:])@mat2[:,i]
                    matrix[s,i] = L1 * L2
            matrices[node] = matrix
        
        elif len(children) == 3:
            children = node.child_nodes()
            child1 = children[0]
            child2 = children[1]
            child3 = children[2]
            branch1 = child1.get_edge_length()
            branch2 = child2.get_edge_length()
            branch3 = child3.get_edge_length()
            Q1 = scipy.linalg.expm(R*branch1/d)
            Q2 = scipy.linalg.expm(R*branch2/d)
            Q3 = scipy.linalg.expm(R*branch3/d)
            
            mat1 = matrices[child1]
            mat2 = matrices[child2]
            mat3 = matrices[child3]
            matrix = np.zeros((4,l))
            
            for i in range(0,l): #for each item in the sequence
                for s in range(0,4): #for each row, aka for each A, G, C, T
                    #print(s)
                    L1 = np.transpose(Q1[s,:])@mat1[:,i]
                    L2 = np.transpose(Q2[s,:])@mat2[:,i]
                    L3 = np.transpose(Q3[s,:])@mat3[:,i]
                    matrix[s,i] = L1 * L2 * L3
            matrices[node] = matrix
        
        else:
            print('SOMETHING WRONG')
            print('SOMETHING WRONG')
            print('SOMETHING WRONG')
            print('SOMETHING WRONG')
            print('SOMETHING WRONG')
        
    prob_A, prob_C, prob_G, prob_T = gtr_probs # You can use these if it's more convenient
    rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG = gtr_rates # You can use these if it's more convenient
    log_likelihood_score = 0. # You can modify this variable, or you can use your own; either way is fine
    # TODO Your code here
    R = np.array([[-(rate_AG*prob_G + rate_AC*prob_C + rate_AT*prob_T), rate_AG*prob_G, rate_AC*prob_C, rate_AT*prob_T],
                  [rate_AG*prob_A, -(rate_AG*prob_A + rate_CG*prob_C + rate_GT*prob_T), rate_CG*prob_C, rate_GT*prob_T],
                  [rate_AC*prob_A, rate_CG*prob_G, -(rate_AC*prob_A + rate_CG*prob_G + rate_CT*prob_T), rate_CT*prob_T],
                  [rate_AT*prob_A, rate_GT*prob_G, rate_CT*prob_C, -(rate_AT*prob_A + rate_GT*prob_G + rate_CT*prob_C)]])
    
    d = -prob_A*R[0,0] -prob_G*R[1,1] -prob_C*R[2,2] -prob_T*R[3,3] #for normalizing the rates matris, R
    
    #list(test_dict.keys())[0] 
    l = len(seqs[list(seqs.keys())[0] ])   #length of sequences
    
    matrices = {}
    
    for node in tree.traverse_postorder():
        if node.is_leaf():popLeafs(node)
        else:
            updateMatrix(node)
            if node.is_root(): 
                root = node
                #print('FOUND THE ROOT')
    
    like = 0
    pis = np.array([prob_A, prob_G, prob_C, prob_T])
    rootMatrix = matrices[root]
    for i in range(0,l):
        temp = np.log(pis@rootMatrix[:,i])
        like = like + temp
        #print(like)
    #print(like)
    #like = np.log(like)
    #print(like)
    return like
    #return log_likelihood_score
def cut(num, eln):
    amt = eln
    value = math.trunc(amt * num) / amt
    #if value == 0.0: value = format(value, '.6f')
    #if value == 0.0: value = '0.000000'
    return value
    
def read_FASTA(filename):
    stream = open(filename); seqs = {}; name = None; seq = ''
    for line in stream:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                assert len(seq) != 0, "Malformed FASTA"
                seqs[name] = seq
            name = l[1:]
            assert name not in seqs, "Duplicate sequence ID: %s" % name
            seq = ''
        else:
            seq += l
    assert name is not None and len(seq) != 0, "Malformed FASTA"
    seqs[name] = seq; stream.close()
    return seqs

def calcP(i, liklihood):
    '''
    pickle = 'distribution_0' + str(i) + '.pickle'
    if i == 10: pickle = 'distribution_10.pickle'
    pickle_in = open(pickle,"rb")
    '''
    if i == 1: pickle_in = open('dist_01.pickle',"rb")
    if i == 2: pickle_in = open('dist_02.pickle',"rb")
    if i == 3: pickle_in = open('dist_03.pickle',"rb")
    if i == 4: pickle_in = open('dist_04.pickle',"rb")
    if i == 5: pickle_in = open('dist_05.pickle',"rb")
    if i == 6: pickle_in = open('dist_06.pickle',"rb")
    if i == 7: pickle_in = open('dist_07.pickle',"rb")
    if i == 8: pickle_in = open('dist_08.pickle',"rb")
    if i == 9: pickle_in = open('dist_09.pickle',"rb")
    if i == 10: pickle_in = open('dist_10.pickle',"rb")
    example_dict = pickle.load(pickle_in)
    #keys = list(example_dict.keys())
    #a = example_dict[keys[0]]
    dist_list = example_dict
        #dist_list = a.tolist()
    nullSize = len(dist_list)
    pVal = stats.t.sf(np.abs((liklihood - sum(dist_list)/nullSize)/stats.tstd(dist_list)), nullSize-1)
    return pVal

if __name__ == "__main__":
    
    #mypath = /Users/femw90/"Google Drive"/"SECKSY (1)"/"Grad School"/"Spring Quarter 2020"/"ECE 208"/"Homework 8"/hw8-sp20-warg15-master/example
    mypath = dir_path = os.path.dirname(os.path.realpath(__file__))
    faspath = join(mypath, 'example')
    fasfiles = [f for f in listdir(faspath) if isfile(join(faspath, f))]
    
    paramspath = join(mypath, 'params')
    paramsfiles = [f for f in listdir(paramspath) if isfile(join(paramspath, f))]
    
    treespath = join(mypath, 'trees')
    treesfiles = [f for f in listdir(treespath) if isfile(join(treespath, f))]
    
    listFas = []
    listParams = []
    listTrees = []
    
    for i in fasfiles:
        if i.endswith(".fas"):
            listFas.append(i)
    listFas.sort()
    
    for i in paramsfiles:
        if i.endswith(".txt"):
            listParams.append(i)
            
    for i in treesfiles:
        if i.endswith(".tre"):
            listTrees.append(i)
    
    myfile = join(treespath, treesfiles[0])
    tree = read_tree_newick(myfile)
    
    '''
    myfile = join(treespath, 'Exon.c3-unpartitioned.tre')
    treeDef = read_tree_newick(myfile)
    myfile = join(faspath, '04.fas')
    seqs = read_FASTA(myfile)
    myfile = join(faspath, 'gtr_params.txt')
    gtr_probs, gtr_rates = [[float(e.strip()) for e in l.strip().split()] for l in open(myfile)]
    value = likelihood(treeDef, seqs, gtr_probs, gtr_rates)
    print(value)
    '''
    storeTrees = dict()
    storeParams = dict()
    storeNormal = dict()
    
    for tree in listTrees:
        storeTrees[tree] = list()
    for param in listParams:
        storeParams[param] = list()
    for fas in listFas:
        storeNormal[fas] = list()
    
    i=0
    for fas in listFas: #go 1-10 in the fas files
        i=i+1
        print(i)
        myfile = join(faspath, 'tree.nwk')
        treeDef = read_tree_newick(myfile)
        myfile = join(faspath, fas)
        seqs = read_FASTA(myfile)
        myfile = join(faspath, 'gtr_params.txt')
        gtr_probs, gtr_rates = [[float(e.strip()) for e in l.strip().split()] for l in open(myfile)]
        value = likelihood(treeDef, seqs, gtr_probs, gtr_rates)
        #value = 69
        pVal = calcP(i, value)
        pVal = cut(pVal,10.0**6)
        storeNormal[fas].append(pVal)
        
        with open('testV6.csv', 'a', newline='') as file:
            writer = csv.writer(file, delimiter = ' ')
            writer.writerow([i, 'gtr_params.txt', 'tree.nwk', pVal ])
            file.close()
        
        for param in listParams:
            myfile = join(paramspath, param)
            gtr_probsTemp, gtr_ratesTemp = [[float(e.strip()) for e in l.strip().split()] for l in open(myfile)]
            value = likelihood(treeDef, seqs, gtr_probsTemp, gtr_ratesTemp)
            #value = 69
            pVal = calcP(i, value)
            pVal = cut(pVal,10.0**6)
            storeParams[param].append(pVal)
            
            with open('testV6.csv', 'a', newline='') as file:
                writer = csv.writer(file, delimiter = ' ')
                writer.writerow([i, param, 'tree.nwk', pVal ])
                file.close()
        
        for tree in listTrees:
            myfile = join(treespath, tree)
            treeTemp = read_tree_newick(myfile)
            value = likelihood(treeTemp, seqs, gtr_probs, gtr_rates)
            #value = 69
            pVal = calcP(i, value)
            pVal = cut(pVal,10.0**6)
            storeTrees[tree].append(pVal)
            
            with open('testV6.csv', 'a', newline='') as file:
                writer = csv.writer(file, delimiter = ' ')
                writer.writerow([i, 'gtr_params.txt', tree, pVal ])
                file.close()
                
                
    
    
    pickle_out = open("storeTrees.pickle","wb")
    pickle.dump(storeTrees, pickle_out)
    pickle_out.close()
    
    pickle_out = open("storeParams.pickle","wb")
    pickle.dump(storeParams, pickle_out)
    pickle_out.close()
    
    pickle_out = open("storeNormal.pickle","wb")
    pickle.dump(storeNormal, pickle_out)
    pickle_out.close()
    
    
    
    
    
    
    
    