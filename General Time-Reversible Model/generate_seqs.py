#!/usr/bin/env python3
ERROR_GTR_PARAMS_FILE = "Invalid GTR parameters file"
ERROR_ROOT_SEQ = "Invalid root sequence. Must be path to FASTA file or integer"

import numpy as np
import treeswift 
import scipy
from scipy.linalg import expm, sinm, cosm

def random_seq(k, gtr_probs):
    '''
    This function generates a random sequence of length k using GTR stationary probabilities
    :param k: The length of the output sequence
    :param gtr_probs: The GTR stationary probabilities as a list [prob_A, prob_C, prob_G, prob_T]
    :return: A random string of length k generated using gtr_probs
    '''
    
    prob_A, prob_C, prob_G, prob_T = gtr_probs # You can use these if it's more convenient
    # TODO Your code here
    lst = ['A', 'C', 'G', 'T']; answer = np.random.choice(lst, size = k, p = gtr_probs); strg = ""
    for i in answer:  
        strg += i
    return strg#""

def evolve(tree, root_seq, gtr_probs, gtr_rates):
    def normalize(matrix):
        matrix[0,:] = matrix[0,:] / sum(matrix[0,:])
        matrix[1,:] = matrix[1,:] / sum(matrix[1,:])
        matrix[2,:] = matrix[2,:] / sum(matrix[2,:])
        matrix[3,:] = matrix[3,:] / sum(matrix[3,:])
        return matrix
        
    def new_Seq(node):
        node_label = node.get_label()
        seq = root_seq
        if isinstance(node_label, str):
            parent = node.get_parent()
            parent_label = parent.get_label()
            seq = ''
            branch = node.get_edge_length()
            newR = normalize(scipy.linalg.expm(R*branch/d))
            for instance in parent_label:
                if instance == 'A': char = np.random.choice(['A', 'G', 'C', 'T'], p = newR[0,:])
                if instance == 'G': char = np.random.choice(['A', 'G', 'C', 'T'], p = newR[1,:])
                if instance == 'C': char = np.random.choice(['A', 'G', 'C', 'T'], p = newR[2,:])
                if instance == 'T': char = np.random.choice(['A', 'G', 'C', 'T'], p = newR[3,:])
                seq = seq + char
        return seq
    '''
    This function simulates the evolution of a root sequence down a given tree under the GTR model
    :param tree: A TreeSwift Tree object representing the phylogenetic tree
    :param root_seq: The root sequence to evolve down the tree
    :param gtr_probs: The GTR stationary probabilities as a list [prob_A, prob_C, prob_G, prob_T]
    :param gtr_rates: The GTR transition rates as a list [rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG]
    :return: A dictionary where keys are the labels of the given tree and values are evolved sequences
    '''
    prob_A, prob_C, prob_G, prob_T = gtr_probs # You can use these if it's more convenient
    rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG = gtr_rates # You can use these if it's more convenient
    seqs = dict() # This will be your output (keys = leaf labels (str) and values = sequences (str))
    # TODO Your code here
    R = np.array([[-(rate_AG*prob_G + rate_AC*prob_C + rate_AT*prob_T), rate_AG*prob_G, rate_AC*prob_C, rate_AT*prob_T],
                  [rate_AG*prob_A, -(rate_AG*prob_A + rate_CG*prob_C + rate_GT*prob_T), rate_CG*prob_C, rate_GT*prob_T],
                  [rate_AC*prob_A, rate_CG*prob_G, -(rate_AC*prob_A + rate_CG*prob_G + rate_CT*prob_T), rate_CT*prob_T],
                  [rate_AT*prob_A, rate_GT*prob_G, rate_CT*prob_C, -(rate_AT*prob_A + rate_GT*prob_G + rate_CT*prob_C)]])
    
    d = -prob_A*R[0,0] -prob_G*R[1,1] -prob_C*R[2,2] -prob_T*R[3,3] #for normalizing the rates matris, R
    
    for node in tree.traverse_preorder():
        if node.is_root(): 
            node.set_label(root_seq)
        else:
            c = new_Seq(node)
            label = node.get_label()
            if node.is_leaf():
                seqs[label] = c
            else:
                node.set_label(c)
    return seqs

if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your functions, and outputting the results
    '''
    # parse user arguments
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=False, type=str, default='stdin', help="Input Tree File (Newick)")
    parser.add_argument('-p', '--gtr_params', required=True, type=str, help="GTR Parameters File")
    parser.add_argument('-r', '--root_seq', required=True, type=str, help="Root Sequence (FASTA or integer length)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Sequence File (FASTA)")
    args = parser.parse_args()

    # load input tree
    from treeswift import read_tree_newick
    if args.tree == 'stdin':
        from sys import stdin
        tree = read_tree_newick(stdin)
    else:
        tree = read_tree_newick(args.tree)

    # load GTR parameters
    try:
        gtr_probs, gtr_rates = [[float(e.strip()) for e in l.strip().split()] for l in open(args.gtr_params)]
    except:
        raise ValueError(ERROR_GTR_PARAMS_FILE)

    # load (or generate) root sequence
    from os.path import isfile
    try:
        if isfile(args.root_seq):
            root_seq = ''.join([l.strip() for l in open(args.root_seq)][1:]).upper()
            assert set(root_seq) == {'A','C','G','T'}
        else:
            root_seq = random_seq(int(args.root_seq), gtr_probs).upper()
    except:
        raise ValueError(ERROR_ROOT_SEQ)

    # run student code and check output
    seqs = evolve(tree, root_seq, gtr_probs, gtr_rates)
    leaves = {leaf.label for leaf in tree.traverse_leaves()}
    for l in leaves:
        if l not in seqs:
            raise RuntimeError("No sequence for leaf: %s" % l)
    for l in seqs:
        if l not in leaves:
            raise RuntimeError("Extra sequence in output: %s" % l)

    # output resulting sequences
    if args.output == 'stdout':
        from sys import stdout as outfile
    else:
        outfile = open(args.output,'w')
    for pair in seqs.items():
        outfile.write('>%s\n%s\n' % pair)
    outfile.close()