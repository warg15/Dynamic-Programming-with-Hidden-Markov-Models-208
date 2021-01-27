#!/usr/bin/env python3
AMINOS   = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
BLOSUM62 = [[ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2],
            [ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
            [-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3],
            [-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2],
            [-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3],
            [ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3],
            [-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2],
            [-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1],
            [-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2],
            [-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1],
            [-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1],
            [-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2],
            [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3],
            [-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1],
            [-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2],
            [ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2],
            [ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2],
            [ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1],
            [-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2],
            [-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2, 7]]


def twoway_align(s1, s2, indel_rate):
    import math
    import numpy as np
    #gap_penalty = 2 * np.log2(indel_rate)
    gap_penalty = math.log2(0.92*indel_rate)
    m = len(s1)
    n = len(s2)

    E1 = np.zeros((m+1,n+1))  
    for i in range(0, m + 1):
        E1[i][0] = gap_penalty * i
    do = 0
    for j in range(0, n + 1):
        E1[0][j] = gap_penalty * j
    for i in range(1, m+1):
        for j in range(1, n+1):
            E1[i][j] = max(E1[i][j-1]+gap_penalty, E1[i-1][j]+gap_penalty, E1[i-1][j-1] + BLOSUM62[AMINOS.index(s1[i-1])][AMINOS.index(s2[j-1])])
            do += E1[i][j]
            
    a1 = '' #s1, the start sequence
    a2 = ''
    i=m
    j=n
    
    gapS = []
    gapE = []
    
    while i>0 and j>0:
        if E1[i][j] == E1[i-1][j-1] + BLOSUM62[AMINOS.index(s1[i-1])][AMINOS.index(s2[j-1])]:
            a1 =a1+ s1[i-1]
            a2 =a2+ s2[j-1]
            i =i- 1
            j =j- 1
        elif E1[i][j] == E1[i-1][j] + gap_penalty:
            a1 =a2+ s1[i-1]
            a2 =a2+ '-'
            gapE.append(len(a2)-1)
            i =i- 1
        elif E1[i][j] == E1[i][j-1] + gap_penalty:
            a1 =a1+ '-'
            gapS.append(len(a1)-1)
            a2 =a2+ s2[j-1]
            j =j- 1
    while i > 0:
        a1 =a1+ s1[i-1]
        a2 =a2+ '-'
        gapE.append(len(a2)-1)
        i =i- 1
    while j > 0:
        a1 =a1+ '-'
        gapS.append(len(a1)-1)
        a2 =a2+ s2[j-1]
        j =j- 1

    #a1 = a1[::-1]
    #a2 = a2[::-1]
    return a1, a2, gapS, gapE

def MSA(seqs,tree,indel_rate):
    '''
    This function computes the MSA for the input sequences
    :param s1: seqs is a dictionary where keys are species names and values are amino acid sequences
    :param s2: tree is a TreeSwift object
    :param s3: indel_rate is a floating point number giving the relative rate of indel to substitutions
    '''
    import math
    
    ##########################################
    ##########################################
    ##########################################
    ##########################################
    
    import copy
    
    gap_penalty = math.log2(0.92*indel_rate)
    AMINOS.append('-')
    for a in BLOSUM62:
        a.append(gap_penalty)
    BLOSUM62.append([gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,gap_penalty,0])
    
    
    unaligned = []
    aligned = []
    
    for node in tree.traverse_leaves():
        unaligned.append(node)
    
    closest_pair = [node, node, 9999999]
    
    numLeaves = 0
    
    
    for nodeS in tree.traverse_leaves():
        numLeaves=numLeaves+1
        for nodeE in tree.traverse_leaves():
            if nodeS != nodeE:
                dist = tree.distance_between(nodeS, nodeE)
                if dist < closest_pair[2]:
                    closest_pair[0] = nodeS
                    closest_pair[1] = nodeE
                    closest_pair[2] = dist
    
    unaligned.remove(closest_pair[0])
    aligned.append(closest_pair[0])
    
    used = []
    
    gaps_dict = {}
    for node in tree.traverse_leaves():
        gaps_dict[node.get_label()] = list()
    
    #call 2-way align on start and end
    
    distances = tree.distance_matrix()
    
    for i in range(numLeaves-1):
        print(i)
        closest_pair = [node, node, 9999999]
        for j in aligned:
            #if j not in used:
                 for k in unaligned:
                     dist = distances[j][k]
                     if dist < closest_pair[2]:
                         closest_pair[0] = j #start
                         closest_pair[1] = k #end
                         closest_pair[2] = dist
        start = seqs[closest_pair[0].get_label()]
        end = seqs[closest_pair[1].get_label()]
        #print(closest_pair[0].get_label(), closest_pair[1].get_label())
        s1, s2, gapS, gapE = twoway_align(start, end, indel_rate)

        if len(gaps_dict[closest_pair[0].get_label()]) > 0:
            gaps_dict[closest_pair[1].get_label()].extend(gaps_dict[closest_pair[0].get_label()][:])
            #gaps_dict[closest_pair[1].get_label()].append(',')
        if len(gapE) > 0:
            gaps_dict[closest_pair[1].get_label()].extend(gapE)
            gaps_dict[closest_pair[1].get_label()].append(',')
        if len(gapS) > 0:
            for node in aligned:
                label = node.get_label()
                gaps_dict[label].extend(gapS)
                gaps_dict[label].append(',')
        
        
        #seqs[closest_pair[0].get_label()] = s1
        #seqs[closest_pair[1].get_label()] = s2
        unaligned.remove(closest_pair[1])
        aligned.append(closest_pair[1])
        used.append(closest_pair[0])
        
    
    gaps_add = copy.deepcopy(gaps_dict)
    for label in gaps_add.keys():
        for j in range(len(gaps_add[label])):
            if gaps_add[label][j] != ',':
                gaps_add[label][j]  = 0
    
    
    
    for label in gaps_dict.keys():
        print (label)
        for i in range(len(gaps_dict[label])):
            if gaps_dict[label][i] != ',':
                insertion_index = gaps_dict[label][i] + gaps_add[label][i]
                seqs[label] = seqs[label][0:insertion_index] + '-' + seqs[label][insertion_index::]
                indicator = 0
                for j in range (i,len(gaps_dict[label])):
                    if gaps_dict[label][j] != ',':
                        if indicator == 1:
                            if insertion_index < gaps_dict[label][j]:
                                gaps_add[label][j] = gaps_add[label][j] + 1
                    if gaps_dict[label][j] == ',':
                        indicator = 1
                
    
    
    
    aln = seqs
    return aln

def read_FASTA(stream):
    '''
    This function reads a FASTA file from a given stream and returns a dictionary mapping identifiers to sequences
    '''
    seqs = {}; name = None; seq = ''
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
    seqs[name] = seq
    return seqs

if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your function, and printing the output
    '''
    import argparse
    from treeswift import *
    
    
    
    
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_sequences', required=False, type=str, default='stdin', help="Input Sequences (FASTA format)")
    parser.add_argument('-r', '--indel_rate', required=True, type=float, help="Indel Rate")
    parser.add_argument('-t', '--guidance_tree', required=True, type=str, help="Guidance tree (Newick format)")
    parser.add_argument('-o', '--output_alignment', required=False, type=str, default='stdout', help="Output Alignment (FASTA format)")
    args = parser.parse_args()
    assert args.indel_rate >= 0, "Relative indel rate must be non-negative"
    
    if args.input_sequences == 'stdin':
        from sys import stdin as infile
    else:
        infile = open(args.input_sequences)
    if args.output_alignment == 'stdout':
        from sys import stdout as outfile
    else:
        outfile = open(args.output_alignment,'w')
    
    seqs = read_FASTA(infile); infile.close()
    tree = read_tree_newick(args.guidance_tree)
    indel_rate = args.indel_rate

    aln = MSA(seqs,tree,indel_rate)
    for ID in aln:
        outfile.write(">%s\n%s\n" % (ID,aln[ID]))
    outfile.close()
    
    
    
    
    
