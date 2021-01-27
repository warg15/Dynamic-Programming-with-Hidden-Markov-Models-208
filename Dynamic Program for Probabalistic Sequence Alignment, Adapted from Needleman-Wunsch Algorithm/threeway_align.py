#!/usr/bin/env python3

import numpy as np
import math

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

def threeway_align(s1, s2, s3, indel_rate):
    '''
    This function computes the threeway sequence alignment between strings s1, s2, and s3 using a relative indel rate of indel_rate
    :param s1: string s1
    :param s2: string s2
    :param s3: string s3
    :param indel_rate: the relative indel rate
    '''
    
    gap_Penalty = int(2*math.log2(indel_rate*19/20))
    #gap_Penalty = -13
    gap_Penalty_two = 2*gap_Penalty
    sz_s1 = len(s1)+1
    sz_s2 = len(s2)+1
    sz_s3 = len(s3)+1
    
    empty = -10101010101010
    cube = np.full(( sz_s1, sz_s2, sz_s3), empty, dtype=int)
    cube[:,0,0] = np.arange(0,gap_Penalty_two*sz_s1, gap_Penalty_two )
    cube[0,:,0] = np.arange(0,gap_Penalty_two*sz_s2, gap_Penalty_two )
    cube[0,0,:] = np.arange(0,gap_Penalty_two*sz_s3, gap_Penalty_two )
    
    #fill in the sides first
    for j in range (1,sz_s1): #s1
        for k in range (1,sz_s2): #s2
            cube[j,k,0] = max(
                    cube[j-1,k-1,0] + BLOSUM62[AMINOS.index(s1[j-1])][AMINOS.index(s2[k-1])]+gap_Penalty_two,
                    cube[j-1,k,0] + gap_Penalty_two,
                    cube[j,k-1,0] + gap_Penalty_two
                    )
    for j in range (1,sz_s1): #s1
        for k in range (1,sz_s3): #s3
            cube[j,0,k] = max(
                    cube[j-1,0, k-1] + BLOSUM62[AMINOS.index(s1[j-1])][AMINOS.index(s3[k-1])]+gap_Penalty_two,
                    cube[j-1,0,k] + gap_Penalty_two,
                    cube[j,0,k-1] + gap_Penalty_two
                    )
    for j in range (1,sz_s2): #s2
        for k in range (1,sz_s3): #s3
            cube[0,j,k] = max(
                    cube[0,j-1,k-1] + BLOSUM62[AMINOS.index(s2[j-1])][AMINOS.index(s3[k-1])]+gap_Penalty_two,
                    cube[0,j-1,k] + gap_Penalty_two,
                    cube[0,j,k-1] + gap_Penalty_two
                    )
    
    diff = 2*max( #set a maximum amount of deviation from the diagonal in any direction
            abs(sz_s1-sz_s2),
            abs(sz_s2-sz_s1),
            abs(sz_s3-sz_s2),
            5
            )
    for i in range (1,sz_s1):
        #print(i)
        for j in range (1,sz_s2):
            if (abs(i-j)< diff):
                for k in range (1,sz_s3):
                    if (abs(i-k)< diff) and (abs(j-k)< diff):
                        cube[i,j,k] = max(
                            cube[i,j-1,k] + gap_Penalty_two, #1
                            cube[i,j-1, k-1] + BLOSUM62[AMINOS.index(s2[j-1])][AMINOS.index(s3[k-1])]+ gap_Penalty_two, #2
                            cube[i, j, k-1] + gap_Penalty_two, #3
                            cube[i-1, j-1, k] + BLOSUM62[AMINOS.index(s1[i-1])][AMINOS.index(s2[j-1])] + gap_Penalty_two, #4
                            cube[i-1,j-1,k-1] + BLOSUM62[AMINOS.index(s1[i-1])][AMINOS.index(s2[j-1])] + BLOSUM62[AMINOS.index(s2[j-1])][AMINOS.index(s3[k-1])] + BLOSUM62[AMINOS.index(s1[i-1])][AMINOS.index(s3[k-1])], #5
                            cube[i-1,j, k-1] + BLOSUM62[AMINOS.index(s1[i-1])][AMINOS.index(s3[k-1])] + gap_Penalty_two,#6
                            cube[i-1, j, k] + gap_Penalty_two, #7
                            )
    x = sz_s1-1
    y = sz_s2-1
    z = sz_s3-1
    #print('fat')
    ##########################################################################################################################################
    #now do backtracking to find alignment
    s1new=""
    s2new=""
    s3new=""
    #check all 7 movement posibilities and figure out which one happened
    while x > 0 and y > 0 and z > 0:
        xyzCell = cube[x,y,z]
        #print(x,y,z)
        if xyzCell == cube[x,y-1,z] + gap_Penalty_two: #1
            s1new = '-'+s1new
            s2new = s2[y-1]+s2new
            s3new = '-'+s3new
            y=y-1
        elif xyzCell == cube[x,y-1, z-1] + BLOSUM62[AMINOS.index(s2[y-1])][AMINOS.index(s3[z-1])]+ gap_Penalty_two: #2
            s1new = '-'+s1new
            s2new = s2[y-1]+s2new
            s3new = s3[z-1]+s3new
            y=y-1
            z=z-1
        elif xyzCell == cube[x, y, z-1] + gap_Penalty_two: #3
            s1new = '-'+s1new
            s2new = '-'+s2new
            s3new = s3[z-1]+s3new
            z=z-1
        elif xyzCell == cube[x-1, y-1, z] + BLOSUM62[AMINOS.index(s1[x-1])][AMINOS.index(s2[y-1])] + gap_Penalty_two: #4
            s1new = s1[x-1]+s1new
            s2new = s2[y-1]+s2new
            s3new = '-'+s3new
            y=y-1
            x=x-1
        elif xyzCell == cube[x-1,y-1,z-1] + BLOSUM62[AMINOS.index(s1[x-1])][AMINOS.index(s2[y-1])] + BLOSUM62[AMINOS.index(s2[y-1])][AMINOS.index(s3[z-1])] + BLOSUM62[AMINOS.index(s1[x-1])][AMINOS.index(s3[z-1])]: #5
            s1new = s1[x-1]+s1new
            s2new = s2[y-1]+s2new
            s3new = s3[z-1]+s3new
            y=y-1
            x=x-1
            z=z-1
        elif xyzCell == cube[x-1,y, z-1] + BLOSUM62[AMINOS.index(s1[x-1])][AMINOS.index(s3[z-1])] + gap_Penalty_two:#6
            s1new = s1[x-1]+s1new
            s2new = '-'+s2new
            s3new = s3[z-1]+s3new
            x=x-1
            z=z-1
        elif xyzCell == cube[x-1, y, z] + gap_Penalty_two: #7
            s1new = s1[x-1]+s1new
            s2new = '-'+s2new
            s3new = '-'+s3new
            x=x-1
    #print('fat')
    while x > 0 and y > 0:
        xyzCell = cube[x,y,z]
        if xyzCell == cube[x,y-1,z] + gap_Penalty_two: #1
            s1new = '-'+s1new
            s2new = s2[y-1]+s2new
            s3new = '-'+s3new
            y=y-1
        elif xyzCell == cube[x-1, y-1, z] + BLOSUM62[AMINOS.index(s1[x-1])][AMINOS.index(s2[y-1])] + gap_Penalty_two: #4
            s1new = s1[x-1]+s1new
            s2new = s2[y-1]+s2new
            s3new = '-'+s3new
            y=y-1
            x=x-1
        elif xyzCell == cube[x-1, y, z] + gap_Penalty_two: #7
            s1new = s1[x-1]+s1new
            s2new = '-'+s2new
            s3new = '-'+s3new
            x=x-1
    while x > 0 and z > 0:
        xyzCell = cube[x,y,z]
        if xyzCell == cube[x, y, z-1] + gap_Penalty_two: #3
            s1new = '-'+s1new
            s2new = '-'+s2new
            s3new = s3[z-1]+s3new
            z=z-1
        elif xyzCell == cube[x-1,y, z-1] + BLOSUM62[AMINOS.index(s1[x-1])][AMINOS.index(s3[z-1])] + gap_Penalty_two:#6
            s1new = s1[x-1]+s1new
            s2new = '-'+s2new
            s3new = s3[z-1]+s3new
            x=x-1
            z=z-1
        elif xyzCell == cube[x-1, y, z] + gap_Penalty_two: #7
            s1new = s1[x-1]+s1new
            s2new = '-'+s2new
            s3new = '-'+s3new
            x=x-1
    while y > 0 and y > 0:
        xyzCell = cube[x,y,z]
        if xyzCell == cube[x,y-1,z] + gap_Penalty_two: #1
            s1new = '-'+s1new
            s2new = s2[y-1]+s2new
            s3new = '-'+s3new
            y=y-1
        elif xyzCell == cube[x,y-1, z-1] + BLOSUM62[AMINOS.index(s2[y-1])][AMINOS.index(s3[z-1])]+ gap_Penalty_two: #2
            s1new = '-'+s1new
            s2new = s2[y-1]+s2new
            s3new = s3[z-1]+s3new
            y=y-1
            z=z-1
        elif xyzCell == cube[x, y, z-1] + gap_Penalty_two: #3
            s1new = '-'+s1new
            s2new = '-'+s2new
            s3new = s3[z-1]+s3new
            z=z-1
    while z > 0:
        s1new = '-'+s1new
        s2new = '-'+s2new
        s3new = s3[z-1]+s3new
        z=z-1
    while y > 0:
        s1new = '-'+s1new
        s2new = s2[y-1]+s2new
        s3new = '-'+s3new
        y=y-1
    while x > 0:
        s1new = s1[x-1]+s1new
        s2new = '-'+s2new
        s3new = '-'+s3new
        x=x-1
    return s1new,s2new,s3new # TODO replace this with your code

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
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_sequences', required=False, type=str, default='stdin', help="Input Sequences (FASTA format)")
    parser.add_argument('-r', '--indel_rate', required=True, type=float, help="Indel Rate")
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
    assert len(seqs) == 3, "Input FASTA file must have exactly 3 sequences"
    k1,k2,k3 = sorted(seqs.keys())
    a1,a2,a3 = threeway_align(seqs[k1], seqs[k2], seqs[k3], args.indel_rate)
    assert len(a1) == len(a2) == len(a3), "Aligned sequences must have equal length"
    assert len(a1) != 0, "Returned empty strings"
    for seq,orig in [(a1,seqs[k1]), (a2,seqs[k2]), (a3,seqs[k3])]:
        assert seq.replace('-','') == orig, "Removing gaps from the aligned sequences must yield the original sequences"
    for ID,seq in [(k1,a1), (k2,a2), (k3,a3)]:
        outfile.write(">%s\n%s\n" % (ID,seq))
    outfile.close()
