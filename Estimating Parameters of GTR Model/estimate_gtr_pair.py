#!/usr/bin/env python3
ERROR_PROB = "Stationary probabilities must be between 0 and 1"
ERROR_SEQS = "Invalid sequence file"
ERROR_TWO_SEQS = "Sequence file must have exactly 2 sequences"
PROB_KEYS = ['A', 'C', 'G', 'T']
RATE_KEYS = ['CT', 'AT', 'GT', 'AC', 'CG', 'AG']

import numpy as np
import treeswift 
import scipy
from scipy.linalg import expm, sinm, cosm, logm

def gtr_params_pair(r, s, d):
    '''
    def normalize(matrix):
        matrix[0,:] = matrix[0,:] / sum(matrix[0,:])
        matrix[1,:] = matrix[1,:] / sum(matrix[1,:])
        matrix[2,:] = matrix[2,:] / sum(matrix[2,:])
        matrix[3,:] = matrix[3,:] / sum(matrix[3,:])
        return matrix
    '''
    '''
    This function estimates GTR model parameters from a tree and sequences
    :param r: A sequence of length ``k``
    :param s: A sequence of length ``k``
    :param d: The pairwise distance between ``r`` and ``s``
    :return: A dictionary ``gtr_probs`` storing the GTR stationary probabilities, and a dictionary ``gtr_rates`` storing the GTR transition rates
    '''
    t=d
    gtr_probs = dict() # keys: {'A', 'C', 'G', 'T'}   values: the corresponding GTR stationary probabilities
    gtr_rates = dict() # keys: {'AC', 'AG', 'AT', 'CG', 'CT', 'GT'}   values: the corresponding GTR transition rates
    # TODO Your code here
    length = len(s)
    amountLetters = [0.0,0.0,0.0,0.0] #A C G T
    for i in r:
        if i == 'A': amountLetters[0] = amountLetters[0]+1
        if i == 'G': amountLetters[1] = amountLetters[1]+1
        if i == 'C': amountLetters[2] = amountLetters[2]+1
        if i == 'T': amountLetters[3] = amountLetters[3]+1
        
    for i in s:
        if i == 'A': amountLetters[0] = amountLetters[0]+1
        if i == 'G': amountLetters[1] = amountLetters[1]+1
        if i == 'C': amountLetters[2] = amountLetters[2]+1
        if i == 'T': amountLetters[3] = amountLetters[3]+1
    normLetters = [float(i)/sum(amountLetters) for i in amountLetters]
    gtr_probs['A'] = normLetters[0]
    gtr_probs['G'] = normLetters[1]
    gtr_probs['C'] = normLetters[2]
    gtr_probs['T'] = normLetters[3]
    
    transCount = dict() # keys: {'AA', 'CC', 'GG', 'TT', 'AC', 'AG', 'AT', 'CG', 'CT', 'GT'}
    transCount['AA'] = 0.0; transCount['AC'] = 0.0; transCount['AG'] = 0.0; transCount['AT'] = 0.0
    transCount['GA'] = 0.0; transCount['GC'] = 0.0; transCount['GG'] = 0.0; transCount['GT'] = 0.0
    transCount['CA'] = 0.0; transCount['CC'] = 0.0; transCount['CG'] = 0.0; transCount['CT'] = 0.0
    transCount['TA'] = 0.0; transCount['TC'] = 0.0; transCount['TG'] = 0.0; transCount['TT'] = 0.0
    for k in range(length): 
        rChar = r[k]
        sChar = s[k]
        keyString = rChar + sChar
        transCount[keyString] = transCount[keyString] + 1.0
        keyString = sChar + rChar
        transCount[keyString] = transCount[keyString] + 1.0
        
    '''
    newR = np.array([[transCount['AA'], transCount['AG'], transCount['AC'], transCount['AT'] ],
                  [transCount['GA'], transCount['GG'], transCount['GC'], transCount['GT'] ],
                  [transCount['CA'], transCount['CG'], transCount['CC'], transCount['CT'] ],
                  [transCount['TA'], transCount['TG'], transCount['TC'], transCount['TT'] ]])
    '''  
    newerR = np.array([[transCount['AA']/amountLetters[0], transCount['AG']/amountLetters[0], transCount['AC']/amountLetters[0], transCount['AT']/amountLetters[0] ],
                  [transCount['GA']/amountLetters[1], transCount['GG']/amountLetters[1], transCount['GC']/amountLetters[1], transCount['GT']/amountLetters[1] ],
                  [transCount['CA']/amountLetters[2], transCount['CG']/amountLetters[2], transCount['CC']/amountLetters[2], transCount['CT']/amountLetters[2] ],
                  [transCount['TA']/amountLetters[3], transCount['TG']/amountLetters[3], transCount['TC']/amountLetters[3], transCount['TT']/amountLetters[3] ]])
    
    #temp = np.copy(newR)
    #newRnorm = normalize(temp)
    R = logm(newerR)/t
    
    gtr_rates['AC'] = ((R[0,2]/gtr_probs['C']) + (R[2,0]/gtr_probs['A']))/2
    gtr_rates['AG'] = ((R[0,1]/gtr_probs['G']) + (R[1,0]/gtr_probs['A']))/2
    gtr_rates['AT'] = ((R[0,3]/gtr_probs['T']) + (R[3,0]/gtr_probs['A']))/2
    gtr_rates['CG'] = ((R[2,1]/gtr_probs['G']) + (R[1,2]/gtr_probs['C']))/2
    gtr_rates['CT'] = ((R[2,3]/gtr_probs['T']) + (R[3,2]/gtr_probs['C']))/2
    gtr_rates['GT'] = ((R[1,3]/gtr_probs['T']) + (R[3,1]/gtr_probs['G']))/2
    
    ag = gtr_rates['AG']
    gtr_rates['AC'] = gtr_rates['AC']/ag
    gtr_rates['AG'] = gtr_rates['AG']/ag
    gtr_rates['AT'] = gtr_rates['AT']/ag
    gtr_rates['CG'] = gtr_rates['CG']/ag
    gtr_rates['CT'] = gtr_rates['CT']/ag
    gtr_rates['GT'] = gtr_rates['GT']/ag
    
    
    a=1
    return gtr_probs,gtr_rates

def read_FASTA(filename):
    '''
    This function reads a FASTA file from a file and returns a dictionary mapping identifiers to sequences
    '''
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

if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your functions, and outputting the results
    '''
    # parse user arguments
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--seqs', required=True, type=str, help="Sequences (FASTA)")
    parser.add_argument('-d', '--distance', required=True, type=str, help="Pairwise Distance")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output GTR Parameters")
    args = parser.parse_args()

    # load sequences
    try:
        seqs = read_FASTA(args.seqs)
    except:
        raise ValueError(ERROR_SEQS)
    if len(seqs) != 2:
        raise ValueError(ERROR_TWO_SEQS)
    r,s = seqs.values()

    # load distance
    from os.path import isfile
    if isfile(args.distance):
        f = open(args.distance); d = float(f.read()); f.close()
    else:
        d = float(args.distance)

    # run student code and output
    gtr_probs, gtr_rates = gtr_params_pair(r,s,d)
    for k in PROB_KEYS:
        assert k in gtr_probs, "Missing GTR stationary probability: %s" % k
        if gtr_probs[k] < 0 or gtr_probs[k] > 1:
            raise ValueError(ERROR_PROB)
    for k in RATE_KEYS:
        assert k in gtr_rates, "Missing GTR transition rate: %s" % k
    if args.output == 'stdout':
        from sys import stdout as outfile
    else:
        outfile = open(args.output,'w')
    outfile.write('%f %f %f %f\n' % tuple(gtr_probs[k] for k in PROB_KEYS))
    outfile.write('%f %f %f %f %f %f\n' % tuple(gtr_rates[k] for k in RATE_KEYS))
    outfile.close()
