#!/usr/bin/env python3
ERROR_PROB = "Stationary probabilities must be between 0 and 1"
ERROR_SEQS = "Invalid sequence file"
PROB_KEYS = ['A', 'C', 'G', 'T']
RATE_KEYS = ['CT', 'AT', 'GT', 'AC', 'CG', 'AG']

import numpy as np
import treeswift 
import scipy
from scipy.linalg import expm, sinm, cosm, logm



def gtr_params(tree, seqs):
    
    def gtr_params_pair(r, s, d):
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
            transCount[keyString] = transCount[keyString] + 1
            keyString = sChar + rChar
            transCount[keyString] = transCount[keyString] + 1
            
        newerR = np.array([[transCount['AA']/amountLetters[0], transCount['AG']/amountLetters[0], transCount['AC']/amountLetters[0], transCount['AT']/amountLetters[0] ],
                      [transCount['GA']/amountLetters[1], transCount['GG']/amountLetters[1], transCount['GC']/amountLetters[1], transCount['GT']/amountLetters[1] ],
                      [transCount['CA']/amountLetters[2], transCount['CG']/amountLetters[2], transCount['CC']/amountLetters[2], transCount['CT']/amountLetters[2] ],
                      [transCount['TA']/amountLetters[3], transCount['TG']/amountLetters[3], transCount['TC']/amountLetters[3], transCount['TT']/amountLetters[3] ]])
        
        R = logm(newerR)/t
        
        gtr_rates['AC'] = ((R[0,2]/gtr_probs['C']) + (R[2,0]/gtr_probs['A']))/2
        gtr_rates['AG'] = ((R[0,1]/gtr_probs['G']) + (R[1,0]/gtr_probs['A']))/2
        gtr_rates['AT'] = ((R[0,3]/gtr_probs['T']) + (R[3,0]/gtr_probs['A']))/2
        gtr_rates['CG'] = ((R[2,1]/gtr_probs['G']) + (R[1,2]/gtr_probs['C']))/2
        gtr_rates['CT'] = ((R[2,3]/gtr_probs['T']) + (R[3,2]/gtr_probs['C']))/2
        gtr_rates['GT'] = ((R[1,3]/gtr_probs['T']) + (R[3,1]/gtr_probs['G']))/2
        
        return gtr_probs,gtr_rates
    
    
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
            children = node.child_nodes()
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
                    L1 = np.transpose(Q1[s,:])@mat1[:,i]
                    L2 = np.transpose(Q2[s,:])@mat2[:,i]
                    matrix[s,i] = L1 * L2
                
            matrices[node] = matrix
        
        #REWRITE USING DICTIONARIES
        #prob_A, prob_C, prob_G, prob_T = gtr_probs # You can use these if it's more convenient
        #rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG = gtr_rates # You can use these if it's more convenient
        prob_A, prob_C, prob_G, prob_T = gtr_probs['A'], gtr_probs['C'], gtr_probs['G'], gtr_probs['T']
        rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG = gtr_rates['CT'], gtr_rates['AT'], gtr_rates['GT'], gtr_rates['AC'], gtr_rates['CG'], gtr_rates['AG']
        
        log_likelihood_score = 0. # You can modify this variable, or you can use your own; either way is fine
        # TODO Your code here
        R = np.array([[-(rate_AG*prob_G + rate_AC*prob_C + rate_AT*prob_T), rate_AG*prob_G, rate_AC*prob_C, rate_AT*prob_T],
                      [rate_AG*prob_A, -(rate_AG*prob_A + rate_CG*prob_C + rate_GT*prob_T), rate_CG*prob_C, rate_GT*prob_T],
                      [rate_AC*prob_A, rate_CG*prob_G, -(rate_AC*prob_A + rate_CG*prob_G + rate_CT*prob_T), rate_CT*prob_T],
                      [rate_AT*prob_A, rate_GT*prob_G, rate_CT*prob_C, -(rate_AT*prob_A + rate_GT*prob_G + rate_CT*prob_C)]])
        
        d = -prob_A*R[0,0] -prob_G*R[1,1] -prob_C*R[2,2] -prob_T*R[3,3] #for normalizing the rates matris, R
        l = len(seqs[list(seqs.keys())[0] ])   #length of sequences
        matrices = {}
        
        for node in tree.traverse_postorder():
            if node.is_leaf():popLeafs(node)
            else:
                updateMatrix(node)
                if node.is_root(): root = node
        
        like = 0
        pis = np.array([prob_A, prob_G, prob_C, prob_T])
        rootMatrix = matrices[root]
        for i in range(0,l):
            temp = np.log(pis@rootMatrix[:,i])
            like = like + temp
        #print(like)
        return like
    
    
    '''
    This function estimates GTR model parameters from a tree and sequences
    :param tree: A TreeSwift Tree object representing the phylogenetic tree
    :param seqs: A dictionary where keys are labels corresponding to the labels of the leaves in ``tree`` and values are sequences (strings)
    :return: A dictionary ``gtr_probs`` storing the GTR stationary probabilities, and a dictionary ``gtr_rates`` storing the GTR transition rates
    '''
    gtr_probs = dict() # keys: {'A', 'C', 'G', 'T'}   values: the corresponding GTR stationary probabilities
    gtr_rates = dict() # keys: {'AC', 'AG', 'AT', 'CG', 'CT', 'GT'}   values: the corresponding GTR transition rates
    # TODO Your code here
    

    
    aProbs = []
    gProbs = []
    cProbs = []
    tProbs = []
    
    acValues = []
    agValues = []
    atValues = []
    cgValues = []
    ctValues = []
    gtValues = []
    
    for i in tree.traverse_leaves():
        #print(i)
        for j in tree.traverse_leaves():
            iLabel = i.get_label()
            jLabel = j.get_label()
            #if iLabel != jLabel:
            d = tree.distance_between(i, j)
            r = seqs[iLabel]
            s = seqs[jLabel]
            if r != s:
                gtr_probs, gtr_rates = gtr_params_pair(r,s,d)
                aProbs.append(gtr_probs['A'])
                gProbs.append(gtr_probs['G'])
                cProbs.append(gtr_probs['C'])
                tProbs.append(gtr_probs['T'])
                
                acValues.append(gtr_rates['AC'])
                agValues.append(gtr_rates['AG'])
                atValues.append(gtr_rates['AT'])
                cgValues.append(gtr_rates['CG'])
                ctValues.append(gtr_rates['CT'])
                gtValues.append(gtr_rates['GT'])

    
    avgAProbs = sum(aProbs)/len(aProbs)
    avgGProbs = sum(gProbs)/len(gProbs)
    avgCProbs = sum(cProbs)/len(cProbs)
    avgTProbs = sum(tProbs)/len(tProbs)
    
    avgAC = sum(acValues)/len(acValues)
    avgAG = sum(agValues)/len(agValues)
    avgAT = sum(atValues)/len(atValues)
    avgCG = sum(cgValues)/len(cgValues)
    avgCT = sum(ctValues)/len(ctValues)
    avgGT = sum(gtValues)/len(gtValues)
    
    
    gtr_rates['AC'] = avgAC
    gtr_rates['AG'] = avgAG
    gtr_rates['AT'] = avgAT
    gtr_rates['CG'] = avgCG
    gtr_rates['CT'] = avgCT
    gtr_rates['GT'] = avgGT
    
    
    rates = dict()
    gt = gtr_rates['GT']
    rates['AC'] = gtr_rates['AC']/gt
    rates['AG'] = gtr_rates['AG']/gt
    rates['AT'] = gtr_rates['AT']/gt
    rates['CG'] = gtr_rates['CG']/gt
    rates['CT'] = gtr_rates['CT']/gt
    rates['GT'] = gtr_rates['GT']/gt
    print(rates)
    
    gtr_probs['A'] = avgAProbs
    gtr_probs['G'] = avgGProbs
    gtr_probs['C'] = avgCProbs
    gtr_probs['T'] = avgTProbs
    

    ####################################################
    ####################################################
    ###########          OPTIMIZE           ############
    ####################################################
    ####################################################
    
    stepAmount = 0.03
    numSteps = 1
    countSteps = 7

    order = ['AG', 'CT','AC','GT','CG','AT']

    
    rates = dict()
    rates['AC'] = gtr_rates['AC']
    rates['AG'] = gtr_rates['AG']
    rates['AT'] = gtr_rates['AT']
    rates['CG'] = gtr_rates['CG']
    rates['CT'] = gtr_rates['CT']
    rates['GT'] = gtr_rates['GT']
    
    def decreaseSteps(smallerRates):
        smallerRates[j] = smallerRates[j]-stepSize
        newLiklihood = likelihood(tree, seqs, gtr_probs, smallerRates)
        #print(newLiklihood)
        return newLiklihood
    
    def increaseSteps(largerRates):
        largerRates[j] = largerRates[j]+stepSize
        newLiklihood = likelihood(tree, seqs, gtr_probs, largerRates)
        #print(newLiklihood)
        return newLiklihood
    
    breakOut = False
    for i in range(numSteps):
        #print(stepYesNo)
        print(likelihood(tree, seqs, gtr_probs, gtr_rates))
        
        #stepAmount = stepAmount*.8
        noImprove = 0
        if breakOut: break
        #for k in range(6):
        for j in order:
            #j = order[k]
            print(j)
            #if stepYesNo[k] == 1:
            stepSize = stepAmount*gtr_rates[j]
            #stepSize = 0.01
            
            improve = [0,0,0] #smaller, same, larger
            
            improve[1] = likelihood(tree, seqs, gtr_probs, gtr_rates)
            
            smallerRates = gtr_rates.copy()
            smallerRates[j] = smallerRates[j]-stepSize
            improve[0] = likelihood(tree, seqs, gtr_probs, smallerRates)
            
            largerRates = gtr_rates.copy()
            largerRates[j] = largerRates[j]+stepSize
            improve[2] = likelihood(tree, seqs, gtr_probs, largerRates)
            
            maxIndex = improve.index(max(improve))
            iteration = 0

            if maxIndex == 0:
                curr_like = improve[0]
                new_like = decreaseSteps(smallerRates)
                while  curr_like < new_like and iteration < countSteps:
                    curr_like = new_like
                    new_like = decreaseSteps(smallerRates)
                    iteration = iteration +1
                smallerRates[j] = smallerRates[j]+stepSize
                rates[j] = smallerRates[j]
                
            if maxIndex == 2:
                curr_like = improve[2]
                new_like = increaseSteps(largerRates)
                while  curr_like < new_like and iteration < countSteps:
                    curr_like = new_like
                    new_like = increaseSteps(largerRates)
                    teration = iteration +1
                largerRates[j] = largerRates[j]-stepSize
                rates[j] = largerRates[j]
                
            if maxIndex == 1:
                #print('did not improve', j)
                noImprove = noImprove + 1
                if noImprove == 6:
                    breakOut = True
        
        gtr_rates['AC'] = rates['AC']
        gtr_rates['AG'] = rates['AG']
        gtr_rates['AT'] = rates['AT']
        gtr_rates['CG'] = rates['CG']
        gtr_rates['CT'] = rates['CT']
        gtr_rates['GT'] = rates['GT']
            
    
    gt = gtr_rates['GT']
    gtr_rates['AC'] = gtr_rates['AC']/gt
    gtr_rates['AG'] = gtr_rates['AG']/gt
    gtr_rates['AT'] = gtr_rates['AT']/gt
    gtr_rates['CG'] = gtr_rates['CG']/gt
    gtr_rates['CT'] = gtr_rates['CT']/gt
    gtr_rates['GT'] = gtr_rates['GT']/gt
    
    
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
    parser.add_argument('-t', '--tree', required=False, type=str, default='stdin', help="Input Tree File (Newick)")
    parser.add_argument('-s', '--seqs', required=True, type=str, help="Multiple Sequence Alignment (FASTA)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output GTR Parameters")
    args = parser.parse_args()

    # load input tree
    from treeswift import read_tree_newick
    if args.tree == 'stdin':
        from sys import stdin
        tree = read_tree_newick(stdin)
    else:
        tree = read_tree_newick(args.tree)
    leaves = {l.label for l in tree.traverse_leaves()}

    # load sequences
    try:
        seqs = read_FASTA(args.seqs)
    except:
        raise ValueError(ERROR_SEQS)
    for k in seqs:
        assert k in leaves, "Sequence ID not in tree: %s" % k

    # run student code and output
    gtr_probs, gtr_rates = gtr_params(tree,seqs)
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
