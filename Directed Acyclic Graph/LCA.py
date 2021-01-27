#!/usr/bin/env python3

def find_LCAs(parent):
    LCA = dict() # This is the nested dictionary you will be modifying
    # TODO You will fill in the following "lca" function
    def lca(u, v):

        #This function computes the Least Common Ancestor between nodes u and v
        #:param u: node u
        #:param v: node v
        #:return: A set containing the LCAs of u and v
        
        #### return true if a is an ancestor of b ####
        def isAncestor(a,b): # is a an ancestor of b?
            if a == b:
                return True
            for i in parent[b]:
                if isAncestor(a,i):
                    return True
            return False
        ##################################################
        
        #### Remove Common ancestors that aren't LCA ####
        def removeNons(a,b):
            toRemove = set()
            for i in LCA[a][b]:
                for j in LCA[a][b]:
                    if i != j:
                        if isAncestor(i,j) is True:
                            toRemove.add(i)
            for l in toRemove:
                LCA[a][b].remove(l)
        ##################################################
        
        
        if u in LCA.keys():
            if v in LCA[u].keys():
                return #LCA[u][v]
            else:
                LCA[u].update({v: set()}) #.update({v:[]})
        else: # else we create dict for u with dict of list for v inside it
            LCA[u] = {v:set()}
        
        if isAncestor(u,v):
            LCA[u][v].add(u)
        
        for j in parent[u]:
            lca(j,v)
            LCA[u][v].update(LCA[j][v])
        removeNons(u,v)
        
    # Now, we will call your recursive "lca" function on all pairs of nodes to populate the "LCA" dictionary
    for u in parent:
        for v in parent:
            lca(u,v)
            
    return LCA



if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your function, and printing the output
    '''
    from sys import argv
    if len(argv) != 2:
        print("USAGE: %s <parent_dictionary>" % argv[0]); exit(1)
    from ast import literal_eval
    parent = literal_eval(open(argv[1]).read())
    print(find_LCAs(parent))

    
    '''
    binaryDict = {'AA': [], 'BB': ['AA', 'CC'], 'CC':['AA'], 'DD':['AA', 'BB', 'CC'], 'F13':['DD','AA', 'BB', 'CC']}
    #ancest = binaryDict
    #ancest['AA'].append(ancest['BB'][0])
    #print(find_LCAs(binaryDict))
    parent = binaryDict
    print(find_LCAs(parent))
    '''
