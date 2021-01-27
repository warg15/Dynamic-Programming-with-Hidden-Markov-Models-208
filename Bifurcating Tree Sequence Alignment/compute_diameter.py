#!/usr/bin/env python3

def compute_diameter(tree): # TODO You will fill in this function
    '''
    This function computes a diameter path (i.e. a longest path between any two nodes) of a tree and its length
    :param tree: the tree
    :return: a list d_path containing the diameter path and a floating point number d_length storing the diameter length
    '''
        
    def gV (lst):
        value = 0
        for owo in lst:
            value = value + lengths[names.index(owo)]
        return value
    
    #lst = ['a','i','e','d']
    #b = gV(lst)
    
    lengths = list()
    names = list()
    
    whichNode = tree.label_to_node(selection='all')
    
    for node in tree.traverse_postorder():
        #print(node)
        names.append( node.get_label())
        if node.is_root():
            print('root')
            lengths.append(0)
        else:
            lengths.append( node.get_edge_length())
    
    ancestors = {}
    values = {}
    curr_Ancest = list()
    curr_value = list()
    
    #rootNode
    
    for node in tree.traverse_preorder(): #find the ancestors and the traverse length 
                                            #from them to root
        if node.is_root():
            rootNode = node
        label = node.get_label()
        if node.is_root() == False:
            parent = node.get_parent()
            parent_label = parent.get_label()
            if len(curr_Ancest) > 0:
                while curr_Ancest != [] and curr_Ancest[0] != parent_label:
                    curr_Ancest.pop(0)
                    curr_value.pop(0)
        index = names.index(label)
        curr_value.insert(0, lengths[index]) 
        curr_Ancest.insert(0, label)
        
        values[label] = curr_value[:]
        ancestors[label] = curr_Ancest[:]
        
        
    start_Node = list(ancestors.keys())[0]
    start_length = sum(values[start_Node])
    
    for node in tree.traverse_preorder(): #find the node to start from
        challenger = node.get_label()
        challenger_length = sum(values[challenger])
        if challenger_length > start_length:
            start_Node = challenger
            start_length = challenger_length
    d_path = []
    d_length = 0
    for node in tree.traverse_preorder(): #find diameter
        curr_node = node.get_label()
        list_curr = ancestors[curr_node][:]
        list_start = ancestors[start_Node][:]
        intersection = list(set(list_curr) & set(list_start))
        list_curr = list(set(list_curr) - set(intersection))
        list_start = list(set(list_start) - set(intersection))
        list_curr.reverse()
        path = list_start + list_curr
        length = gV(path)
        
        if length > d_length:
            d_length = length
            d_path = path[:]
        
    #order the path
    curr_node = whichNode[start_Node]
    good_node = whichNode[start_Node]
    d_path.append(rootNode.get_label())
    copy_path = d_path[:]
    new_path = []
    new_path.append(start_Node)
    copy_path.remove(start_Node)
    for j in range(len(copy_path)):
        dist = 999999999
        for k in range(len(copy_path)):
            testNode = whichNode[copy_path[k]]
            if testNode.get_parent() == curr_node:
                new_path.append(testNode.get_label())
                good_node = testNode
            elif curr_node.get_parent() == testNode:
                new_path.append(testNode.get_label())
                good_node = testNode
        copy_path.remove(good_node.get_label())
        curr_node = good_node
    
    d_length = length
    d_path = new_path[:]
    
    #for node in tree.traverse_postorder():
        # TODO: replace with your code!
     #   raise Exception("NOT IMPLEMENTED ERROR!")
    print(d_path)
    print(d_length)
    return d_path,d_length

if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your function, and printing the output
    '''
    from sys import argv
    from treeswift import read_tree_newick
    
    if len(argv) != 3:
        print("USAGE: %s <newick_tree> <output_file>" % argv[0]); exit(1)
    
    tree = read_tree_newick(argv[1])
    outfile = argv[2]
    d_path,d_length = compute_diameter(tree)
    
    with open(outfile,'w') as fout:
        for x in d_path:
            fout.write(x + " ")
        fout.write("\n" + str(d_length))    
