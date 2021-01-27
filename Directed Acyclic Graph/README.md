# ECE 208 Homework 1: Warm-Up
In this homework assignment, you will be doing some basic warm-up for this course. This assignment includes programming as well as theoretical components to ensure you are comfortable with both sides of the prerequisites of this course. As with all homework assignments in this class, be sure to follow the instructions carefully to ensure that you receive full credit. All coding for this assignment must be done in **Python 3**.

## Problem 1: Least Common Ancestor
An **unrooted tree** is an undirected graph in which any two nodes are connected by exactly one path. Designating one node as the root and orienting all edges away from the root gives a directed graph that represents a **rooted tree**. More generally, any directed graph with no directed cycles is called a **directed acyclic graph (DAG)**. A rooted tree is called **binary** if the degree of all nodes is either one or three, except the root, which has a degree of two; otherwise, it is called **multifurcating**. A node *v* is a **descendant** of a node *u* if and only if (iff) there is a path from *u* to *v*. A node *u* is an **ancestor** of a node *v* iff *v* is a descendant of *u*. Note that, because it is valid to have a single-node path with no edges, every node *u* is considered its own ancestor. A **common ancestor** of two nodes *u* and *v* in a rooted tree or in a DAG is a node that has both *u* and *v* as descendants. A node *x* is called a **least common ancestor (LCA)** of *u* and *v* in a tree or a DAG iff *x* is a common ancestor of *u* and *v* and no other common ancestor of *u* and *v* is a descendant of *x*. Note that two nodes have only one LCA in a tree but can have multiple LCA in a general DAG. Also note that the LCA of a node *u* and itself is itself.

The goal of this problem is to develop a **dynamic programming algorithm** that, given a tree/DAG, produces the LCAs of all pairs of nodes in the tree/DAG.

### Input
The input is a Python [dictionary](https://docs.python.org/3/tutorial/datastructures.html#dictionaries) `p` representing an input DAG (which may simply be a tree). The keys of the dictionary represent nodes of the DAG (`p` will have *n* keys for a DAG with *n*). For a given node `u`, `p[u]` is a Python [list](https://docs.python.org/3/tutorial/datastructures.html#more-on-lists) containing the parents of `u` (for a tree, this list will have only one element). Note that one (for a tree) or more (for a general DAG) node(s) can have an empty list for `p[u]`; these are the root(s). We suggest that you start by solving the case where the input is a tree and then extend it to a general DAG.

### Output
The output is a Python [nested dictionary](https://www.geeksforgeeks.org/python-nested-dictionary/) (i.e., a dictionary of dictionaries) called `LCA` such that the keys of `LCA` are the keys of `p`, the keys of `LCA[u]` are also the keys of `p`, and the value of `LCA[u][v]` should be a Python [set](https://docs.python.org/3/tutorial/datastructures.html#sets) containing all LCAs of `u` and `v`.

### Important Note
Make sure your algorithm is actually dynamic programming. Ask yourself if you are dividing the problem into sub-problems and if you are saving the solutions to sub-problems for future use. If the answer is no, your solution is not considered a dynamic programming. This question has many solutions that are not based on a dynamic programming. Those are not what we are looking for. For example, although the following solution is a recursive algorithm to find LCAs in a tree (not a DAG), it is not a dynamic programming solution because solutions to subproblems are never used in future (i.e., `LCA` is never read from).

```python
def lca(tree):
    LCA = dict()
    for n in tree.nodes_postorder_traverse():
        if n.is_leaf():
            n.parent.leaves = [n]
        else:
            for a in n.left_child.leaves:
                for b in n.right_child.leaves:
                    LCA[a][b] = n
                    LCA[b][a] = n
        n.leaves = n.left_child.leaves + n.right_child.leaves
    return LCA
```
Thus, this simple and elegant algorithm is not acceptable for this homework.

Also note that we do not require you to find the most efficient solution to this problem. As long as your solution is polynomial in time and space complexity, it will be sufficient for this assignment. We will be testing your code on DAGs with 10 ≤ *n* ≤ 100 nodes.

## Homework Deliverables
* **[writeup.pdf](writeup.pdf):** A write-up that includes the following pieces of information:
    1. The recursive formula of the dynamic programming
    2. The base case of the dynamic programming
    3. The asymptotic running time of the algorithm
    * Note that, for each of the above three items, you must not only provide the correct answer, but you must derive or otherwise explain how you obtained that answer

* **[LCA.py](LCA.py):** Your Python 3 code for solving the problem
    * Usage: `python3 LCA.py <parent_dictionary>`
    * We have provided starter code in the file, but you must fill out the `lca` function (labeled with `TODO`)

## Example Input/Output
You can find example input/output in the [examples](examples) folder. The `*.in` files are the input files, and the `*.out` files are their corresponding output files. All files are plaintext and have been beautified.

### Example Input
```python
{
    'A': [],
    'B': [],
    'C': ['A', 'B'],
    'D': ['A', 'B']
}
```

The above dictionary corresponds to the following graph:

![Alt text](https://g.gravizo.com/svg?digraph%20G%20{A->C;A->D;B->C;B->D;})

### Example Output
```python
{
    'A': {
        'A': {'A'},
        'B': set(),
        'C': {'A'},
        'D': {'A'}
    },

    'B': {
        'A': set(),
        'B': {'B'},
        'C': {'B'},
        'D': {'B'}
    },

    'C': {
        'A': {'A'},
        'B': {'B'},
        'C': {'C'},
        'D': {'A', 'B'}
    },

    'D': {
        'A': {'A'},
        'B': {'B'},
        'C': {'A', 'B'},
        'D': {'D'}
    }
}
```

**Note:** The actual output will be on a single line and will thus be much harder to read. You may want to use a [beautifier](https://codebeautify.org/python-formatter-beautifier) to make it more human-readable.

## Grade Breakdown (100 Points)
* **[writeup.pdf](writeup.pdf): 40 Points**
    * *Correct recursive formula:* 20 points
    * *Correct base case:* 10 points
    * *Correct asymptotic running time:* 10 points

* **[LCA.py](LCA.py): 60 Points**
    * *Binary Tree (small):* 10 points
    * *Binary Tree (large):* 10 points
    * *Multifurcating Tree (small):* 10 points
    * *Multifurcating Tree (large):* 10 points
    * *DAG (small):* 10 points
    * *DAG (large):* 10 points
