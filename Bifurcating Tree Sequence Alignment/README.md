# ECE 208 Homework 3: Tree-based Multiple Sequence Alignment 

## Part 1: Working with Trees: Compute Tree Diameter
As a warm-up to help you get used to algorithms in trees, we first consider a classical problem of computing tree diameter. In an edge-weighted rooted tree, a diameter is a longest path between any of its two nodes, where the length of a path is defined as the sum of all of its edges. Design an algorithm to compute a diameter of an edge-weighted rooted tree in linear-time.

### Input: 
A rooted tree T with positive branch lengths.

### Output:
A diameter path of T and its length.

You can use any programming language of your choice to complete this assignment, although Python3 is recommended. If you choose to use Python3, a starting code ```compute_diameter.py``` has been provided to you which will handle input/output, as well as a skeleton of the function ```compute_diameter```. The starter code uses [TreeSwift](https://niemasd.github.io/TreeSwift/), a Python library for manipulating trees, so you may find it important to familiarize yourself with the package. If instead you want to use another programming language, you have to modify the file ```compute_diameter.py``` to make it a binder to your code. You have the freedom to modify the given files however you like. However, you **MUST** maintain the sample input/output convention: the only input is an edge-weighted rooted tree in [Newick format](https://en.wikipedia.org/wiki/Newick_format) and the output file consists of two lines where the first line is a diameter path of the tree and the second line is the length of the path. An example of input/output is given in folder ```Example_diameter```.

```
python compute_diameter.py [INPUT_TREE.newick] [OUTPUT_FILE]
```

## Part 2: Multiple Sequence Alignment using a guidance tree
You have learned in the lecture that the complexity of dynamic programming (DP) for pairwise and threeway sequence alignment are O(n<sup>2</sup>) and O(n<sup>3</sup>), respectively, where n is maximum length of the input sequences. Generalizing this approach, the cost to align k sequences is O(n<sup>k</sup>), which means the algorithm's complexity grows exponentially with the numer of species. Despite the intensive computation, multiple sequence alignment (MSA) has a large number of applications in bioinformatics. Therefore, MSA computation is always an active research in the field. In this part of the assignment, you will design an algorithm to align multiple sequences (e.g. up to 50 species) given a guidance tree. The topology of guidance tree reflects the evolutionary relationship of the sequences and the length of each branch is the expected number of substitutions occured on that branch. Using the guidance tree, you need to come up with an algorithm that can infer an alignment more accurately and more efficiently than the DP algorithm. You can assume that the guidance tree is a binary rooted that has non-negative branch lengths.
As in the previous homework, you are again given the BLOSUM62 score matrix, the homologous probabilies, and the indel-rate. You will implement your algorithm in a programming language of your choice and analyze the complexity and accuracy of your algorithm on a simulated dataset. 

### Input:
A rooted binary tree, a set of taxon - sequence pairs, where each taxon is mapped exactly to one sequence, and an indel rate.

### Output:
An MSA of the input sequences.

Similar to part 1, if you choose Python3, the file ```MSA.py``` provides starter code to handle inputs/outputs. If you use any other language, you have to modify ```MSA.py``` to call your code. You **MUST** keep the input/output structure of ```MSA.py``` as follow:

```
python MSA.py -i [INPUT_FILE.fas] -t [TREE_FILE.newick] -o [OUTPUT_FILE.fas] -r [INDEL_RATE]
```
Both input and output sequence files must be in FASTA format, and the tree must be in Newick format. You may find the threeway alignment algorithm in homework 2 applicable here. You are allowed, but are not required, to use the ```threway_align.py``` provided to you or your own implementation from homework 2.

### Testing your code
Your program will be tested on 10 simulated data sets with 1 guidance tree. In folder ```DATA_MSA```, you will find ```guidance.nwk```, ```guidance_rooted.nwk```, and 10 subfolders ```REP_00``` to ```REP_09```, inside each there is a ```input.fas```. You need to run your program on each ```input.fas``` for each replicate using the same ```guidance_rooted.nwk```, and produce an alignment named ```output.fas``` for each. We provide the true alignment (i.e. ```TRUE.fas```) for ```REP_00``` that you can use to compute the accuracy of your algorithm.

## Part 3: Emperical analyses
Once you have your algorithm implemented, you are now ready to analyze the accuracy of your algorithm. You will use the tree ```guidance_rooted.nwk```, the true alignment ```TRUE.fas```, and the output of your algorithm ```output.fas``` for ```REP_00``` to do the following 2 experiments

### Running time
Subsample the set of the given species in the ```input.fas```, restrict the guidance tree to the same set, and run your algorithm for each sample. Make a plot showing the running time versus the number of species. Make at least 10 replicates for each sample size. Is the plot consistent with the theoretical complexity of your algorithm?

### Effect of the diameter on the accuracy
Run your algorithm on the full set. Then create a large number of subsets of size 3. For each subset, restrict your alignment and the guidance tree to that set. Make a plot showing the accuracy of your alignment (using either SP-score or Modeler computed by FastSP) versus the length of the restricted tree's diameter. In addition, run the DP threeway alignment on each subset and compare the accuracy with your alignment. You can show both of them in the same plot or make separate plots for each of them. Do you see any effect of the diameter on alignment accuracy? Do the results match your expectation? Describe what you observe and explain. 

## Deliverables:
* [compute_diameter.py](compute_diameter.py): implementation of part 1
* The 10 output files, ```output.fas```, for each replicate under [DATA_MSA](DATA_MSA)
* [MSA.py](MSA.py): implementation of part 2
* [writeup.pdf](writeup.pdf): description, proof of correctness, and complexity analyses of the algorithms in part 1 and part 2, and all emperical analyses in part 3.

## Grade Breakdown (100 Points)
* **Part 1: 25 Points**
    * Correctness of [compute_diameter.py](compute_diameter.py): 15 Points
    * Complexity analysis: 10 Points
      * *Proof of correctness:* 5 Points
      * *Proof of linear-time complexity:* 5 Points
      
To get full credit, your algorithm must run in linear-time. Higher-order complexity may get partial credits on a case-by-case basis.        
* **Part 2: 50 Points**
   * Accuracy (25 Points): the accuracy of the 10 MSAs, ```output.fas```, that you will produce. Your grade will be based on the relative accuracy of your algorithm comparing to your classmates. 
   * Runtime (10 Points): the runtime of your algorithm will be measured and compared to your classmates.   
   * Algorithm description and analysis (15 Points): you must describe your algorithm to compute the MSA and analyze the time complexity with respect to n and k, where n is the (largest) sequence length and k is the number of species.
   
* **Part 3: 25 Points**
   * Emperical runtime with varying species number: 5 Points
   * Accuracy of tree-based alignment with changing tree diamter: 10 Points
   * Comparison of tree-based and DP alignment with changing tree diameter: 10 Points
   
## Restrictions on the usage of external software/library
* For part 1, you are required to implement the algorithm by your own. You can use any library or software package that provides functionality and API for working with trees and graphs. However, no library functionality other than Newick reading, Newick writing, and tree traversal, is accepted for this part.
* For part 2, you can use any of the built-in functions in TreeSwift. If you do not use TreeSwift, you are only allowed to use the library functions that TreeSwift has an equivalent counterpart. 
* For part 3, you have the freedom to use any library or software to facilitate the data analysis process. For instance, to compute the tree diameter, you can use your implementation in part 1, the Treeswift's built-in function, or any other software/library. The same applied for any other parts of your pipeline. You do not need to turn in the code for part 3. All results and figures must be included in the writeup.
* Although we allow the usage of a programming language other than Python3, you are totally responsible for the deliverable of your code. You should think of your submission as a mini Software that is ready to be downloaded and run. All parts of your program must be compiled, built, and included in the submission. Your code must be wrapped inside the two files ```compute_diameter.py``` and ```MSA.py``` so that the user can execute your program using only these two files. You are also responsible for the differences in OS system and make sure that your program can run on multiple OS. In principle, we should be able to test your work by just calling ```compute_diameter.py``` and ```MSA.py```. If extra installation is needed, you have to include the detailed instructions on ```install.txt```. Note that, *you may lose partial or all points for any part of the assignment that needs extra installations*, so plan accordingly. Thus, if you don't have experience with software development and distribution, we highly recommend using Python3 and TreeSwift for this assignment.   

## Final Remarks:
All the MSAs that you turn in must be the actual outputs of your program. We will perform random tests and if your results cannot be reproduced, your case will be considered a violation of Academic Integrity. If you implement a randomized algorithm, you **MUST** use a **fixed** seed number and include it in ```MSA.py```. 
