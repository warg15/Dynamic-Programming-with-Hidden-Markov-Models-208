# ECE 208 Homework 6: Computing the Likelihood of Trees
In this homework assignment, we will be implementing [Felsenstein's tree-pruning algorithm](https://en.wikipedia.org/wiki/Felsenstein%27s_tree-pruning_algorithm)<sup>1</sup> to compute the (log-)likelihood of observing DNA sequences *S* under the [General Time Reversible (GTR) model](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#GTR_model_(Tavar%C3%A9_1986))<sup>2</sup> given a tree *T* with known branch lengths and topology. This is a Dynamic Programming algorithm that makes computation of the likelihood tractable.

## Problem
The goal of this problem is to implement Felsenstein's tree-pruning algorithm. Given a multiple sequence alignment ***S***, a phylogenetic tree ***T*** (with branch lengths), stationary probabilities **(*π*<sub>A</sub>, *π*<sub>C</sub>, *π*<sub>G</sub>, *π*<sub>T</sub>)**, and transition rates ***α* = *r*(A→G) = *r*(G→A)**, ***β* = *r*(A→C) = *r*(C→A)**, ***γ* = *r*(A→T) = *r*(T→A)**, ***δ* = *r*(G→C) = *r*(C→G)**, ***ε* = *r*(G→T) = *r*(T→G)**, and ***η* = *r*(C→T) = *r*(T→C)**, you will compute the log-likelihood of *T*.

### Input
[compute_likelihood.py](compute_likelihood.py) will take as input the following items:
* `-t`: A phylogenetic tree *T* (either binary or multifurcating) in the [Newick format](https://en.wikipedia.org/wiki/Newick_format)
* `-p`: A file containing the GTR parameters (see below)
* `-s`: The path to a [FASTA file](https://en.wikipedia.org/wiki/FASTA_format) containing a multiple sequence alignment *S* such that each sequence in *S* corresponds to a leaf in *T*
* `-o`: The output file, which will contain the resulting log-likelihood score (default: standard output)

You will only be modifying the `likelihood(tree, seqs, gtr_probs, gtr_rates)` function, which takes as input the following parameters:
* `tree`: A rooted [TreeSwift](https://github.com/niemasd/TreeSwift) `Tree` object (either binary or multifurcating)
* `seqs`: A dictionary where keys are labels corresponding to the labels of the leaves in `tree` and values are sequences (strings)
* `gtr_probs`: A list of floats in the format `[prob_A, prob_C, prob_G, prob_T]` representing the stationary probabilities (sums to 1)
* `gtr_rates`: A list of floats in the format `[rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG]` representing the transition rates

#### GTR Parameter File
The GTR parameter file will be in the following format (an example can be found in [example/gtr_params.txt](example/gtr_params.txt):

```
prob_A prob_C prob_G prob_T
rate_CT rate_AT rate_GT rate_AC rate_CG rate_AG
```

**Note:** Recall that units of GTR rate matrix and units of branch length that you are scoring should be the same, or else the calculations will be wrong. One important note about the input GTR parameter file is that rates in this file **are NOT necessarily normalized to be in substitution units**. Thus, **normalization is required** and is left to your code. In the write-up, you will also need to describe how you normalize the GTR transition rates.

### Output
[compute_likelihood.py](compute_likelihood.py) (and thus `likelihood(tree, seqs, gtr_probs, gtr_rates)`) will output the resulting log-likelihood score. Your output should be a natural logarithm (i.e., base *e*).

### Examples
We have provided 10 multiple sequence alignments along with the phylogenetic tree from which they were all generated, and we include the true log-likelihood scores for each in the [example](example) folder.

## Homework Deliverables
* **[compute_likelihood.py](compute_likelihood.py):** Your Python 3 code for implementing Felsenstein's tree-pruning algorithm
    * Usage: `python3 compute_likelihood.py -t <tree> -p <gtr_params> -s <seqs> -o <out_file>`
    * We have provided starter code in the file, but you must fill out the `likelihood(tree, seqs, gtr_probs, gtr_rates)` function
        * **Note:** Feel free to add import statements to NumPy and SciPy within this function as you need!
* **[writeup.pdf](writeup.pdf):** A write-up that includes the following pieces of information:
    1. The recursive formula of the dynamic programming
    2. The base case of the dynamic programming
    3. The asymptotic running time of the algorithm

## Grade Breakdown (100 Points)
* **[writeup.pdf](writeup.pdf): 20 Points**
    * *Correct recursive formula:* 10 points
    * *Correct base case:* 5 points
    * *Correct asymptotic running time:* 5 points

* **[compute_likelihood.py](compute_likelihood.py): 80 Points**
    * *Log-Likelihood Scores:* 80 points
        * 20 replicates, each worth 4 points

## References
1. Felsenstein J. Maximum Likelihood and Minimum-Steps Methods for Estimating Evolutionary Trees from Data on Discrete Characters. *Systematic Biology*, 22(3):240–249, 1973.
2. Tavaré S. Some probabilistic and statistical problems in the analysis of DNA sequences. *Lectures on Mathematics in the Life Sciences*, 17:57–86, 1986.
