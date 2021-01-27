# ECE 208 Homework 2: Three-Way Sequence Alignment
In this programming assignment, we extend the pairwise sequence alignment algorithm to three sequences. You will be given three amino acid sequences *s*<sub>1</sub>, *s*<sub>2</sub>, and *s*<sub>3</sub> with sequences lengths *m*<sub>1</sub>, *m*<sub>2</sub>, and *m*<sub>3</sub>, respectively. The sequences are homologous, but have substitutions, insertions, and deletions, and therefore, differ in lengths.

## Problem
Your goal is to design a **dynamic programming algorithm** to align 3 sequences *s*<sub>1</sub>, *s*<sub>2</sub>, and *s*<sub>3</sub> such that **the sum of pairwise alignment scores** is maximized. 

## Sequence Evolution Model
You will use **BLOSUM62**<sup>1</sup> as the substitution model, which will be given to you. To make the problem more interesting and also, closer to the reality, the **gap penalty** will be hiden from you. Instead, we give you the **indel rate**, from which you can compute **gap penalty** using some logic and reasoning.

### Indel Rate
**Indel rate** is the probability of an indel relative to the probability of a substitution. For example, an indel rate of 0.1 means that, on average, for every 10 substitutions, you would expect 1 insertion or deletion. The indel rate will be given to you as input, and you need to find the best settings of gap penalty with respect to the given indel rate and the **BLOSUM62** model.

### Affine versus non-affine gap penalty
In reality, a **single indel event** can involve one or more characters; that is, a single event can insert or delete a sequence of length l>=1 of consecutive characters to the ancestor sequence. In such a setting, an **affine** gap penalty is usually more suitable than the **non-affine** counterpart to align the sequences. 
However, for the sake of simplicity, you can choose to implement either **non-affine** gap penalty as taught in class, or the **affine** version which requires some updates in the algorithm. The choice is yours. Note that,
finding the best way to penalize gaps is a part of this assignment, and there is no single correct solution for this.

## Hints
If you compare this problem to the pairwise alignment problem, you will see that instead of a 2D matrix as in the pairwise case, you will need to compute a 3D matrix of dimensions (|*s*<sub>1</sub>|+1) × (|*s*<sub>2</sub>|+1) × (|*s*<sub>3</sub>|+1). The recursive formula will also be similar, just with more cases involved.

When you want to find gap penalty, you may find that you want to know the *p*<sub>*ij*</sub> values used in the BLOSUM62 matrix, which we have included in [blosum62sym.csv](blosum62sym.csv). The other thing you may need to know is how this is turned into the BLOSUM62 matrix. Convince yourself that what they used is the following (important information here is that log is in base-2 and there is a multiplication by 2):

![equation](https://latex.codecogs.com/svg.latex?2%5Clog_%7B2%7D%5Cleft%28%5Cfrac%7Bp_%7Bij%7D%7D%7Bq_%7Bi%7Dq_%7Bj%7D%7D%5Cright%20%29)

## Code Structure
A python file `threeway_align.py` is given to you as the **starting** code. Type `python threeway_align.py -h` to see the input/output options. We will use the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) for both inputs and outputs. For your convenience, we provided the code for reading and writing FASTA files. Your task is to write the function `threeway_align` inside file `threeway_align.py`.

### Input
The `threeway_align` function takes as input the following parameters:
* Three Python [strings](https://docs.python.org/3/library/stdtypes.html#textseq) `s1`, `s2`, and `s3`, representing the three amino acid sequences *s*<sub>1</sub>, *s*<sub>2</sub>, and *s*<sub>3</sub>
* A positive real-valued number `indel_rate`, which is the ratio of the indel rate to an average substitution rate of 1

### Output
The `threeway_align` function returns three strings `aligned_s1`, `aligned_s2`, and `aligned_s3` with the following properties:
* All three have equal length
    * `len(aligned_s1) == len(aligned_s2) == len(aligned_s3)`
* Removing all gap characters from `aligned_s1`, `aligned_s2`, and `aligned_s3` yields `s1`, `s2`, and `s3`, respectively
    * `aligned_s1.replace('-','') == s1 and aligned_s2.replace('-','') == s2 and aligned_s3.replace('-','') == s3`

### Important Note
You are highly adviced **against** changing any part of the starting code except for the body of the `threeway_align` function. 
If you do, you are totally responsible for the behavior of the entire file `threeway_align.py`. You **must** make sure that the file `threeway_align.py` can handle FASTA files and have the same input/output settings as it was originally given to you. This is necessary for us to test your code.

## Testing and Evaluation
After you implemented the `threeway_align` function, you would be able to call `threeway_align.py` as follow:

`python threeway_align.py -i <input_file>.fas -o <output_file>.fas`

In the `Examples` folder, there are three tests given to you, each of which includes a `.fas` file - the unaligned sequences, and a `.aln` file - the true alignment. 
These tests, as well as all the tests that your program will be graded on, were simulated following a model that is **hidden** to you. You know that substitutions follow BLOSUM62, but you do not know much about the indel model other than its indel rate. Therefore, you are not expected to have a program that gives the **perfect** alignments for all the tests. To evaluate the accuracy of your program, we will compare the alignment that your program produces with the true alignment. To do that, we use **FastSP**<sup>2</sup> to compute the SP-score. You are encouraged to read the paper to understand what SP-score is and how FastSP works.

We include the file `FastSP.jar` that you can run to compute the SP-score of your alignment (in addition to other accuracy measurements that will not be used here). If you have java installed, you can use the following command to run FastSP:

`java -jar FastSP.jar -r Examples/Test1/01.aln -e <YOUR_ALIGNMENT.aln>`

SP-score is a similarity measurement of two alignments. It takes values from 0 to 1, and the higher the SP-score, the more similar the two alignments are to each other. As mentioned before, you are not expected to produce perfect alignments. A good program will give alignments with high SP-scores, but not neccessarily obtain the perfect SP-score of 1 for all tests. To give you a sense of what a good SP-score can be, we include the `.sp` file for each of the three tests. Although you are not required to run FastSP, you are encouraged to do so and compare your SP-scores with those of the provided `.sp` files. In all the provided tests, the indel rate is 0.01.

## Homework Deliverables
* **[writeup.pdf](writeup.pdf):** A write-up that includes the following pieces of information:
    1. The recursive formula of the dynamic programming
    2. The base case(s) of the dynamic programming
    3. The asymptotic running time of the algorithm
    4. Your logic behind setting the gap penalty
    * For each of the above items, you must derive or explain how you obtained the answer.

* **[threeway_align.py](threeway_align.py):** Your Python3 code for solving the problem
    * We have provided starter code in the file, but you must fill out the `threeway_align` function (labeled with `TODO`)
    * Usage: `python3 threeway_align.py -i <input_sequences> -r <indel_rate> -o <output_alignment>`
        * Both `<input_sequences>` and `<output_alignment>` will be in the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format)

## Grade Breakdown (100 Points)
* **[writeup.pdf](writeup.pdf): 40 Points**
    * *Correct recursive formula:* 10 points
    * *Correct base case(s):* 10 points
    * *Correct asymptotic running time:* 10 points
    * *Valid logic behind setting the gap penalty:* 10 points

* **[threeway_align.py](threeway_align.py): 60 Points**

    You will be graded based on the accuracy and execution time of your programming comparing to other classmates.
    
    * *Accuracy of the program:* 40 points
    * *Execution time:* 20 points
    
## Final Remarks
You are encouraged to think out-of-the-box and search the literatures for or come up with other methods to do alignment other than our "traditional" dynamic programming algorithm. **However**, as the goal of this homework is for you to apply dynamic programming to sequence alignment, you **MUST** implement a dynamic programming algorithm to **maximize the sum of pairwise alignment scores**. Thus, you **MUST** fulfill all components listed in the "Homework Deliverables" section, including the formulation of the dynamic programming algorithm (base cases and recursive formula) and the derivation of the gap penalty that you use, to get full credit for this assignment.

## References
1. Henikoff S, Henikoff JG. Amino acid substitution matrices from protein blocks. *Proceedings of the National Academy of Sciences*, 89(22):10915–10919, 1992. [doi:10.1073/pnas.89.22.10915](https://doi.org/10.1073/pnas.89.22.10915)
2. Mirarab S, Warnow T. FASTSP: linear time calculation of alignment accuracy. *Bioinformatics*, 27(23):3250–3258, 2011. [doi:10.1093/bioinformatics/btr553](https://doi.org/10.1093/bioinformatics/btr553)
