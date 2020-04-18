# codon-opt

## Codon optimization problem with motif engineering

**Given:**
A target protein *P* (i.e., a sequence of amino acids), a set of desired motifs *D* and a set of undesired motifs *U*

**Find:** A sequence *S* of basis triplets such that:
1. The *k*-triplet of *S* encodes the *k*-th amino acid of *P*
2. *S* contains the minimum number of substrings *u* of *U*
3. *S* contains the maximum number of substrings *d* of *D*
4. The *Codon Adaptation Index* (CAI) of *S* is maximized

This code constructs and solves an ILP model 
by the commercial ILP solver GUROBI (www.gurobi.com)

Further details can be found on

```
@article{ArbibPinarRossiTessitore2020,
author = {Arbib, C., Pinar, M., Rossi, F., Tessitore A.},
title = {Codon Optimization by 0-1 Linear Programming},
journal = {Computers and Operations Research},
volume = {119},
number = {},
pages = {},
year = {2020},
doi = {https://doi.org/10.1016/j.cor.2020.104932},
}
```

## Input data

A ```.fasta``` file containing a list of target proteins
(see https://en.wikipedia.org/wiki/FASTA_format for format description
or the sample file ```target.fasta```).

Two (optional) text files containing the lists of forbidden and desired motifs
(see, for instance, ```forbidden.txt``` and ```desired.txt``` files).

## Requirements

python=3.7.6
pandas=0.25
gurobipy=8.1.1

## Usage

```python codon-opt.py target.fasta [-f forbidden_file.txt]  [-d desired_file.txt] [-o output_file]```






