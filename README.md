# codon-opt

## Codon optimization problems with motif engineering

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
journal = {Computers and Operations Research (submitted)},
volume = {},
number = {},
pages = {},
year = {},
doi = {},
}
```

## Input data

A ```.fasta``` file containing a list of target proteins
(see https://en.wikipedia.org/wiki/FASTA_format for format description)

Two (optional) text files ```forbidden.txt``` and ```desired.txt```
containing a list of forbidden and desired motifs

## Python modules

```python
    import math
    import pandas
    import argparse
    import gurobipy
```

## Usage

```codon-opt target.fasta [-f forbidden_file.txt]  [-d desired_file.txt] ```






