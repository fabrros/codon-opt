# codon-opt

**Codon optimization problems with motif engineering (COME)**

**Given:**
A target protein *P* (i.e., a sequence of amino acids), a set of desired motifs *D* and a set of undesired motifs *U*,
two positive integers *d* and *u*

**Find:** A sequence *S* of basis triplets such that:
1. The *k*-triplet of *S* encodes the *k*-th amino acid of *P*
2. *S* contains at most *u* substrings of *U*
3. *S* contains at least *d* substrings of *D*
4. The *Codon Adaptation Index* (CAI) of *S* is maximized

This code solves COME by constructing an ILP model 
and by solving the model 
with the commercial ILP solver GUROBI (www.gurobi.com)

