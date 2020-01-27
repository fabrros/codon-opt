import gurobipy as gb
import pandas as pd
import math

#
# Inverse table for standard genetic code
# Amino acid -> DNA codons list
# Stop codons are denoted by '#'
#

InvTableDNA = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'K': ['AAA', 'AAG'],
    'N': ['AAT', 'AAC'],
    'M': ['ATG'],
    'D': ['GAT', 'GAC'],
    'F': ['TTT', 'TTC'],
    'C': ['TGC', 'TGT'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'E': ['GAA', 'GAG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'W': ['TGG'],
    'H': ['CAT', 'CAC'],
    'Y': ['TAT', 'TAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    '#': ['TAA', 'TGA', 'TAG']
}

#
# Reverse InvTableDNA
#

TableDNA = {i: k for k, v in InvTableDNA.items() for i in v}

#
# Codon normalized frequencies
#

FrequencyNorm = {'GGG': 0.48,
                 'GGA': 0.48,
                 'GGT': 1.00,
                 'GGC': 1.00,
                 'GAG': 0.52,
                 'GAA': 1.00,
                 'GAT': 1.00,
                 'GAC': 0.58,
                 'GTG': 1.00,
                 'GTA': 0.54,
                 'GTT': 0.91,
                 'GTC': 0.62,
                 'GCG': 1.00,
                 'GCA': 0.83,
                 'GCT': 0.68,
                 'GCC': 0.90,
                 'AGG': 0.17,
                 'AGA': 0.28,
                 'AGT': 0.76,
                 'AGC': 1.00,
                 'AAG': 0.39,
                 'AAA': 1.00,
                 'AAT': 1.00,
                 'AAC': 0.87,
                 'ATG': 1.00,
                 'ATA': 0.32,
                 'ATT': 1.00,
                 'ATC': 0.76,
                 'ACG': 0.66,
                 'ACA': 0.57,
                 'ACT': 0.55,
                 'ACC': 1.00,
                 'TGG': 1.00,
                 'TGA': 0.55,
                 'TGT': 0.92,
                 'TGC': 1.00,
                 'TAG': 0.15,
                 'TAA': 1.00,
                 'TAT': 1.00,
                 'TAC': 0.64,
                 'TTG': 0.29,
                 'TTA': 0.34,
                 'TTT': 1.00,
                 'TTC': 0.68,
                 'TCG': 0.56,
                 'TCA': 0.72,
                 'TCT': 0.77,
                 'TCC': 0.64,
                 'CGG': 0.39,
                 'CGA': 0.24,
                 'CGT': 1.00,
                 'CGC': 0.96,
                 'CAG': 1.00,
                 'CAA': 0.50,
                 'CAT': 1.00,
                 'CAC': 0.70,
                 'CTG': 1.00,
                 'CTA': 0.11,
                 'CTT': 0.29,
                 'CTC': 0.23,
                 'CCG': 1.00,
                 'CCA': 0.48,
                 'CCT': 0.46,
                 'CCC': 0.32}

#
# Codon Adaptation Index costs to optimize CAI
#

Costs = {k: math.log(FrequencyNorm[k]) for k in FrequencyNorm.keys()}

#
# Functions
#

def read_fasta(fastafile:str)->pd.DataFrame:
    """
    Read fasta file using pandas
    :param: filename
    :return: pandas DataFrame with target proteins
    """
    df = pd.read_csv(fastafile, header=None, comment='>', engine='python',
                     names=['Protein'])

    df['Protein'] = df['Protein'].apply(lambda x: x[:len(x) - len(x) % 3])

    df['Length'] = df['Protein'].str.len()

    print(fastafile, ('has %d protein(s)' % df.shape[0]))

    return df


def read_motifs(forbiddenfile: str,
                desiredfile: str)->(list, list):
    """
    Read forbidden and desired motifs
    :param: forbidden and desired file names
    :return: lists of forbidden and desired motifs
    """

    forbidden = []
    desired = []

    if forbiddenfile is not None:
        with open(forbiddenfile) as file:
            forbidden = file.readlines()

        forbidden = [x.strip() for x in forbidden]

    if desiredfile is not None:
        with open(desiredfile) as file:
            desired = file.readlines()

        desired = [x.strip() for x in desired]

    return forbidden, desired


def find_amino_from_codon(codon: str,
                          table: dict)->dict:

    """
    Find all possible aminos that can be codified by a (partly empty) codon
    :param: a 3-chars string representing a codon. May contain '*' as don't care
    either at the beginning or at the end of the sequence
    :return: a dictionary with all possible amino acids that the sequence codifies
    """

    aminos = {}

    if codon.isalpha():
        if table[codon] in aminos:
            aminos[table[codon]].append(codon)
        else:
            aminos[table[codon]] = [codon]

    if codon[0] == '*':
        partialcodon = codon.lstrip('*')
        offset = len(codon) - len(partialcodon)
        for i in table:
            if partialcodon == i[offset:]:
                if table[i] in aminos:
                    aminos[table[i]].append(i)
                else:
                    aminos[table[i]] = [i]

    if codon[2] == '*':
        partialcodon = codon.rstrip('*')
        offset = len(codon) - len(partialcodon)

        for i in table:
            if partialcodon == i[:3 - offset]:
                if table[i] in aminos:
                    aminos[table[i]].append(i)
                else:
                    aminos[table[i]] = [i]

    return aminos


def find_amino_sequence_in_protein(sequence: str,
                                   aminos: str)->list:
    """
    Check if a non-empty sequence of amino acids is
    contained in the aminos string
    :param: two sequences (strings) of amino acids
    :return: an (eventually empty) list of positions
    """
    feasiblepos = []

    pos = aminos.find(sequence)

    while pos != -1:
        restart = pos
        start = pos
        end = (pos + len(sequence) - 1)
        feasiblepos.append((start, end))
        pos = aminos.find(sequence, restart + 1)

    return feasiblepos


def find_motif_position(motif:str,
                        aminos:str)->(dict,dict):
    """
    Given a sequence of amino acids find positions in which
    a motif may be located
    :param:  motif, amino acids sequense
    :return: a dictionaries with feasible positions and
    corresponding codon encoding
    Positions starts from 0 and refer to basis position
    Codon encoding at first and last level may contain
    more than one alternative
    """

    feasiblePositions = {}
    feasiblePositionsBasis = {}

    #
    # Test three possible positions
    #

    for offset in range(3):

        sequence = ''.join('*' * offset) + motif
        startoffset = offset

        if len(sequence) % 3 != 0:
            endoffset = (3 - (len(sequence) % 3))
            sequence += ('*' * (3 - (len(sequence) % 3)))
        else:
            endoffset = 0

        firstlevel = find_amino_from_codon(sequence[0:3], TableDNA)

        if len(sequence) > 3:
            lastlevel = find_amino_from_codon(sequence[-3:], TableDNA)
        else:
            raise ValueError('Motif must have length > 3')

        middlelevels = []

        for idx in range(3, len(sequence) - 3, 3):
            codon = sequence[idx: idx + 3]
            middlelevels.append(find_amino_from_codon(codon, TableDNA))

        #
        # Build the sequences starting from possible
        # encodings
        #

        for head in firstlevel:
            sequence = [{head: firstlevel[head]}]

            sequence += middlelevels

            aminosequence = ''.join([[*i.keys()][0] for i in middlelevels])
            aminosequence = head + aminosequence

            for tail in lastlevel:
                sequence += [{tail: lastlevel[tail]}]
                aminosequence += tail

                feasiblePoslist = find_amino_sequence_in_protein(aminosequence, aminos)

                for pos in feasiblePoslist:
                    feasiblePositions[pos] = sequence.copy()
                    startbasispos = pos[0] * 3 + startoffset
                    endbasispos = (pos[1] + 1) * 3 - endoffset - 1
                    feasiblePositionsBasis[startbasispos, endbasispos] = sequence.copy()

                del sequence[-1]
                aminosequence = aminosequence[:-1]

    return feasiblePositionsBasis



def protein_from_solution(x:gb.tupledict,
                          protein:str)->str:

    """
    Build the protein from a solver's solution (x variables only)
    :param: solution x
    :return: protein
    """
    for i in x:
        if x[i] > 0.5:
            position = i[0]
            amino = i[1]
            codon = i[2]
            protein[position] = codon

    return ''.join(protein)


def build_y_variables_index(desired:list,
                            aminos:str)->(list,dict,dict):
    """
    Build the index set of y variables
    :param: list of desired motifs, amino acids sequence
    :return list of y indexes, feasible positions for the basis, index of desired motifs
    """

    indexlist = []
    feasible_positions_desired = {}
    feasible_positions_desired_basis = {}
    index_desired = {}

    count = 1

    for motif in desired:
        index_desired[motif] = count
        count += 1

        aux = find_motif_position(motif, aminos)
        if aux:

            #feasible_positions_desired[motif] = aux
            feasible_positions_desired_basis[motif] = aux

            for pos in feasible_positions_desired_basis[motif]:
                indexlist.append((index_desired[motif], str(pos).replace(" ", "")))

    return indexlist, feasible_positions_desired_basis, index_desired


def build_z_variables_index(forbidden:list,
                            aminos:str)->(list, dict, dict):
    """
    Build the index set of z variables
    :param: list of forbidden motifs, amino acids sequence
    :return list of z indexes, feasible positions for the basis, index of forbidden motifs
    """

    indexlist = []
    feasible_positions_forbidden = {}
    feasible_positions_forbidden_basis = {}
    index_forbidden = {}

    count = 1

    for motif in forbidden:
        index_forbidden[motif] = count
        count += 1

        aux = find_motif_position(motif, aminos)
        if aux:

            #feasible_positions_forbidden[motif] = aux
            feasible_positions_forbidden_basis[motif] = aux

            for pos in feasible_positions_forbidden_basis[motif]:
                indexlist.append((index_forbidden[motif], str(pos).replace(" ", "")))
        # else:
        # print('No feasible positions')
        # print ('-'*50)

    return indexlist, feasible_positions_forbidden_basis, index_forbidden


def count_forbidden(protein:str,
                    forbidden:list)->int:
    """
    Count the number of forbidden motifs in protein
    """
    totalforb = 0

    for motif in forbidden:
        found = protein.find(motif)
        count = 0
        while found != -1:
            count += 1
            start = found + 1
            found = protein.find(motif, start)
        totalforb += count

    return totalforb


def count_desired(protein, desired):
    """
    Count the number of desired motifs in
    """

    totaldes = 0

    positions = set()

    for motif in desired:
        count = 0
        found = protein.find(motif)
        while found != -1:
            positions.add(found)
            count += 1
            start = found + 1
            found = protein.find(motif, start)
        totaldes += count

    return totaldes


def gen_protein_hierarchical_objectives(target_protein, forbidden, desired):
    """
    Generate a protein through the MIP Model
    by optimizing hierarchically:
        1. The number of forbidden motifs (min)
        2. The number of desired motifs (max)
        3. The Codom Adaptation Index (CAI)
    """

    # Setup the GUROBI Model

    model = gb.Model()
    model.Params.OutputFlag = 0
    model.Params.Threads = 4

    indexlist = []
    costs = {}

    protein_length = len(target_protein)

    aminos = ''

    # Build index and costs of x variables

    for h in range(protein_length // 3):
        amino = TableDNA[target_protein[h * 3: h * 3 + 3]]
        aminos += amino
        for i in enumerate(InvTableDNA[amino]):
            indexlist.append((h, TableDNA[target_protein[h * 3: h * 3 + 3]], i[1]))
            costs[indexlist[-1]] = Costs[i[1]]

    x = model.addVars(indexlist, vtype='B', name='x')

    # CAI objective has index 2 and priority 0 (lowest priority)

    model.setObjectiveN(x.prod(costs), 2, 0, name='CAI')
    model.ModelSense = -1

    # Assignment constraints

    model.addConstrs((x.sum(i, '*', '*') == 1 for i in range(protein_length // 3)), name='Pos')

    # Build index of y variables

    indexlist, feasible_positions_desired, index_desired = \
        build_y_variables_index(desired, aminos)

    # y variables

    if indexlist:

        y = model.addVars(indexlist, vtype='B', name='y')

        for motif in feasible_positions_desired:
            for pos in feasible_positions_desired[motif]:
                aux = [codon[h] for codon in feasible_positions_desired[motif][pos] for h in codon]
                for idx, codlist in enumerate(aux):
                    model.addConstr(y[index_desired[motif], str(pos).replace(" ", "")] <=
                                    gb.quicksum([x[pos[0] // 3 + idx, TableDNA[cod], cod] for cod in codlist]),
                                    name='DesMotif[' + str(index_desired[motif]) + ']' + str(pos).replace(" ",
                                                                                                          "") + 'Pos[' + str(
                                        pos[0] // 3 + idx) + ']')

        # Desired motifs objective has index 1 and priority 1 (medium priority)

        model.setObjectiveN(y.prod({i: 1.0 for i in indexlist}), 1, 1, name='Desired')

    #
    # Build index of z variables
    #

    indexlist, feasible_positions_forbidden, index_forbidden = \
        build_z_variables_index(forbidden, aminos)

    # z variables

    if indexlist:
        z = model.addVars(indexlist, vtype='B', name='z')

        for motif in feasible_positions_forbidden:

            for pos in feasible_positions_forbidden[motif]:
                lhs = gb.LinExpr()
                rhs = 0.0
                aux = [codon[h] for codon in feasible_positions_forbidden[motif][pos] for h in codon]
                for idx, codlist in enumerate(aux):
                    lhs += gb.quicksum([x[pos[0] // 3 + idx, TableDNA[cod], cod] for cod in codlist])
                    rhs += 1.0
                model.addConstr(lhs - z[index_forbidden[motif], str(pos).replace(" ", "")]
                                <= rhs - 1,
                                name='ForbMotif[' + str(index_forbidden[motif]) + ']' + str(pos).replace(" ", ""))

        # Desired motifs objective has index 0 and priority 2 (highest priority)

        model.setObjectiveN(z.prod({i: -1.0 for i in indexlist}), 0, 2, name='Forbidden')

    # Attach data to Gurobi model

    model._vars = x
    model._tentativeprotein = ['***' for i in range(protein_length // 3)]

    model._currentprotein = ''
    model._forbids = forbidden

    model.optimize()

    assert model.Status == gb.GRB.Status.OPTIMAL

    xsol = []

    for s in range(model.SolCount):
        model.params.SolutionNumber = s
        xsol.append(model.getAttr('Xn', x))

    #  Get the best solution for CAI evaluation

    model.params.SolutionNumber = 0
    model.params.ObjNumber = 2
    CAI_exp = model.ObjNVal

    finalprotein = protein_from_solution(xsol[0], model._tentativeprotein)
    CAI = math.pow(math.e, CAI_exp / (len(finalprotein) // 3))
    n_forbidden = count_forbidden(finalprotein, forbidden)
    n_desired = count_desired(finalprotein, desired)

    return finalprotein, model.Runtime, model.NodeCount, model.Status, n_forbidden, n_desired, CAI

def main():

    dataset = read_fasta('gencode_filtered.fasta')
    forbidden, desired = read_motifs('forbidden.cpg', 'desirable.cpg')


    count = 0
    print ('Num. prot.  | Length  | Forb. bef. | Des. bef. | CAI bef. | Forb. af.| Des. af.| CAI af. | Time')
    for index, row in dataset.iterrows():
        count += 1

        Protein = row['Protein']

        CAI_before = math.pow(math.e,
                              sum(Costs[Protein[i * 3:i * 3 + 3]] for i in range(len(Protein) // 3)) / (
                                          len(Protein) // 3))

        n_forbidden_start = count_forbidden(Protein, forbidden)
        n_desired_start = count_desired(Protein, desired)

        print('%11d' % count, '|%8d' % len(Protein), '|%11d' % n_forbidden_start, '|%10d' % n_desired_start, ('|   %1.4f' % CAI_before), end=' ')

        dataset.loc[index, 'Forb. before'] = n_forbidden_start
        dataset.loc[index, 'Des. before'] = n_desired_start
        dataset.loc[index, 'CAI before'] = CAI_before

        finalprotein, time, nodes, status, n_forbidden, n_desired, CAI = \
            gen_protein_hierarchical_objectives(Protein, forbidden,
                                                desired)

        dataset.loc[index, 'Forb. opt'] = count_forbidden(finalprotein, forbidden)
        dataset.loc[index, 'Des. opt'] = count_desired(finalprotein, desired)
        dataset.loc[index, 'CAI opt'] = CAI
        dataset.loc[index, 'Time'] = time

        print('|%9d' % dataset.loc[index, 'Forb. opt'], '|%8d' % dataset.loc[index, 'Des. opt'],
              ('| %1.4f' % dataset.loc[index, 'CAI opt']), ('| %4.4f' % time))

    dataset.to_csv('optimized.csv')


if __name__ == '__main__':
    main()