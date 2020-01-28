import gurobipy as gb
import pandas as pd
import math
import argparse


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


def read_fasta(fasta_file: str) -> pd.DataFrame:
    """
    Read fasta file using pandas
    :param: filename
    :return: pandas DataFrame with target proteins
    """
    df = pd.read_csv(fasta_file, header=None, comment='>', engine='python',
                     names=['Protein'])

    df['Protein'] = df['Protein'].apply(lambda x: x[:len(x) - len(x) % 3])

    df['Length'] = df['Protein'].str.len()

    print(fasta_file, ('has %d protein(s)' % df.shape[0]))

    return df


def read_motifs(forbidden_file: str,
                desired_file: str) -> (list, list):
    """
    Read forbidden and desired motifs
    :param: forbidden and desired file names
    :return: lists of forbidden and desired motifs
    """

    forbidden = []
    desired = []

    if forbidden_file is not None:
        with open(forbidden_file) as file:
            forbidden = file.readlines()

        forbidden = [x.strip() for x in forbidden]

    if desired_file is not None:
        with open(desired_file) as file:
            desired = file.readlines()

        desired = [x.strip() for x in desired]

    return forbidden, desired


def find_amino_from_codon(codon: str,
                          table: dict) -> dict:

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
        partial_codon = codon.lstrip('*')
        offset = len(codon) - len(partial_codon)
        for i in table:
            if partial_codon == i[offset:]:
                if table[i] in aminos:
                    aminos[table[i]].append(i)
                else:
                    aminos[table[i]] = [i]

    if codon[2] == '*':
        partial_codon = codon.rstrip('*')
        offset = len(codon) - len(partial_codon)

        for i in table:
            if partial_codon == i[:3 - offset]:
                if table[i] in aminos:
                    aminos[table[i]].append(i)
                else:
                    aminos[table[i]] = [i]

    return aminos


def find_amino_sequence_in_protein(sequence: str,
                                   aminos: str) -> list:
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


def find_motif_position(motif: str,
                        aminos: str) -> (dict, dict):
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

    feasible_positions = {}
    feasible_positions_basis = {}

    #
    # Test three possible positions
    #

    for offset in range(3):

        sequence = ''.join('*' * offset) + motif
        start_offset = offset

        if len(sequence) % 3 != 0:
            end_offset = (3 - (len(sequence) % 3))
            sequence += ('*' * (3 - (len(sequence) % 3)))
        else:
            end_offset = 0

        first_level = find_amino_from_codon(sequence[0:3], TableDNA)

        if len(sequence) > 3:
            last_level = find_amino_from_codon(sequence[-3:], TableDNA)
        else:
            raise ValueError('Motif must have length > 3')

        middle_levels = []

        for idx in range(3, len(sequence) - 3, 3):
            codon = sequence[idx: idx + 3]
            middle_levels.append(find_amino_from_codon(codon, TableDNA))

        #
        # Build the sequences starting from possible
        # encodings
        #

        for head in first_level:
            sequence = [{head: first_level[head]}]

            sequence += middle_levels

            amino_sequence = ''.join([[*i.keys()][0] for i in middle_levels])
            amino_sequence = head + amino_sequence

            for tail in last_level:
                sequence += [{tail: last_level[tail]}]
                amino_sequence += tail

                feasible_pos_list = find_amino_sequence_in_protein(amino_sequence, aminos)

                for pos in feasible_pos_list:
                    feasible_positions[pos] = sequence.copy()
                    start_basis_pos = pos[0] * 3 + start_offset
                    end_basis_pos = (pos[1] + 1) * 3 - end_offset - 1
                    feasible_positions_basis[start_basis_pos, end_basis_pos] = sequence.copy()

                del sequence[-1]
                amino_sequence = amino_sequence[:-1]

    return feasible_positions_basis


def protein_from_solution(x: gb.tupledict,
                          protein: list) -> str:

    """
    Build the protein from a solver's solution (x variables only)
    :param: solution x
    :return: protein
    """
    for i in x:
        if x[i] > 0.5:
            position = i[0]
            codon = i[2]
            protein[position] = codon

    return ''.join(protein)


def build_y_variables_index(desired: list,
                            aminos: str) -> (list, dict, dict):
    """
    Build the index set of y variables
    :param: list of desired motifs, amino acids sequence
    :return list of y indexes, feasible positions for the basis, index of desired motifs
    """

    index_list = []
    feasible_positions_desired_basis = {}
    index_desired = {}

    count = 1

    for motif in desired:
        index_desired[motif] = count
        count += 1

        aux = find_motif_position(motif, aminos)
        if aux:

            feasible_positions_desired_basis[motif] = aux

            for pos in feasible_positions_desired_basis[motif]:
                index_list.append((index_desired[motif], str(pos).replace(" ", "")))

    return index_list, feasible_positions_desired_basis, index_desired


def build_z_variables_index(forbidden: list,
                            aminos: str) -> (list, dict, dict):
    """
    Build the index set of z variables
    :param: list of forbidden motifs, amino acids sequence
    :return list of z indexes, feasible positions for the basis, index of forbidden motifs
    """

    index_list = []
    feasible_positions_forbidden_basis = {}
    index_forbidden = {}

    count = 1

    for motif in forbidden:
        index_forbidden[motif] = count
        count += 1

        aux = find_motif_position(motif, aminos)
        if aux:

            feasible_positions_forbidden_basis[motif] = aux

            for pos in feasible_positions_forbidden_basis[motif]:
                index_list.append((index_forbidden[motif], str(pos).replace(" ", "")))

    return index_list, feasible_positions_forbidden_basis, index_forbidden


def count_forbidden(protein: str,
                    forbidden: list) -> int:
    """
    Count the number of forbidden motifs in protein
    """
    total_forb = 0

    for motif in forbidden:
        found = protein.find(motif)
        count = 0
        while found != -1:
            count += 1
            start = found + 1
            found = protein.find(motif, start)
        total_forb += count

    return total_forb


def count_desired(protein, desired):
    """
    Count the number of desired motifs in
    """

    total_des = 0

    positions = set()

    for motif in desired:
        count = 0
        found = protein.find(motif)
        while found != -1:
            positions.add(found)
            count += 1
            start = found + 1
            found = protein.find(motif, start)
        total_des += count

    return total_des


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

    index_list = []
    costs = {}

    protein_length = len(target_protein)

    aminos = ''

    # Build index and costs of x variables

    for h in range(protein_length // 3):
        amino = TableDNA[target_protein[h * 3: h * 3 + 3]]
        aminos += amino
        for i in enumerate(InvTableDNA[amino]):
            index_list.append((h, TableDNA[target_protein[h * 3: h * 3 + 3]], i[1]))
            costs[index_list[-1]] = Costs[i[1]]

    x = model.addVars(index_list, vtype='B', name='x')

    # CAI objective has index 2 and priority 0 (lowest priority)

    model.setObjectiveN(x.prod(costs), 2, 0, name='CAI')
    model.ModelSense = -1

    # Assignment constraints

    model.addConstrs((x.sum(i, '*', '*') == 1 for i in range(protein_length // 3)), name='Pos')

    # Build index of y variables

    index_list, feasible_positions_desired, index_desired = \
        build_y_variables_index(desired, aminos)

    # y variables

    if index_list:

        y = model.addVars(index_list, vtype='B', name='y')

        for motif in feasible_positions_desired:
            for pos in feasible_positions_desired[motif]:
                aux = [codon[h] for codon in feasible_positions_desired[motif][pos] for h in codon]
                for idx, codlist in enumerate(aux):
                    model.addConstr(y[index_desired[motif], str(pos).replace(" ", "")] <=
                                    gb.quicksum([x[pos[0] // 3 + idx, TableDNA[cod], cod] for cod in codlist]),
                                    name='DesMotif[' + str(index_desired[motif]) + ']' +
                                         str(pos).replace(" ", "") + 'Pos[' + str(pos[0] // 3 + idx) + ']')

        # Desired motifs objective has index 1 and priority 1 (medium priority)

        model.setObjectiveN(y.prod({i: 1.0 for i in index_list}), 1, 1, name='Desired')

    #
    # Build index of z variables
    #

    index_list, feasible_positions_forbidden, index_forbidden = \
        build_z_variables_index(forbidden, aminos)

    # z variables

    if index_list:
        z = model.addVars(index_list, vtype='B', name='z')

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

        model.setObjectiveN(z.prod({i: -1.0 for i in index_list}), 0, 2, name='Forbidden')

    # Attach data to Gurobi model

    model._vars = x
    model._tentative_protein = list('***' for i in range(protein_length // 3))

    model._current_protein = ''
    model._forbids = forbidden

    model.optimize()

    assert model.Status == gb.GRB.Status.OPTIMAL

    x_sol = []

    for s in range(model.SolCount):
        model.params.SolutionNumber = s
        x_sol.append(model.getAttr('Xn', x))

    #  Get the best solution for CAI evaluation

    model.params.SolutionNumber = 0
    model.params.ObjNumber = 2
    cai_exp = model.ObjNVal

    final_protein = protein_from_solution(x_sol[0], model._tentative_protein)
    cai = math.pow(math.e, cai_exp / (len(final_protein) // 3))
    n_forbidden = count_forbidden(final_protein, forbidden)
    n_desired = count_desired(final_protein, desired)

    return final_protein, model.Runtime, model.NodeCount, model.Status, n_forbidden, n_desired, cai


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('target_file', help='Target file (FASTA format)')
    parser.add_argument('--forbidden', '-f', help='File of forbidden motifs')
    parser.add_argument('--desired', '-d', help='File of desired motifs')
    parser.add_argument('--output', '-o', help='Output file', default='optimized')

    args = parser.parse_args()

    dataset = read_fasta(args.target_file)

    forbidden, desired = read_motifs(args.forbidden, args.desired)

    count = 0
    print('Num. prot.  | Length  | Forb. bef. | Des. bef. | CAI bef. | Forb. af.| Des. af.| CAI af. | Time')
    for index, row in dataset.iterrows():
        count += 1

        protein = row['Protein']

        cai_before = math.pow(math.e,
                              sum(Costs[protein[i * 3:i * 3 + 3]] for i in range(len(protein) // 3)) / (
                                          len(protein) // 3))

        n_forbidden_start = count_forbidden(protein, forbidden)
        n_desired_start = count_desired(protein, desired)

        print('%11d' % count, '|%8d' % len(protein), '|%11d' % n_forbidden_start,
              '|%10d' % n_desired_start, ('|   %1.4f' % cai_before), end=' ')

        dataset.loc[index, 'Forb. before'] = n_forbidden_start
        dataset.loc[index, 'Des. before'] = n_desired_start
        dataset.loc[index, 'CAI before'] = cai_before

        final_protein, time, nodes, status, n_forbidden, n_desired, cai = \
            gen_protein_hierarchical_objectives(protein, forbidden,
                                                desired)

        dataset.loc[index, 'Forb. opt'] = count_forbidden(final_protein, forbidden)
        dataset.loc[index, 'Des. opt'] = count_desired(final_protein, desired)
        dataset.loc[index, 'CAI opt'] = cai
        dataset.loc[index, 'Time'] = time

        print('|%9d' % dataset.loc[index, 'Forb. opt'], '|%8d' % dataset.loc[index, 'Des. opt'],
              ('| %1.4f' % dataset.loc[index, 'CAI opt']), ('| %4.4f' % time))

    dataset.to_csv(args.output + '.csv')


if __name__ == '__main__':
    main()

