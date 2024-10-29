from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
import random
import pandas as pd


def calculate_differences(seq1, seq2):
    """
    Calculate the different positions between the two provided sequences, and 
    identify characters to mutate to at those positions.

    Params: 
        seq1: SeqRecord object containing the origin sequence
        seq2: SeqRecord object containing the target sequence

    Returns:
        List[tuple(int, char)] lsit of tuples (position index, target character)
    """

    pos_mutation = []

    for i, character in enumerate(seq1.seq):
        # print(character, seq2.seq[i])
        if character != seq2.seq[i]:

            pos_mutation.append((i, seq2.seq[i]))

    return pos_mutation


def get_specified_mutations(seq1, seq2, positions):
    """
    Retrieve the characters to mutate to at specific positions, if they differ 
    from the origin sequence.

    Params: 
        seq1: SeqRecord object containing the origin sequence
        seq2: SeqRecord object containing the target sequence

    Returns:
        List[tuple(int, char)] lsit of tuples (position index, target character)
    """
    pos_mutation = []

    for position in positions:
        if seq1.seq[position] != seq2.seq[position]:
            pos_mutation.append((position, seq2.seq[position]))

    return pos_mutation


def remove_gaps(sequence):
    """
    Removes gaps from a sequence and returns the ungapped sequence
    """

    ungapped = ''
    for character in sequence.seq:
        if character == '-':
            continue
        else:
            ungapped += character

    record = SeqRecord(
        Seq(ungapped),
        id=f"{sequence.id}",
        description=''
    )

    return record

def remove_common_gaps(seq1, seq2):
    """Removes gaps which are common to both sequences, preserving alignment."""

    seq1_out, seq2_out = '', ''
    removed_at = []


    for i, character in enumerate(seq1.seq):
        if character == '-' and seq2.seq[i] == '-':
            removed_at.append(i)
            continue
        else:
            seq1_out += character
            seq2_out += seq2.seq[i]

    seqrec1 = SeqRecord(
        Seq(seq1_out),
        id=f"{seq1.id}",
        description='')
    
    seqrec2 = SeqRecord(
        Seq(seq2_out),
        id=f"{seq2.id}",
        description='')
    
    return seqrec1, seqrec2, removed_at



def get_nonconservative_mutations(seq1, seq2):
    """
    Get mutations along the whole sequence only if the mutation causes
    the amino acide to change types, e.g. basic to acidic etc. 
    """
    AMINO_TYPES = {'A':'aliphatic', 'G':'aliphatic', 'I':'aliphatic', 
                   'L':'aliphatic', 'P':'aliphatic', 'V':'aliphatic', 
                   'F':'aromatic', 'W':'aromatic', 'Y':'aromatic', 
                   'R':'basic', 'H':'basic', 'K':'basic', 
                   'D':'acidic', 'E':'acidic', 'S':'', 
                   'T':'polar', 'C':'polar', 'M':'polar', 
                   'N':'polar', 'Q':'polar', '-':'gap'}

    pos_mutation = []

    for i, seq1_char in enumerate(seq1.seq):
        seq2_char = seq2.seq[i]

        char1_type = AMINO_TYPES[seq1_char]
        char2_type = AMINO_TYPES[seq2_char]

        if char1_type != char2_type:
            pos_mutation.append((i, seq2_char))

    return pos_mutation



def assign_priority(row): 
    """Assign the mutation priorities to the row of the dataframe based on compared conservation."""  

    # Extract values for comparison
    vc1 = set(row['very_conserved_positions_origin'])
    vc4 = set(row['very_conserved_positions_target'])
    c1 = set(row['conserved_positions_origin'])
    c4 = set(row['conserved_positions_target'])

    if (vc1 and vc4): # if both have a highly conserved site
        # the two are the same - case 1
        if vc1 == vc4:
            return 10
            # don't care about anything else, this is lowest priority

        # the two are different 
        else:
            # one overlaps with the other's conserved - mid priority - case 2
            if (vc1.intersection(c4) or vc4.intersection(c1)):
                return 2

            # no overlap - highest priority - case 3
            else:
                return 1

    elif (not vc1 and not vc4): # neither have highly conserved site

        # both have conserved
        if (c1 and c4):
            # conserved are the same - low priority - case 4
            if c1 == c4:
                return 8

            # conserved intersect - mid priority - case 5
            elif c1.intersection(c4):
                return 6

            # conserved are unique - higher priority - case 6
            else:
                return 4

        # only one has conserved - case 7
        elif c1 or c4:
            return 7
            # mid-low priority

        # neither have conserved - case 8 (neither have anything)
        else:
            return 10
            # this position literally does not matter because there's no conservation


    else: # one or the other has a highly conserved site (vc1 or vc4)

        # both have conserved
        if c1 and c4:

            # conserved is the same - case 9
            if c1 == c4:
                return 6 # means same position exists in c and vc but one has it higher
                # so this is reasonably low priority

            # conserved intersects
            
                # with very conserved - case 10
            elif c1.intersection(vc4) or c4.intersection(vc1):
                # if conserved is same as other very conserved - conservation has mostly remained so low priority
                if (c1 == vc4) or (c4 == vc1):
                    return 6
                # else conserved just intesects with very conserved
                else:
                    return 5

                # with conserved - case 11
            elif c1.intersection(c4):
                return 4

            # conserved is unique - case 12
            else:
                return 3


        # neither have conserved
            # literally impossible 

        # only one has conserved - case 13
            # has to be same one as highly conserved so not helpful/important
            # this is gonna be more important going one direction vs the other
            # leaving it higher for now because going from vc to nothing is a big change
        else:
            return 3
        
def parse_filename(inputfile):
        """helper to get the dataset and sequences out of filename"""
        removedpath = inputfile.split('/')[-1]
        removedextension = removedpath.split('.')[0]
        name_parts = removedextension.split('_')

        dataset = name_parts[0]
        originseq = name_parts[2]
        targetseq = name_parts[3]
        return dataset, originseq, targetseq



def get_probabilistic_mutations(inputfile, removed_at):
    dataset, originseq, targetseq = parse_filename(inputfile)

    # because Index is 1 indexed and removed at was 0 indexed, add 1 to all of them
    removed_at = [x+1 for x in removed_at]

    datapath = '/'.join(inputfile.split('/')[:-1])
    originseq_file = f'{datapath}/{dataset}_{originseq}_marginal.tsv'
    targetseq_file = f'{datapath}/{dataset}_{targetseq}_marginal.tsv'

    df_origin = pd.read_csv(originseq_file, sep='\t')
    df_target = pd.read_csv(targetseq_file, sep='\t')

    """
    nan_rows_origin = df_origin[df_origin.drop(columns=['Index']).isna().all(axis=1)]
    nan_rows_target = df_target[df_target.drop(columns=['Index']).isna().all(axis=1)]

    common_indices_origin = nan_rows_origin['Index'].isin(nan_rows_target['Index'])
    common_indices_target = nan_rows_target['Index'].isin(nan_rows_origin['Index'])

    # this removes the common gaps so that everything lines up with my other sequences
    # Filter out the rows from origin that match the common 'Index' values
    df_origin_cleaned = df_origin[~(df_origin['Index'].isin(nan_rows_origin[common_indices_origin]['Index']))]
    # filter target similarly
    df_target_cleaned = df_target[~(df_target['Index'].isin(nan_rows_target[common_indices_target]['Index']))]

    """

    # a workaround for the marginal distribution having different shit to the actual seq
    df_origin_cleaned = df_origin[~df_origin['Index'].isin(removed_at)]
    df_target_cleaned = df_target[~df_target['Index'].isin(removed_at)]

    # now fix index in both dfs to be 0 indexed and sequential based on sequence
    i = 0

    # # Create a single string by concatenating the highest-value column names across all rows
    # result_string = ''.join(
    #     df_origin.drop(columns=['Index']).apply(
    #     lambda row: '-' if row.isna().all() else row.idxmax(), axis=1
    #     )
    # )
    # Print the final string
    # print(result_string)
    # print('okay so its seeing the input correctly')

    for index, row in df_origin_cleaned.iterrows():
        df_origin_cleaned.at[index, 'Index'] = i
        df_target_cleaned.at[index, 'Index'] = i
        i += 1

    df_origin = df_origin_cleaned
    df_target = df_target_cleaned

    # grab all the conserved/very conserved position aas
    df_origin['conserved_positions'] = df_origin.apply(lambda row: [col for col in df_origin.columns if col != 'Index' 
                                                    and row[col] >= 0.2], axis=1)
    df_target['conserved_positions'] = df_target.apply(lambda row: [col for col in df_target.columns if col != 'Index' 
                                                        and row[col] >= 0.2], axis=1)

    df_origin['very_conserved_positions'] = df_origin.apply(lambda row: [col for col in df_origin.columns if col != 'Index' 
                                                            and col != 'conserved_positions' 
                                                            and row[col] >= 0.85], axis=1)
    df_target['very_conserved_positions'] = df_target.apply(lambda row: [col for col in df_target.columns if col != 'Index' 
                                                            and col != 'conserved_positions'
                                                            and row[col] >= 0.85], axis=1)

    # merge the dataframes with just the position/conservation values
    df_combined = pd.merge(
    df_origin[['Index', 'conserved_positions', 'very_conserved_positions']],
    df_target[['Index', 'conserved_positions', 'very_conserved_positions']],
    on='Index',
    suffixes=('_origin', '_target')
    )

    # now assign the priorities to each position
    df_combined['Priority'] = df_combined.apply(assign_priority, axis=1)

    # invert the priority to use it as a weight
    df_combined['Weight'] = 1 / df_combined['Priority']
    # then normalise the weights
    df_combined['Weight'] = df_combined['Weight'] / df_combined['Weight'].sum()

    # select order of mutation 'randomly' but using the probability weights
    # print('length of df', len(df_combined['Index']))
    mutation_order = df_combined.sample(n=len(df_combined), weights='Weight', replace=False)['Index'].tolist()

    # print(df_combined.sort_values(by='Priority', ascending=True))
    # print(df_combined['Index'])
    # print(mutation_order)

    return mutation_order



def generate_mutations(inputfile, outputfile, mutation_position_output, method_type, positions=None, seed=42):
    """
    Assumes inputfile is fasta of two aligned sequences, 
    first is origin and second is target.

    outputfile is file to write mutated sequences to

    mutation_position_output is the file containing each sequence's mutated positions

    method_type is a string representing which mutation method to perform
    """


    random.seed(random.random())

    mutated_seqs = []

    records = list(SeqIO.parse(inputfile, 'fasta'))
    origin = records[0]
    target = records[1]

    if len(origin.seq) != len(target.seq):
        raise ValueError("These sequences are not aligned")


    # get the possible mutations for the specified method
    if method_type == 'random':
        # remove the gaps that are common (i.e. gaps in both sequences)
        origin, target, _ = remove_common_gaps(origin, target)

        possible_mutations = calculate_differences(origin, target)

        # shuffle the mutations
        random.shuffle(possible_mutations)


    elif method_type == 'specified':
        if positions is None:
            raise ValueError("positions were not specified")

        # remove the gaps that are common (i.e. gaps in both sequences)
        origin, target, _ = remove_common_gaps(origin, target)

        possible_mutations = get_specified_mutations(origin, target, positions)

        # shuffle the mutations
        random.shuffle(possible_mutations)


    elif method_type == 'nonconservative':
        # remove the gaps that are common (i.e. gaps in both sequences)
        origin, target, _ = remove_common_gaps(origin, target)

        possible_mutations = get_nonconservative_mutations(origin, target)

        # shuffle the mutations
        random.shuffle(possible_mutations)


    elif method_type == 'marginal_weights':
        origin, target, removed_at = remove_common_gaps(origin, target)
        # print(f'length origin {len(origin)}, length target {len(target)}')

        mutation_positions = get_probabilistic_mutations(inputfile, removed_at)
        possible_mutations = get_specified_mutations(origin, target, mutation_positions)


    
    # insert first (unmutated) sequence
    first = SeqRecord(
        Seq(str(origin.seq)),
        id=f"{origin.id}_{target.id}_0",
        description=''
    )

    # remove gaps
    record = remove_gaps(first)
    # record = first

    # store sequence as a mutable for alteration
    mutableseq = MutableSeq(str(origin.seq))

    mutated_seqs.append(record)

    i = 1

    cumulative_positions = []
    mutationfile = open(mutation_position_output, 'w')

    # for each mutation possible
    for pos, mutation in possible_mutations:

        # insert a blank line to represent no mutations for the first sequence
        # and to help format later sequences
        mutationfile.write('\n')

        # mutate the working sequence
        mutableseq[pos] = mutation

        # add the position to the cumulative list
        cumulative_positions.append(pos)

        # create a new seqrecord with the mutated sequence
        record = SeqRecord(
            Seq(str(mutableseq)),
            id=f"{origin.id}_{target.id}_{i}",
            description=''
        )

        record = remove_gaps(record)

        mutated_seqs.append(record) 
        
        # write the cumulative sequences to the tracking file
        # with open(mutation_position_output, 'w') as mutationfile:
        toprint = ','.join(str(pos) for pos in cumulative_positions)
        mutationfile.write(toprint)
            
        i += 1
        

    # write all sequences to output file
    SeqIO.write(mutated_seqs, outputfile, 'fasta')




# inputfile = snakemake.input
# outputfile = snakemake.output


# generate_mutations('../data/NR1_NR4_ancestors.fasta', '../data/testoutput.fasta')
# generate_mutations('../data/NR1_NR4_ancestors.fasta', '../data/testoutput.fasta', 'testposlist.txt', 
#                    'specified', [0,1,2,3,4])


# generate_mutations('../../data/reportdata/cd70_NR1toNR4_N6_N81.fasta', '../data/testoutput.fasta', 'testposlist.txt', 
#                    'marginal_weights')



# generate_mutations(snakemake.input.fasta, snakemake.output.fasta, snakemake.wildcards.method_name)

# this can take snakemake.wildcards.method_name as an extra input (make this the method type as a string)
# put this into the generate_mutations signature and do an if else statement for how to get the list of positions


generate_mutations(snakemake.input.fasta, 
                   snakemake.output.generated_sequences, 
                   snakemake.output.mutation_positions, 
                   snakemake.wildcards.method_name)