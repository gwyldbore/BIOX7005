from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
import random


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

    seq1_out, seq2_out = ''

    for i, character in enumerate(seq1.seq):
        if character == '-' and seq2.seq[i] == '-':
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
    
    return seqrec1, seqrec2




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
                   'N':'polar', 'Q':'polar'}

    pos_mutation = []

    for i, seq1_char in enumerate(seq1.seq):
        seq2_char = seq2.seq[i]

        char1_type = AMINO_TYPES[seq1_char]
        char2_type = AMINO_TYPES[seq2_char]

        if char1_type != char2_type:
            pos_mutation.append((i, seq2_char))

    return pos_mutation




def generate_mutations(inputfile, outputfile, method_type, positions=None, seed=42):
    """
    Assumes inputfile is fasta of two aligned sequences, 
    first is origin and second is target.
    """


    random.seed(random.random())

    mutated_seqs = []

    records = list(SeqIO.parse(inputfile, 'fasta'))
    origin = records[0]
    target = records[1]

    if len(origin.seq) != len(target.seq):
        raise ValueError("These sequences are not aligned")
    

    # remove the gaps that are common (i.e. gaps in both sequences)
    origin, target = remove_common_gaps(origin, target)



    # get the possible mutations for the specified method
    if method_type == 'random':
        possible_mutations = calculate_differences(origin, target)

    elif method_type == 'specified':
        if positions is None:
            raise ValueError("positions were not specified")
        possible_mutations = get_specified_mutations(origin, target, positions)

    elif method_type == 'nonconservative':
        possible_mutations = get_nonconservative_mutations(origin, target)

    # # calculate allowed characters to mutate to
    # if positions is None:
    #     possible_mutations = calculate_differences(origin, target)
    # else:
    #     possible_mutations = get_specified_mutations(origin, target, positions)



    # shuffle the mutations
    random.shuffle(possible_mutations)

    
    # insert first (unmutated) sequence
    first = SeqRecord(
        Seq(str(origin.seq)),
        id=f"{origin.id}_{target.id}_0",
        description=''
    )

    # remove gaps
    # record = remove_gaps(first)
    record = first

    # store sequence as a mutable for alteration
    mutableseq = MutableSeq(str(origin.seq))

    mutated_seqs.append(record)

    i = 1
    for pos, mutation in possible_mutations:
        # mutate the working sequence
        mutableseq[pos] = mutation

        # create a new seqrecord with the mutated sequence
        record = SeqRecord(
            Seq(str(mutableseq)),
            id=f"{origin.id}_{target.id}_{i}",
            description=''
        )

        # record = remove_gaps(record)

        mutated_seqs.append(record) 
        i += 1

    # write all sequences to output file
    SeqIO.write(mutated_seqs, outputfile, 'fasta')




# inputfile = snakemake.input
# outputfile = snakemake.output


# generate_mutations('../data/NR1_NR4_ancestors.fasta', '../data/testoutput.fasta')
# generate_mutations('NR1_NR4_ancestors.fasta', 'testoutput.fasta', [0,1,2,3,4])
generate_mutations(snakemake.input.fasta, snakemake.output.generated_sequences, 'random')

# generate_mutations(snakemake.input.fasta, snakemake.output.fasta, snakemake.wildcards.method_name)

# this can take snakemake.wildcards.method_name as an extra input (make this the method type as a string)
# put this into the generate_mutations signature and do an if else statement for how to get the list of positions
