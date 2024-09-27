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



def generate_mutations(inputfile, outputfile, positions=None, seed=42):
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


    # calculate allowed characters to mutate to
    if positions is None:
        possible_mutations = calculate_differences(origin, target)
    else:
        possible_mutations = get_specified_mutations(origin, target, positions)



    # shuffle the mutations
    random.shuffle(possible_mutations)

    
    # insert first (unmutated) sequence
    first = SeqRecord(
        Seq(str(origin.seq)),
        id=f"{origin.id}_{target.id}_0",
        description=''
    )

    # remove gaps
    record = remove_gaps(first)

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

        record = remove_gaps(record)

        mutated_seqs.append(record) 
        i += 1

    # write all sequences to output file
    SeqIO.write(mutated_seqs, outputfile, 'fasta')




# inputfile = snakemake.input
# outputfile = snakemake.output


# generate_mutations('../data/NR1_NR4_ancestors.fasta', '../data/testoutput.fasta')
# generate_mutations('NR1_NR4_ancestors.fasta', 'testoutput.fasta', [0,1,2,3,4])
generate_mutations(snakemake.input.fasta, snakemake.output.generated_sequences, [0,1,2])
