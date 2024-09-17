from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
import random


def calculate_differences(seq1, seq2):
    positions = []

    for i, character in enumerate(seq1.seq):
        # print(character, seq2.seq[i])
        if character != seq2.seq[i]:

            positions.append((i, seq2.seq[i]))

    return positions

def generate_mutations(inputfile, outputfile, positions=None, seed=42):
    """
    Assumes inputfile is fasta of two aligned sequences, first is origin and second is target.
    """
    random.seed(seed)

    mutated_seqs = []

    records = list(SeqIO.parse(inputfile, 'fasta'))
    origin = records[0]
    target = records[1]

    if len(origin.seq) != len(target.seq):
        raise ValueError("These sequences are not aligned")

    possible_mutations = calculate_differences(origin, target)


    random.shuffle(possible_mutations)

    
    record = SeqRecord(
        Seq(str(origin.seq)),
        id=f"{origin.id}_{target.id}_0",
        description=''
    )
    mutableseq = MutableSeq(str(origin.seq))

    mutated_seqs.append(record)


    i = 1
    for pos, mutation in possible_mutations:
        mutableseq[pos] = mutation

        record = SeqRecord(
            Seq(str(mutableseq)),
            id=f"{origin.id}_{target.id}_{i}",
            description=''
        )

        mutated_seqs.append(record) 

        i += 1




    SeqIO.write(mutated_seqs, outputfile, 'fasta')




# inputfile = snakemake.input
# outputfile = snakemake.output


generate_mutations('NR1_NR4_ancestors.fasta', 'testoutput.fasta')