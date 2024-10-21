from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
import random

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



def main():
    inputfile = snakemake.input.generated_sequences
    outputfile = snakemake.output.generated_sequences_ungapped

    records = list(SeqIO.parse(inputfile, 'fasta'))
    all_seqs = []

    for sequence in records:
        current = SeqRecord(
        Seq(str(sequence.seq)),
        id=f"{sequence.id}",
        description=''
    )
        
        record = remove_gaps(current)
        all_seqs.append(record)


        # write all sequences to output file
        SeqIO.write(all_seqs, outputfile, 'fasta')




if __name__ == "__main__":
    main()