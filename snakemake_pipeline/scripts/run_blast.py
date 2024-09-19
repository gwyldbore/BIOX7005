import random
import seq_utils

def dummy_blast_results(seq_df, output_df):
    # Randomly assign either 'dummy_NR1' or 'dummy_NR4' to each sequence
    seq_df['blast_results'] = [random.choice(['dummy_NR1', 'dummy_NR4']) for _ in range(len(seq_df))]

    seq_df.to_csv(output_df, index=False)


def main():
    fasta = snakemake.input.generated_sequences
    output_df = snakemake.output.blast_df

    seq_df = seq_utils.get_sequence_df(fasta)

    dummy_blast_results(seq_df, output_df)


if __name__ == "__main__":
    main()
