import pandas as pd

def main():
    # load in the merged df
    merged_df = pd.read_csv(snakemake.input.merged_df)

    results_df = merged_df[['info', 'sequence', 'num_mutation', 'blast_prediction', 'interproscan_prediction', 'embedding_prediction', 'mutated_positions' ]].copy()

    results_df.to_csv(snakemake.output.results_df)


if __name__ == "__main__":
    main()

