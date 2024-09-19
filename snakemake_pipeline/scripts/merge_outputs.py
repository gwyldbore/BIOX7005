import pandas as pd
import pickle

def load_interproscan_df(path):

    column_names = [
        "info",
        "label",
        "sub_label",
        "description",
        "start",
        "stop",
        "score",
        "status",
        "date",
        "extended_description"
    ]
    selected_columns = [0, 3, 4, 5, 6, 7, 8, 9, 10, 12]

    interpro_df = pd.read_csv(path, sep='\t', usecols=selected_columns, names=column_names)
    return interpro_df


def main():
    # Load the dataframes
    blast_df = pd.read_csv(snakemake.input.blast_df)
    interproscan_df = load_interproscan_df(snakemake.input.interproscan_df)

    with open(snakemake.input.embedding_df, "rb") as input_file:
        embedding_df = pickle.load(input_file)

    merged_df = embedding_df.merge(interproscan_df, on='info', how='left')

    pivot_df = merged_df.pivot_table(index='info', columns='label', values='extended_description',
                                     aggfunc=lambda x: '; '.join(x))
    pivot_df = pivot_df.reset_index()

    # Merge back with embedding dataframe
    final_df = embedding_df.merge(pivot_df, on=['info'], how='left')

    # Merge with blast results
    final_df = final_df.merge(blast_df[['info', 'blast_results']], on='info', how='left')

    final_df['has_subfamily_4'] = final_df['PRINTS'].str.contains('subfamily 4', na=False)
    final_df['has_subfamily_1'] = final_df['PRINTS'].str.contains('subfamily 1', na=False)

    # Save the merged dataframe
    final_df.to_csv(snakemake.output.merged_df, index=False)


if __name__ == "__main__":
    main()
