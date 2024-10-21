import pickle
import seq_utils


def load_dataframe(inputfile):

    # with open(snakemake.input.embedding_df, "rb") as input_file:
    #     embedding_df = pickle.load(input_file)

    # rb is read bytes
    # this loads it into a pandas df
    with open(inputfile, "rb") as input_file:
        embedding_df = pickle.load(input_file)

    # apply the family labels to the nodes
    embedding_df['Clade'] = embedding_df['info'].apply(seq_utils.tag_node)

    return embedding_df



df = load_dataframe('../../familyprediction/Nuclear_Receptor_Notebooks/embdding_df.csv')

