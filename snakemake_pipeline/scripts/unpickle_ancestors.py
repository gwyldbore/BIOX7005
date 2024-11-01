import pickle


with open('../data/ancestor_embedding_combined_df.csv', 'rb') as input_file:
    embeddings_df = pickle.load(input_file)


if isinstance(embeddings_df, pd.DataFrame):
    embeddings_df.to_csv("ancestor_embedding_combined_nonpickle.csv", index=False)
else:
    # Convert to DataFrame if needed
    df = pd.DataFrame(embeddings_df)
    df.to_csv("ancestor_embedding_combined_nonpickle.csv", index=False)