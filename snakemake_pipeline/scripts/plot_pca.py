import pickle
import pandas as pd
import seq_utils
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


def plot_pca(all_embeddings_df, nodes_to_label, outpath, col_name='protbert_cls_embedding'):

    embeddings = np.vstack(all_embeddings_df[col_name].values)


    # Apply PCA
    num_components = 2
    pca = PCA(n_components=num_components)
    pca_result = pca.fit_transform(embeddings)

    # Add PCA results to the DataFrame
    all_embeddings_df['pca1'] = pca_result[:, 0]
    all_embeddings_df['pca2'] = pca_result[:, 1]

    # Get unique clades
    clades_with_color = all_embeddings_df['Clade'].dropna().unique()
    num_clades = len(clades_with_color)

    # Define color map for the clades
    colors = plt.cm.get_cmap('Set1', num_clades).colors

    plt.figure(figsize=(20, 14))

    # Plot all points in gray first to show entries with no clade
    no_clade_df = all_embeddings_df[all_embeddings_df['Clade'].isna()]
    # probably a way to set colour=blue to be a gradient??
    plt.scatter(no_clade_df['pca1'], no_clade_df['pca2'], color='blue', alpha=0.5, label='No Clade')

    # Plot points with clades in different colors
    for clade, color in zip(clades_with_color, colors):
        subset = all_embeddings_df[all_embeddings_df['Clade'] == clade]
        plt.scatter(subset['pca1'], subset['pca2'], label=clade, color=color)



    # sequence_no_for_colour = all_embeddings_df['num_mutation'].dropna().unique()
    
    # for seq_no in sequence_no_for_colour:
    #     mutated_subset = all_embeddings_df[all_embeddings_df['num_mutation'] == seq_no]
    #     plt.scatter(mutated_subset['pca1'], mutated_subset['pca2'], c=mutated_subset['num_mutation'], cmap='Blues')

    mutation_df = all_embeddings_df.dropna(subset=['num_mutation'])
    plt.scatter(mutation_df['pca1'], mutation_df['pca2'], 
                c=[int(x) for x in mutation_df['num_mutation']], cmap='cool')
    
    cbar = plt.colorbar()

    # Set plot titles and labels
    plt.title("PCA by Clade")
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.legend()
    plt.savefig(outpath)

def plot_pca_colour_by_predicted(all_embeddings_df, nodes_to_label, outpath, col_name='protbert_cls_embedding'):

    embeddings = np.vstack(all_embeddings_df[col_name].values)


    # Apply PCA
    num_components = 2
    pca = PCA(n_components=num_components)
    pca_result = pca.fit_transform(embeddings)

    # Add PCA results to the DataFrame
    all_embeddings_df['pca1'] = pca_result[:, 0]
    all_embeddings_df['pca2'] = pca_result[:, 1]

    # Get unique clades
    clades_with_color = all_embeddings_df['Clade'].dropna().unique()
    num_clades = len(clades_with_color)

    # Define color map for the clades
    colors = plt.cm.get_cmap('Set1', num_clades).colors

    plt.figure(figsize=(20, 14))

    # Plot all points in gray first to show entries with no clade
    no_clade_df = all_embeddings_df[all_embeddings_df['Clade'].isna()]
    # probably a way to set colour=blue to be a gradient??
    plt.scatter(no_clade_df['pca1'], no_clade_df['pca2'], color='blue', alpha=0.5, label='No Clade')

    # Plot points with clades in different colors
    for clade, color in zip(clades_with_color, colors):
        subset = all_embeddings_df[all_embeddings_df['Clade'] == clade]
        plt.scatter(subset['pca1'], subset['pca2'], label=clade, color=color)



    # Overlay predictions with new colors
    prediction_df = all_embeddings_df.dropna(subset=['overall_prediction'])
    unique_predictions = prediction_df['overall_prediction'].unique()

    # Define a new colormap for predictions
    prediction_colors = plt.cm.cool(np.linspace(0, 1, len(unique_predictions)))

    for prediction, color in zip(unique_predictions, prediction_colors):
        pred_subset = prediction_df[prediction_df['overall_prediction'] == prediction]
        plt.scatter(pred_subset['pca1'], pred_subset['pca2'], 
                    color=color, marker='x', s=100, alpha=0.8, 
                    label=f'Prediction: {prediction}')



    # Set plot titles and labels
    plt.title("PCA by Clade")
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.legend()
    plt.savefig(outpath)


def main():
    # Load the df with the mutated sequences
    with open(snakemake.input.embedding_df, "rb") as input_file:
        embedding_df = pickle.load(input_file)

    # load the df with the prediction info
    with open(snakemake.input.predictions_df, "rb") as prediction_input:
        prediction_df = pd.read_csv(prediction_input)
    embedding_predictions = pd.merge(embedding_df, prediction_df[['info', 'overall_prediction']], 
                                     on='info', how='left')

    nodes_to_label = embedding_df['info'].values
    # print('Nodes to label:', nodes_to_label)

    # Load previously calculated ancestor embeddings
    with open("./data/ancestor_embedding_df.csv", "rb") as input_file:
        ancestor_embedding_df = pickle.load(input_file)
    # with open(snakemake.input.ancestor_embeddings, "rb") as input_file:
    #     ancestor_embedding_df = pickle.load(input_file)

    ancestor_embedding_df['Clade'] = ancestor_embedding_df['info'].apply(seq_utils.tag_node)

    # Filter for only NR1 or NR4 clades
    specific_ancestor_embedding_df = ancestor_embedding_df[ancestor_embedding_df['Clade'].isin(['NR1', 'NR4'])]

    # Concatenate the embeddings and ancestor embeddings
    all_embeddings_df = pd.concat([embedding_df, specific_ancestor_embedding_df])

    # Plot PCA
    plot_pca(all_embeddings_df, nodes_to_label, snakemake.output.plot_mutation)


    all_embeddings_prediction_df = pd.concat([embedding_predictions, specific_ancestor_embedding_df])
    plot_pca_colour_by_predicted(all_embeddings_prediction_df, nodes_to_label, snakemake.output.plot_prediction)

if __name__ == "__main__":
    main()
