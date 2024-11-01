import pickle
import pandas as pd
import seq_utils
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from sklearn.decomposition import PCA
from matplotlib.colors import LinearSegmentedColormap, ListedColormap


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
    # colors = plt.cm.get_cmap('Set1', num_clades).colors
    # colors = plt.colormaps['PiYG'].resampled(num_clades)
    # clade_cmap = ListedColormap(['#d3d3d3', '#a9a9a9'])

    clade_cmap = ListedColormap(['#d3d3d3', '#ffca8a'])
    colors = plt.get_cmap(clade_cmap, num_clades).colors

    # plt.figure(figsize=(20, 14))
    fig, ax = plt.subplots(figsize=(20,14))

    # Plot all points in gray first to show entries with no clade
    no_clade_df = all_embeddings_df[all_embeddings_df['Clade'].isna()]
    # probably a way to set colour=blue to be a gradient??
    # plt.scatter(no_clade_df['pca1'], no_clade_df['pca2'], color='blue', alpha=0.5, label='No Clade')
    ax.scatter(no_clade_df['pca1'], no_clade_df['pca2'], color='blue', alpha=0.5, label='No Clade')

    # Plot points with clades in different colors
    for clade, color in zip(clades_with_color, colors):
        subset = all_embeddings_df[all_embeddings_df['Clade'] == clade]
        # plt.scatter(subset['pca1'], subset['pca2'], label=clade, color=color)
        ax.scatter(subset['pca1'], subset['pca2'], label=clade, color=color)


    mutation_df = all_embeddings_df.dropna(subset=['num_mutation'])

    # viridis = plt.get_cmap('viridis')
    # new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    #     f"trunc({'viridis'},{0.2},{0.8})",
    #     viridis(np.linspace(0.2, 0.8, 256))
    # )

    # scatter = plt.scatter(mutation_df['pca1'], mutation_df['pca2'], 
    #             c=[int(x) for x in mutation_df['num_mutation']], cmap='cool')
    scatter = ax.scatter(mutation_df['pca1'], mutation_df['pca2'], 
                c=[int(x) for x in mutation_df['num_mutation']], cmap='cool')
    
    cax = ax.inset_axes([0.05, 0.05, 0.3, 0.05])
    fig.colorbar(scatter, cax=cax, orientation='horizontal')
    
    # cbar = plt.colorbar()
    # cbar_ax = fig.add_subplot(gs[1])
    # fig.colorbar(scatter, cax=cbar_ax, orientation='vertical')
    # cbar_ax.set_ylabel('Colorbar')

    # Set plot titles and labels
    plt.title("PCA by Clade")
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.legend()
    plt.savefig(outpath)

def plot_pca_ancestors_static(mutations_df, ancestors_df, nodes_to_label, outpath, col_name='protbert_cls_embedding'):
    for _, row in ancestors_df:
        if row['Clade'] == 'NR4':
            print(row)
    ancestor_embeddings = np.vstack(ancestors_df[col_name].values)

    pca = PCA(n_components=2)
    pca_result = pca.fit(ancestor_embeddings)


    # Transform both ancestors_df and mutations_df using the fitted PCA
    ancestors_df[['pca1', 'pca2']] = pca.transform(ancestor_embeddings)
    mutations_df[['pca1', 'pca2']] = pca.transform(np.vstack(mutations_df[col_name].values))


    # Set up the plot
    fig, ax = plt.subplots(figsize=(20, 14))

    # Plot mutations_df entries in blue (No Clade)
    # ax.scatter(mutations_df['pca1'], mutations_df['pca2'], color='blue', alpha=0.5)

    # Plot entries from ancestors_df with clades in different colors
    clades_with_color = ancestors_df['Clade'].unique()
    num_clades = len(clades_with_color)

    clade_cmap = ListedColormap(['#d3d3d3', '#ffca8a'])
    colors = plt.get_cmap(clade_cmap, num_clades).colors

    # Plot points with clades in different colors
    for clade, color in zip(clades_with_color, colors):
        subset = ancestors_df[ancestors_df['Clade'] == clade]
        ax.scatter(subset['pca1'], subset['pca2'], label=clade, color=color)


    mutation_df = mutations_df.dropna(subset=['num_mutation'])
    scatter = ax.scatter(mutation_df['pca1'], mutation_df['pca2'], 
                c=[int(x) for x in mutation_df['num_mutation']], cmap='cool')
    
    cax = ax.inset_axes([0.05, 0.05, 0.3, 0.05])
    fig.colorbar(scatter, cax=cax, orientation='horizontal')

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
    # colors = plt.cm.get_cmap('Set1', num_clades).colors
    # clade_cmap = ListedColormap(['royalblue', 'green'])
    clade_cmap = ListedColormap(['#d3d3d3', '#ffca8a'])
    colors = plt.get_cmap(clade_cmap, num_clades).colors

    plt.figure(figsize=(20, 14))
    # fig, ax = plt.subplots(figsize=(20, 14))

    # Plot all points in gray first to show entries with no clade
    no_clade_df = all_embeddings_df[all_embeddings_df['Clade'].isna()]
    # probably a way to set colour=blue to be a gradient??
    plt.scatter(no_clade_df['pca1'], no_clade_df['pca2'], color='blue', alpha=0.5, label='No Clade')

    # Plot points with clades in different colors
    for clade, color in zip(clades_with_color, colors):
        subset = all_embeddings_df[all_embeddings_df['Clade'] == clade]
        plt.scatter(subset['pca1'], subset['pca2'], label=f'Clade: {clade}', color=color)
        # ax.scatter(subset['pca1'], subset['pca2'], label=f'Clade: {clade}', color=color)



    # Overlay predictions with new colors
    prediction_df = all_embeddings_df.dropna(subset=['overall_prediction'])
    unique_predictions = prediction_df['overall_prediction'].unique()
    unique_predictions.sorted()
    # print(unique_predictions)

    # Define a new colormap for predictions
    prediction_cmap = ListedColormap(['mediumorchid', 'red', 'royalblue', 'forestgreen'])
    # if unique_predictions[0] == 'NR1':
    #     prediction_cmap = ListedColormap(['mediumorchid', 'red', 'royalblue', 'forestgreen'])
    # else:
    #     prediction_cmap = ListedColormap(['forestgreen', 'royalblue', 'red', 'mediumorchid'])
        

    prediction_colors = plt.get_cmap(prediction_cmap).colors

    for prediction, color in zip(unique_predictions, prediction_colors):
        pred_subset = prediction_df[prediction_df['overall_prediction'] == prediction]
        plt.scatter(pred_subset['pca1'], pred_subset['pca2'], 
                    color=color, label=f'Prediction: {prediction}')

    # Set plot titles and labels
    plt.title("PCA by Clade")
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.legend()
    plt.savefig(outpath)




def plot_pca_colour_by_predicted_ancestors_static(mutations_df, ancestors_df, nodes_to_label, outpath, col_name='protbert_cls_embedding'):
    ancestor_embeddings = np.vstack(ancestors_df[col_name].values)

    pca = PCA(n_components=2)
    pca_result = pca.fit(ancestor_embeddings)


    # Transform both ancestors_df and mutations_df using the fitted PCA
    ancestors_df[['pca1', 'pca2']] = pca.transform(ancestor_embeddings)
    mutations_df[['pca1', 'pca2']] = pca.transform(np.vstack(mutations_df[col_name].values))


    # Set up the plot
    fig, ax = plt.subplots(figsize=(20, 14))

    # Plot mutations_df entries in blue (No Clade)
    # ax.scatter(mutations_df['pca1'], mutations_df['pca2'], color='blue', alpha=0.5)

    # Plot entries from ancestors_df with clades in different colors
    clades_with_color = ancestors_df['Clade'].unique()
    num_clades = len(clades_with_color)

    clade_cmap = ListedColormap(['#d3d3d3', '#ffca8a'])
    colors = plt.get_cmap(clade_cmap, num_clades).colors

    # Plot points with clades in different colors
    for clade, color in zip(clades_with_color, colors):
        subset = ancestors_df[ancestors_df['Clade'] == clade]
        ax.scatter(subset['pca1'], subset['pca2'], label=clade, color=color)


    # Overlay predictions with new colors
    prediction_df = mutations_df.dropna(subset=['overall_prediction'])
    unique_predictions = prediction_df['overall_prediction'].unique()

    # define new cmap
    prediction_cmap = ListedColormap(['mediumorchid', 'red', 'royalblue', 'forestgreen'])
    if unique_predictions[0] == 'NR4':
        prediction_cmap = ListedColormap(['forestgreen', 'royalblue', 'red', 'mediumorchid'])

    prediction_colors = plt.get_cmap(prediction_cmap, len(unique_predictions)).colors

    for prediction, color in zip(unique_predictions, prediction_colors):
        pred_subset = prediction_df[prediction_df['overall_prediction'] == prediction]
        plt.scatter(pred_subset['pca1'], pred_subset['pca2'], 
                    color=color, label=f'Prediction: {prediction}')

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
    # with open("./data/ancestor_embedding_df.csv", "rb") as input_file:
    #     ancestor_embedding_df = pickle.load(input_file)
    with open(snakemake.input.ancestor_embeddings, "rb") as input_file:
        ancestor_embedding_df = pickle.load(input_file)

    dataset_name = snakemake.wildcards.dataset_name

    # ancestor_embedding_df['Clade'] = ancestor_embedding_df['info'].apply(seq_utils.tag_node, dataset=dataset_name)
    ancestor_embedding_df['Clade'] = ancestor_embedding_df['info'].apply(seq_utils.tag_node, dataset='combined')
    print(ancestor_embedding_df['Clade'])
    # Filter for only NR1 or NR4 clades
    # specific_ancestor_embedding_df = ancestor_embedding_df[ancestor_embedding_df['Clade'].isin(['NR1', 'NR4'])]
    specific_ancestor_embedding_df = ancestor_embedding_df


    for _, row in specific_ancestor_embedding_df.iterrows():
        if row['Clade'] == 'NR4':
            print(row)

    # Concatenate the embeddings and ancestor embeddings
    # all_embeddings_df = pd.concat([embedding_df, specific_ancestor_embedding_df])

    # Plot PCA
    # plot_pca(all_embeddings_df, nodes_to_label, snakemake.output.plot_mutation)
    plot_pca_ancestors_static(embedding_df, specific_ancestor_embedding_df, nodes_to_label, snakemake.output.plot_mutation)


    # all_embeddings_prediction_df = pd.concat([embedding_predictions, specific_ancestor_embedding_df])
    # mutation_prediction_df = pd.concat([embedding_predictions, embedding_df])
    # plot_pca_colour_by_predicted(all_embeddings_prediction_df, nodes_to_label, snakemake.output.plot_prediction)
    plot_pca_colour_by_predicted_ancestors_static(embedding_predictions, specific_ancestor_embedding_df, nodes_to_label, snakemake.output.plot_prediction)






if __name__ == "__main__":
    main()
