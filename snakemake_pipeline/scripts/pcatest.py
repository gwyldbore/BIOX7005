import pickle
import pandas as pd
import seq_utils
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from sklearn.decomposition import PCA
from matplotlib.colors import LinearSegmentedColormap, ListedColormap





def plot_pca_ancestors_static(mutations_df, ancestors_df, nodes_to_label, outpath, col_name='protbert_cls_embedding'):
    print('inside plot')
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

    marker_dict = {0: 'o',
                1: 'x',
                2: '*',
                3: '^',
                4: 's'}
    

    mutation_df = mutations_df.dropna(subset=['num_mutation'])
    scatter = ax.scatter(mutation_df['pca1'], mutation_df['pca2'], 
                c=[int(x) for x in mutation_df['num_mutation']], cmap='cool',
                marker=marker_dict[mutation_df['filegroup']])
    
    
    # Plot points with different markers based on 'filegroup' column
    for group, marker in marker_dict.items():
        group_subset = mutation_df[mutation_df['filegroup'] == group]
        scatter = ax.scatter(
            group_subset['pca1'], group_subset['pca2'],
            c=[int(x) for x in group_subset['num_mutation']],
            cmap='cool', marker=marker, alpha=0.8
        )
    
    cax = ax.inset_axes([0.05, 0.05, 0.3, 0.05])
    fig.colorbar(scatter, cax=cax, orientation='horizontal')

    # Set plot titles and labels
    plt.title("PCA by Clade")
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.legend()
    plt.savefig(outpath)





def main():
    input_embeddings = [f'../workflows/cd80/marginal_weights/embeddings/NR1toNR4_N7_N186_{rep}.csv' for rep in range(1, 6)]
    input_predictions = [f'../workflows/cd80/marginal_weights/predictions/NR1toNR4_N7_N186_{rep}.csv' for rep in range(1, 6)]

    emb_list = []
    pred_list = []

    def adjust_info(row, index):
        row['info'] += f'_{index}'

    def adjust_filegroup(row, index):
        return index

    i = 0
    print('opening embeddings')
    for file in input_embeddings:
        with open(file, 'rb') as input_file:
            df = pickle.load(input_file)
        # df['filegroup'] = 0
        df.apply(adjust_info, index=i, axis=1)
        df['filegroup'] = df.apply(adjust_filegroup, index=i, axis=1)
        emb_list.append(df)
        
        i += 1
    embedding_dfs = pd.concat(emb_list)


    # i = 0
    # for file in input_predictions:
    #     df = pd.read_csv(input_file)
    #     df.apply(adjust_info, index=i, axis=1)
        # df.apply(adjust_filegroup, index=i, axis=1)
        # i += 1
    #     pred_list.append(df)
    # prediction_dfs = pd.concat(pred_list)

    nodes_to_label = embedding_dfs['info'].values

    print('opening ancestor')
    with open("../data/ancestor_embedding_df.csv", "rb") as input_file:
        ancestor_embedding_df = pickle.load(input_file)

    print('tagging nodes')
    ancestor_embedding_df['Clade'] = ancestor_embedding_df['info'].apply(seq_utils.tag_node, dataset='cd80')
    specific_ancestor_embedding_df = ancestor_embedding_df[ancestor_embedding_df['Clade'].isin(['NR1', 'NR4'])]

    plot_pca_ancestors_static(embedding_dfs, specific_ancestor_embedding_df, nodes_to_label, '../workflows/cd80/groupedpcaplot_marginalsymbol.svg')




if __name__ == "__main__":
    main()
