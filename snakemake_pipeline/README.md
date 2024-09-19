Notes 

Running snakemake

- Need snakemake installed + all the python moduels it calls
- run 'snakemake --cores 1' from the directory with snakefile in it - can specify different number of cores

- Currently interproscan is hardcoded to run on the server (so this part needs to be run on the server)


- Needs 'ancestor_embedding_df.csv' to be placed in data folder (this was the embedding_df.csv, just renamed for clarity)
- Only running Gene3D and PRINTS interproscan atm (can update this in the call to interproscan in snakefile)
- Running the local version of interproscan
- The embeddings are saved weirdly in the merged df
- Currently needs the NR1_ids.txt / NR4_ids.txt and uses these to classify the original ancestor df
- Quite a bit of hardcoded stuff, especially in plot_pca + merge_outputs
- run_blast is just a dummy script for now
- The random seed is just set new each time, but should either be specified or saved somewhere
- Not currently deleting the /temp folder that is created after running interproscan
- Can't open the ancestor_embedding_df.csv on the server??