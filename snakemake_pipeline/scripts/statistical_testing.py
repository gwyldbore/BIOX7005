import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import scikit_posthocs as sp






def main():
    input_files = snakemake.input
    method_name = snakemake.wildcards.method_name

    replicates = [pd.read_csv(file) for file in input_files]
    combined_df = pd.concat(replicates, ignore_index=True)


    # make a boxplot of the stats
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='method', y='num_mutation', data=combined_df)
    plt.title(f'Mutation Counts by {method_name}')
    plt.savefig(snakemake.output.boxplot)
    plt.close() # close to save memory











if __name__ == "__main__":
    main()