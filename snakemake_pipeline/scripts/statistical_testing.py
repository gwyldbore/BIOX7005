import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import scikit_posthocs as sp






def main():
    input_files = snakemake.input
    method_name = snakemake.params.method_name

    replicates = [pd.read_csv(file) for file in input_files]
    combined_df = pd.concat(replicates, ignore_index=True)


    # make a boxplot of the stats
    # plt.figure(figsize=(10, 6))
    # sns.boxplot(x='method', y='num_mutation', data=combined_df)
    # plt.title(f'Mutation Counts by {method_name}')
    # plt.ylabel("Number of mutations")
    # plt.savefig(snakemake.output.boxplot)
    # plt.close() # close to save memory

    # Get unique methods for plotting
    methods = combined_df['method'].unique()

    # Calculate the number of rows and columns for the grid layout
    num_methods = len(methods)
    cols = 3  # Define the number of columns in the grid
    rows = math.ceil(num_methods / cols)  # Calculate the number of rows needed

    fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows), squeeze=False)

    # Generate boxplots for each method
    for i, method in enumerate(methods):
        row, col = divmod(i, cols)  # Calculate grid position
        sns.boxplot(
            x='method', y='num_mutation', 
            data=combined_df[combined_df['method'] == method], ax=axes[row, col]
        )
        axes[row, col].set_title(f'Mutation Counts for {method}')

    plt.ylabel("Number of mutations")
    plt.savefig(snakemake.output.boxplot)
    plt.close() # close to save memory



    # 3.1 Shapiro-Wilk test for normality
    for method, data in combined_df.groupby('method'):
        stat, p = stats.shapiro(data['num_mutation'])
        print(f'{method} - Shapiro-Wilk test p-value: {p}')

    # 3.2 Levene’s test for equal variances
    stat, p = stats.levene(
        *[group['num_mutation'].values for name, group in combined_df.groupby('method')]
    )
    print(f'Levene’s test p-value: {p}')

    # 3.3 Hypothesis Testing (ANOVA or Kruskal-Wallis)
    if p > 0.05:  # If variances are equal
        stat, p = stats.f_oneway(
            *[group['num_mutation'].values for name, group in combined_df.groupby('method')]
        )
        print(f'ANOVA test p-value: {p}')
        if p < 0.05:
            tukey = pairwise_tukeyhsd(combined_df['num_mutation'], combined_df['method'])
            print(tukey)
    else:
        stat, p = stats.kruskal(
            *[group['num_mutation'].values for name, group in combined_df.groupby('method')]
        )
        print(f'Kruskal-Wallis test p-value: {p}')
        if p < 0.05:
            dunn = sp.posthoc_dunn(combined_df, val_col='num_mutation', group_col='method', p_adjust='bonferroni')
            print(dunn)

    # Step 4: Chi-square Test for Mutated Positions
    positions_df = pd.crosstab(combined_df['Position'], combined_df['method'])
    chi2, p, _, _ = stats.chi2_contingency(positions_df)
    print(f'Chi-square test p-value: {p}')









if __name__ == "__main__":
    main()