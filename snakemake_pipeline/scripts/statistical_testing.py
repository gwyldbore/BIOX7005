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
    plt.ylabel("Number of mutations")
    plt.savefig(snakemake.output.boxplot)
    plt.close() # close to save memory



    # 3.1 Shapiro-Wilk test for normality
    for method, data in combined_df.groupby('Method'):
        stat, p = stats.shapiro(data['Mutation_Count'])
        print(f'{method} - Shapiro-Wilk test p-value: {p}')

    # 3.2 Levene’s test for equal variances
    stat, p = stats.levene(
        *[group['Mutation_Count'].values for name, group in combined_df.groupby('Method')]
    )
    print(f'Levene’s test p-value: {p}')

    # 3.3 Hypothesis Testing (ANOVA or Kruskal-Wallis)
    if p > 0.05:  # If variances are equal
        stat, p = stats.f_oneway(
            *[group['Mutation_Count'].values for name, group in combined_df.groupby('Method')]
        )
        print(f'ANOVA test p-value: {p}')
        if p < 0.05:
            tukey = pairwise_tukeyhsd(combined_df['Mutation_Count'], combined_df['Method'])
            print(tukey)
    else:
        stat, p = stats.kruskal(
            *[group['Mutation_Count'].values for name, group in combined_df.groupby('Method')]
        )
        print(f'Kruskal-Wallis test p-value: {p}')
        if p < 0.05:
            dunn = sp.posthoc_dunn(combined_df, val_col='Mutation_Count', group_col='Method', p_adjust='bonferroni')
            print(dunn)

    # Step 4: Chi-square Test for Mutated Positions
    positions_df = pd.crosstab(combined_df['Position'], combined_df['Method'])
    chi2, p, _, _ = stats.chi2_contingency(positions_df)
    print(f'Chi-square test p-value: {p}')











if __name__ == "__main__":
    main()