import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import scikit_posthocs as sp
import math



def find_first_prediction_changes(df):
    """
    Identify the first occurrence of each new `overall_prediction` category 
    (excluding the first one) within a replicate.
    """
    # Shift the `overall_prediction` column to compare with the previous row.
    previous_predictions = df['overall_prediction'].shift(1)

    # Identify where the prediction changes.
    changes = df[df['overall_prediction'] != previous_predictions]

    # Exclude the first category by skipping the first row.
    first_changes = changes.iloc[1:]

    # Keep only the first instance of each unique `overall_prediction` change.
    first_unique_changes = first_changes.drop_duplicates(subset=['overall_prediction'], keep='first')

    return first_unique_changes



def process_methods_and_replicates(methods):
    """
    Process all methods, each containing multiple replicates, and aggregate results.
    """
    aggregated_results = []

    # Iterate through each method and its replicates
    for method_name, replicates in methods.items():
        for replicate_name, df in replicates.items():
            # Extract first prediction changes for the replicate
            changes = find_first_prediction_changes(df)

            # Add method and replicate information to the dataframe
            changes['method'] = method_name

            # Append the results to the aggregated list
            aggregated_results.append(changes)

    # Combine all results into a single dataframe
    return pd.concat(aggregated_results, ignore_index=True)


def main():
    input_files = snakemake.input


    grouped_results = []

    for file in input_files:
        df = pd.read_csv(file)
        changes = find_first_prediction_changes(df)

        grouped_results.append(changes)
        
    grouped_df = pd.concat(grouped_results, ignore_index=True)

    grouped_df.to_csv('TESTFILE.csv')



    # Create a box plot for each overall_prediction category
    g = sns.catplot(
        data=grouped_df,
        x='method',
        y='num_mutation',
        col='overall_prediction',  # Create separate plots for each category
        kind='box',
        height=5,  # Adjust the height of each plot
        aspect=1.2  # Control the aspect ratio of each plot
    )

    # Adjust the title and labels
    g.figure.suptitle('Comparison of Number of Mutations Across Methods by Prediction Category', 
                y=1.05)  # Adjust the title position
    g.set_axis_labels('Method', 'Number of Mutations')

    plt.savefig(snakemake.output.boxplot)
    plt.close() # close to save memory









    # # method_name = snakemake.wildcards.method_name

    # replicates = [pd.read_csv(file) for file in input_files]
    # combined_df = pd.concat(replicates, ignore_index=True)


    # # make a boxplot of the stats
    # # plt.figure(figsize=(10, 6))
    # # sns.boxplot(x='method', y='num_mutation', data=combined_df)
    # # plt.title(f'Mutation Counts by {method_name}')
    # # plt.ylabel("Number of mutations")
    # # plt.savefig(snakemake.output.boxplot)
    # # plt.close() # close to save memory

    # # Get unique methods for plotting
    # methods = combined_df['method'].unique()

    # # Calculate the number of rows and columns for the grid layout
    # num_methods = len(methods)
    # cols = 3  # Define the number of columns in the grid
    # rows = math.ceil(num_methods / cols)  # Calculate the number of rows needed

    # fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows), squeeze=False)

    # # Generate boxplots for each method
    # for i, method in enumerate(methods):
    #     row, col = divmod(i, cols)  # Calculate grid position
    #     sns.boxplot(
    #         x='method', y='num_mutation', 
    #         data=combined_df[combined_df['method'] == method], ax=axes[row, col]
    #     )
    #     axes[row, col].set_title(f'Mutation Counts for {method}')

    # plt.ylabel("Number of mutations")
    # plt.savefig(snakemake.output.boxplot)
    # plt.close() # close to save memory



    # # 3.1 Shapiro-Wilk test for normality
    # for method, data in combined_df.groupby('method'):
    #     stat, p = stats.shapiro(data['num_mutation'])
    #     print(f'{method} - Shapiro-Wilk test p-value: {p}')

    # # 3.2 Levene’s test for equal variances
    # stat, p = stats.levene(
    #     *[group['num_mutation'].values for name, group in combined_df.groupby('method')]
    # )
    # print(f'Levene’s test p-value: {p}')

    # # 3.3 Hypothesis Testing (ANOVA or Kruskal-Wallis)
    # if p > 0.05:  # If variances are equal
    #     stat, p = stats.f_oneway(
    #         *[group['num_mutation'].values for name, group in combined_df.groupby('method')]
    #     )
    #     print(f'ANOVA test p-value: {p}')
    #     if p < 0.05:
    #         tukey = pairwise_tukeyhsd(combined_df['num_mutation'], combined_df['method'])
    #         print(tukey)
    # else:
    #     stat, p = stats.kruskal(
    #         *[group['num_mutation'].values for name, group in combined_df.groupby('method')]
    #     )
    #     print(f'Kruskal-Wallis test p-value: {p}')
    #     if p < 0.05:
    #         dunn = sp.posthoc_dunn(combined_df, val_col='num_mutation', group_col='method', p_adjust='bonferroni')
    #         print(dunn)

    # # Step 4: Chi-square Test for Mutated Positions
    # positions_df = pd.crosstab(combined_df['mutated_positions'], combined_df['method'])
    # chi2, p, _, _ = stats.chi2_contingency(positions_df)
    # print(f'Chi-square test p-value: {p}')









if __name__ == "__main__":
    main()