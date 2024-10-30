import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import scikit_posthocs as sp
import math
from scipy.stats import kruskal
import statsmodels.api as sm
from statsmodels.formula.api import ols


def get_initial_category(df):
    return df['overall_prediction'].iloc[0]

def find_first_prediction_changes(df):
    """
    Return the first row where each new `overall_prediction` occurs,
    excluding the initial category and any reversion to it.
    """
    # Identify the initial category from the first row.
    initial_category = get_initial_category(df)
    
    # Shift the `overall_prediction` column to compare with the previous row.
    previous_predictions = df['overall_prediction'].shift(1)

    # Identify where the prediction changes.
    changes = df[df['overall_prediction'] != previous_predictions]

    # Exclude the first row and any rows where the category reverts to the initial one.
    valid_changes = changes[changes['overall_prediction'] != initial_category].iloc[1:]

    # Keep only the first instance of each unique `overall_prediction` change.
    first_unique_changes = valid_changes.drop_duplicates(subset=['overall_prediction'], keep='first')

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


def clean_name(name):
    return name.replace('_', ' ')


def get_prediction_order(initial_category):
    base_order = ['NR1', 'NR1-like', 'NR4-like', 'NR4', 'other']

    if initial_category == 'NR1':
        return base_order[1:]
    elif initial_category == 'NR4':
        base_order = base_order[:-1:-1]
        base_order.append('other')
        return base_order
    else:
        # default to base
        return base_order
    

    # Define a function to perform ANOVA or Kruskal-Wallis and post-hoc testing
def perform_statistical_tests(df, category):
    """Perform statistical tests for a given category."""
    # Filter data for the current category
    category_data = df[df['overall_prediction'] == category]

    # Check if there are enough unique methods for comparison
    if category_data['method'].nunique() < 2:
        print(f"Not enough methods for comparison in category: {category}")
        return None

    # Perform Kruskal-Wallis test (non-parametric)
    kruskal_result = kruskal(
        *[group['num_mutation'].values for name, group in category_data.groupby('method')]
    )
    print(f"Kruskal-Wallis result for {category}: p-value = {kruskal_result.pvalue}")

    if kruskal_result.pvalue < 0.05:
        # Perform post-hoc Tukey HSD test if significant
        tukey = pairwise_tukeyhsd(
            endog=category_data['num_mutation'],
            groups=category_data['method'],
            alpha=0.05
        )
        print(f"Tukey HSD post-hoc test for {category}:\n{tukey}\n")
    else:
        print(f"No significant differences found in {category}.\n")




def plot_qq_grid(df, outpath):
    """Generate a grid of Q-Q plots with categories as rows and methods as columns."""
    # Get unique categories and methods
    print(df)
    categories = df['overall_prediction'].dropna().unique()
    methods = df['method'].dropna().unique()
    print(f'methods: {methods}')
    print(f'categories: {categories}')

    # Dynamically adjust the grid size based on valid combinations
    num_categories = len(categories)
    num_methods = len(methods)

    if num_categories == 0 or num_methods == 0:
        print("No valid data to plot.")
        return

    # Create the grid: categories in rows, methods in columns
    fig, axes = plt.subplots(
        num_categories, num_methods, figsize=(6 * num_methods, 6 * num_categories),
        squeeze=False  # Ensure we always get a 2D array of axes
    )
    fig.suptitle('Q-Q Plots for All Categories and Methods', fontsize=18, fontweight='bold')

    # Iterate over categories and methods to populate the grid
    for i, category in enumerate(categories):
        for j, method in enumerate(methods):
            ax = axes[i, j]  # Select the correct axis

            # Filter data for the current category and method
            method_data = df[(df['overall_prediction'] == category) & (df['method'] == method)]['num_mutation']

            if method_data.empty:
                # If no data, disable the axis and display a message
                ax.axis('off')
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center', fontsize=12)
            else:
                # Generate the Q-Q plot
                stats.probplot(method_data, dist="norm", plot=ax)
                ax.set_title(f'{category} - {method}', fontsize=12)

    # Adjust layout to prevent overlap and save the plot
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(outpath)
    plt.close()


def main():
    input_files = snakemake.input


    grouped_results = []

    for file in input_files:
        df = pd.read_csv(file)
        changes = find_first_prediction_changes(df)

        grouped_results.append(changes)

    # can just run this on the last one as it'll be the same for all of them
    initial_category = get_initial_category(df)
    
    
    grouped_df = pd.concat(grouped_results, ignore_index=True)
    grouped_df.to_csv('TESTFILE.csv')


    order = get_prediction_order(initial_category)
    grouped_df['method'] = grouped_df['method'].apply(clean_name)
    grouped_df['overall_prediction'] = pd.Categorical(
        grouped_df['overall_prediction'], categories=order, ordered=True
    )


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

    # Remove individual x-axis labels from each subplot
    for ax in g.axes.flat:
        # Extract the original category name from the title
        original_title = ax.get_title().split(' = ')[1]  # Extract the prediction value
        # Set the cleaned title (with underscores replaced by spaces)
        ax.set_title(clean_name(original_title), fontsize=14, fontweight='bold', pad=10)
        ax.set_xlabel('')

    # Adjust the title and labels
    g.set_axis_labels('Method', 'Number of Mutations')
    # Add a single x-axis label for the entire plot
    plt.savefig(snakemake.output.boxplot)
    plt.close() # close to save memory




    # Get the unique categories from the data
    categories = grouped_df['overall_prediction'].unique().dropna()

    plot_qq_grid(grouped_df, snakemake.output.qqplot)

    # # Run statistical tests for each category
    # for category in categories:
    #     perform_statistical_tests(grouped_df, category)











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