import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
from scipy.stats import shapiro, kruskal
import os
import glob



# Define helper functions
def get_initial_category(df):
    return df['overall_prediction'].iloc[0]

def find_first_prediction_changes(df):
    initial_category = get_initial_category(df)
    previous_predictions = df['overall_prediction'].shift(1)

    changes = df[df['overall_prediction'] != previous_predictions]

    valid_changes = changes[changes['overall_prediction'] != initial_category].iloc[1:]

    first_unique_changes = valid_changes.drop_duplicates(subset=['overall_prediction'], keep='first')
    return first_unique_changes


def find_overall_prediction_changes(df):
    initial_category = get_initial_category(df)
    previous_predictions = df['overall_prediction'].shift(1)
    changes = df[df['overall_prediction'] != previous_predictions]
    valid_changes = changes[changes['overall_prediction'] != initial_category].iloc[1:]

    filtered_changes = []
    for idx, row in valid_changes.iterrows():
        current_prediction = row['overall_prediction']
        if (idx + 4 < len(df) and 
            df.loc[idx + 1, 'overall_prediction'] == current_prediction and 
            df.loc[idx + 2, 'overall_prediction'] == current_prediction and
            df.loc[idx + 3, 'overall_prediction'] == current_prediction and
            df.loc[idx + 4, 'overall_prediction'] == current_prediction):
            filtered_changes.append(row)
    
    return pd.DataFrame(filtered_changes)


def clean_name(name):
    return name.replace('_', ' ')

# Load the data for all three datasets
def load_and_combine_data(file_paths):
    combined_data = []
    for path, dataset_name in file_paths.items():
        df = pd.read_csv(path)
        df['dataset'] = dataset_name  # Add a column to distinguish datasets
        combined_data.append(df)
    return pd.concat(combined_data, ignore_index=True)



# Function to create boxplots
def create_boxplots(data, title, output_path):
    g = sns.catplot(
        data=data,
        x='method',
        y='num_mutation',
        col='overall_prediction',
        kind='box',
        height=7,
        aspect=0.8
    )
    for ax in g.axes.flat:
        original_title = ax.get_title().split(' = ')[1]
        ax.set_title(clean_name(original_title), fontsize=12, fontweight='bold', pad=10)
        ax.set_xlabel('')
        for label in ax.get_xticklabels():
            label.set_rotation(30)
            label.set_ha('right')
    g.set_axis_labels('Method', 'Number of Mutations')
    g.figure.suptitle(title, fontsize='x-large', fontweight='bold')
    # g.set(ylim=(0, 150))
    g.figure.subplots_adjust(bottom=0.25, top=0.88)
    plt.savefig(output_path)
    plt.close()

def get_prediction_order(initial_category):
    base_order = ['NR1', 'NR1-like', 'NR4-like', 'NR4', 'other']

    if initial_category == 'NR1':
        return base_order[1:]
    elif initial_category == 'NR4':
        return ['NR4-like', 'NR1-like', 'NR1', 'other']
    else:
        # default to base
        return base_order
    
# Shapiro-Wilk Test
def run_shapiro_tests(df, output_path):
    categories = df['overall_prediction'].dropna().unique()
    methods = df['method'].dropna().unique()
    results = []

    for category in categories:
        for method in methods:
            method_data = df[(df['overall_prediction'] == category) & (df['method'] == method)]['num_mutation']
            if len(method_data) >= 3:  # At least 3 samples required for Shapiro-Wilk test
                stat, p_value = shapiro(method_data)
                results.append({'Category': category, 'Method': method, 'W-Statistic': stat, 'p-value': p_value})
            else:
                results.append({'Category': category, 'Method': method, 'W-Statistic': 'N/A', 'p-value': 'Insufficient data'})

    # Save results to a file
    with open(output_path, 'w') as f:
        f.write("Shapiro-Wilk Test Results\n")
        f.write("=" * 40 + "\n")
        for result in results:
            f.write(f"Category: {result['Category']}, Method: {result['Method']}, W-Statistic: {result['W-Statistic']}, p-value: {result['p-value']}\n")

# Kruskal-Wallis Test
def run_kruskal_wallis(df, output_path):
    categories = df['overall_prediction'].unique().dropna()
    with open(output_path, 'w') as f:
        f.write("Kruskal-Wallis Test Results\n")
        f.write("=" * 40 + "\n")
        for category in categories:
            category_data = df[df['overall_prediction'] == category]
            grouped_data = [group['num_mutation'].values for _, group in category_data.groupby('method') if len(group) > 0]
            if len(grouped_data) < 2:
                f.write(f"Not enough methods for comparison in category: {category}\n")
                continue
            stat, p_value = kruskal(*grouped_data)
            f.write(f"Category: {category}, Statistic: {stat:.4f}, p-value: {p_value:.4e}\n")

# QQ Plot
def plot_qq_grid(df, output_path):
    categories = df['overall_prediction'].dropna().unique()
    methods = df['method'].dropna().unique()
    num_categories = len(categories)
    num_methods = len(methods)
    fig, axes = plt.subplots(num_categories, num_methods, figsize=(6 * num_methods, 6 * num_categories), squeeze=False)
    fig.suptitle('Q-Q Plots for All Categories and Methods', fontsize=18, fontweight='bold')

    for i, category in enumerate(categories):
        for j, method in enumerate(methods):
            ax = axes[i, j]
            method_data = df[(df['overall_prediction'] == category) & (df['method'] == method)]['num_mutation']
            if method_data.empty:
                ax.axis('off')
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center', fontsize=12)
            else:
                stats.probplot(method_data, dist="norm", plot=ax)
                ax.set_title(f'{category} - {method}', fontsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(output_path)
    plt.close()


# # Load, aggregate, and combine data for all datasets and methods
# def load_and_aggregate_data(base_path, datasets, methods, replicates):
#     combined_data = []
#     for dataset in datasets:
#         for method in methods:
#             # Pattern to match all replicate files for each dataset and method
#             file_pattern = os.path.join(base_path, dataset, method, "results", "*.csv")
#             files = glob.glob(file_pattern)
#             aggregated_df = pd.concat([pd.read_csv(file) for file in files], ignore_index=True)
#             aggregated_df['dataset'] = dataset
#             aggregated_df['method'] = method
#             combined_data.append(aggregated_df)
#     return pd.concat(combined_data, ignore_index=True)



def main():
    # Define parameters
    datasets = ['cd70', 'cd80', 'cd85']
    methods = ['random', 'nonconservative', 'grantham_distances', 'marginal_weights', 'ConSurf']
    reps = range(1, 51)
    datafile_names = {
        "cd70": ["NR1toNR4_N6_N81", "NR4toNR1_N81_N6"],
        "cd80": ["NR1toNR4_N7_N186", "NR4toNR1_N186_N7"],
        "cd85": ["NR1toNR4_N7_N299", "NR4toNR1_N299_N7"]
    }

    # Dictionary to store the data
    dataframes = {}

    # Load only the specified files
    for dataset in datasets:
        for method in methods:
            for datafile_name in datafile_names[dataset]:
                for rep in reps:
                    # Construct file path
                    file_path = f"workflows/{dataset}/{method}/results/{datafile_name}_{rep}.csv"
                    
                    # Check if file exists before loading
                    if os.path.exists(file_path):
                        # Initialize storage for each unique combination if not already present
                        key = (dataset, datafile_name, method)
                        if key not in dataframes:
                            dataframes[key] = []
                        
                        # Read the CSV and append to the appropriate list
                        df = pd.read_csv(file_path)
                        df['dataset'] = dataset
                        df['datafile_name'] = datafile_name
                        df['method'] = method
                        dataframes[key].append(df)

    # Concatenate all dataframes for each key
    for key in dataframes:
        dataframes[key] = pd.concat(dataframes[key], ignore_index=True)


    # Initialize lists to store aggregated results
    first_changes_list = []
    overall_changes_list = []

    # Process each dataframe in dataframes dictionary
    for (dataset, datafile_name, method), df in dataframes.items():
        # Find first and overall prediction changes
        first_changes = find_first_prediction_changes(df)
        overall_changes = find_overall_prediction_changes(df)
        
        # Annotate and add to lists
        first_changes['method'] = method
        first_changes['dataset'] = dataset
        first_changes_list.append(first_changes[['dataset', 'method', 'overall_prediction', 'num_mutation']])

        overall_changes['method'] = method
        overall_changes['dataset'] = dataset
        overall_changes_list.append(overall_changes[['dataset', 'method', 'overall_prediction', 'num_mutation']])

    # Concatenate aggregated data
    first_changes_df = pd.concat(first_changes_list, ignore_index=True)
    overall_changes_df = pd.concat(overall_changes_list, ignore_index=True)


    print(f'first changes: \n {first_changes_df}')
    print(f'overall changes: \n {overall_changes_df}')


    # Generate boxplots
    create_boxplots(first_changes_df, "First Prediction Changes by Method - Combined Data", "results/first_changes_boxplot.png")
    create_boxplots(overall_changes_df, "Overall Prediction Changes by Method - Combined Data", "results/overall_changes_boxplot.png")




if __name__ == "__main__":
    main()

